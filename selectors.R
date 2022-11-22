library('tidyverse')
library('stringr')
library('DescTools')
library('dplyr')
library('jsonlite')


data_path <- function() { paste0(getwd(), '/datafiles/') }

columnToRealVector <- function(column) {
  r <- as.numeric(unlist(column))
  r[(!is.na(r) & is.finite(r))]
}

legacyVariantsAll <- function() {
  c('WildType',
    'D614G',
    'B.1.1.7',
    'B.1.351',
    'B.1.617.2',
    'B.1.1.529')
}

omiStudyVariants <- function() { c('B.1.1.7', 'B.1.617.2', 'B.1.1.529') }

variantFromGreekNames <- function(variants) {
  ans <- rep(NA, length(variants))

  for (k in seq_len(length(variants))) {
    variant <- variants[k]
    ans[k] <- switch(tolower(variant),
                     'delta' = 'B.1.617.2',
                     'alpha' = 'B.1.1.7',
                     'beta' = 'B.1.351',
                     'wild type' = 'WildType',
                     'd614g' = 'D614G',
                     'gamma' = 'P.1',
                     'omicron' = 'B.1.1.529',
                     'omicron ba.1' = 'B.1.1.529',
                     'omicron ba.2' = 'B.1.1.529',
                     variant
    )
  }
  ans
}

guessAge <- function(ageBand) {
  as.numeric(substr(ageBand, 1, 2)) + 2
}

vaccineName <- function(vaccine) {
  if (str_to_title(vaccine) != 'Comirnaty') {
    'Vaxzevria'
  } else { 'Comirnaty' }
}

readFile <- function(fileName) {
  file <- paste0(data_path(), fileName, '.csv')
  df <- as_tibble(read.csv(file)) %>%
    mutate(across(where(is.character), str_trim))
}

getEarlyVariantTitres <- function(vaccine = 'Comirnaty', noPriorCovid = F, exclude_wild_type = T) {

  vaccine <- vaccineName(vaccine)
  df <- readFile(paste0('Crick_', vaccine)) %>%
    filter(sampleOrderInVaccStatus == 1, COVID_vaccStatus %in% c(1, 2))

  ans <- tibble(
    Vaccine = vaccine,
    Dose = df$COVID_vaccStatus,
    DoseInterval = df$COVID_daysBetweenJabs,
    DaysSinceJab1 = df$COVID_daysSinceJab1,
    DaysSinceJab2 = df$COVID_daysSinceJab2,
    Site = df$site,
    Sex = df$sex,
    Age = guessAge(df$ageRange_by5),
    WildType = df$Wildtype_ic50,
    D614G = df$D614G_ic50,
    B.1.1.7 = df$B.1.1.7_ic50,
    B.1.351 = df$B.1.351_ic50,
    B.1.617.2 = df$B.1.617.2_ic50,
    Participant_ID = df$bc_participant_id,
    Prior_Covid = df$COVID_symptoms)

  if (exclude_wild_type) {
    ans$WildType <- NULL
  }
  if (noPriorCovid) {
    ans %>% filter(Prior_Covid == 0)
  } else { ans }

}

###### NB! Refactor Use left_join to exand/add a dataframe, instead of looping
combineEfficaciesTitres <- function(efficacies, titres) {
  df <- NULL

  for (k in seq_len(nrow(efficacies))) {
    row <- efficacies[k,]
    reduced <- titres[(titres$Dose == row$Dose),]
    current <- as.numeric(unlist(reduced[row$Variant]))
    currentGood <- !is.na(current)
    current <- current[(currentGood)]
    ages <- as.numeric(unlist(reduced$Age))[(currentGood)]
    sexes <- unlist(reduced$Sex)[(currentGood)]

    add <- tibble(
      Vaccine = row$Vaccine,
      Variant = row$Variant,
      Dose = row$Dose,
      Efficacy_Definition = row$Efficacy_Definition,
      Efficacy = row$Efficacy,
      Probability = row$Efficacy / 100,
      Efficacy_Age_Bands = row$Age,
      Country = row$Country,
      LCL = row$LCL,
      UCL = row$UCL,
      GMT = Gmean(current),
      Sex = list(sexes),
      Age = list(ages),
      Titres = list(current))
    df <- rbind(df, add)
  }

  df
}

# Comirnaty or Vaxzevria. Excluding the extreme efficacy of 10.4% for Vaxzevria against Beta in ZA; Mehdi et al 2021
getEfficacies <- function(vaccine = 'Comirnaty', isOnlyUK = F) {
  df <- readFile(paste0('efficacies_', vaccineName(vaccine))) %>%
    filter(Age == 'all', Efficacy_Definition == 'SYM')

  if (isOnlyUK) {
    df %>% filter(Country %in% c('England', 'England_August', 'UK', 'Scotland', 'AZ Trial'))
  } else { df }
}

titresToCensoredDF <- function(vector, low = 40, up = 2560) {
  vec <- vector[(!is.na(vector))]
  df <- tibble(
    left = rep(NA, length(vec)),
    right = rep(NA, length(vec))
  )

  for (k in seq_len(nrow(df))) {
    val <- vec[k]
    if (val < low) {
      df$right[k] <- low
    }
    else if (val > up) {
      df$left[k] <- up
    }
    else {
      df$left[k] <- val
      df$right[k] <- val
    }
  }
  df$left <- log10(unlist(df$left))
  df$right <- log10(unlist(df$right))

  as.data.frame(df)
}

#'Diabetes', 'Kidney' or 'Dialysis' or 'Disease' means not healthy.
comirnatyClinicalSubgroupTitres <- function(disease,
                                            excludeKnownCovid = T, exclude_wild_type = T) {
  disease <- str_to_title(disease)

  if (disease == 'Diabetes')
  { diabetic <- 'Y' }
  else if (disease == 'Kidney')
  { diabetic <- 'N' }
  else { diabetic <- c('Y', 'N') }

  df <- readFile('Crick_HD') %>%
    filter(vaccine == 'BNT162b2', diabetic %in% diabetic)

  # function adapted from https://github.com/EdjCarr/Crick-HD-AZD-BNT-VOCs-2021-07/
  getSeronegatives <- function() {
    base <- filter(df, time_point == "baseline_bleed" & (outcome == "NOT_DETECTED" &
      (!is.na(wildtype_ic50)) &
      wildtype_ic50 == 5 &
      D614G_ic50 == 5 &
      (pre_vaccine_positive_test == FALSE | is.na(pre_vaccine_positive_test))))

    base$research_identifier
  }

  if (excludeKnownCovid) {
    df <- filter(df, research_identifier %in% getSeronegatives())
  }
  df <- df %>% filter(time_point !=
                        'baseline_bleed')

  guessAge <- function(ageBand) {
    first <- substr(ageBand, 1, 1)
    if (first == '>') { 75 } else { 42 }
  }

  guessDose <- function(time_point) {
    as.numeric(substr(time_point, 4, 4))
  }

  ans <- tibble(
    Vaccine = 'Comirnaty',
    Dose = guessDose(df$time_point),
    Site = df$dialysis_centre_code,
    Elisa = df$outcome,
    Gender = df$gender,
    Age = sapply(df$age, guessAge),
    WildType = df$wildtype_ic50,
    D614G = df$D614G_ic50,
    B.1.1.7 = df$alpha_ic50,
    B.1.351 = df$beta_ic50,
    B.1.617.2 = df$delta_ic50,
    Participant_ID = df$research_identifier,
    Prior_Covid = df$pre_vaccine_positive_test,
    Diabetic = unname(sapply(df$diabetic, (\(x) x == 'Y')))
  )

  if (exclude_wild_type) {
    ans$WildType <- NULL
  }
  ans
}

# Diabetes, Kidney, Dialysis  or Healthy. NB! Estimates exclude prior Covids anyway
comirnatyClinicalSubgroupEfficacies <- function(disease, exclude_healthy_Whittaker = T) {

  df <- readFile('Subgroup_VE') %>% filter(Vaccine == 'Comirnaty', Disease == disease, Efficacy_Definition == 'SYM')
  # Negative VE, huge CI. We exclude it in estimation. But add it in a VE data table
  if (exclude_healthy_Whittaker)
  { filter(df, !(Ref == 'Whitaker2022' & Disease == 'Healthy')) }
  else { df }
}


toPredictionTitres <- function(titres) {
  variants <- intersect(colnames(titres), legacyVariantsAll())
  doses <- sort(unique(titres$Dose))

  ans <- NULL
  for (variant in variants) {
    for (dose in doses) {
      reduced <- filter(titres, Dose == dose)
      current <- as.numeric(unlist(reduced[, (variant)]))
      current <- current[(!is.na(current))]

      add <- tibble(
        Variant = variant,
        Dose = dose,
        GMT = Gmean(current),
        Titres = list(current)

      )
      ans <- rbind(ans, add)
    }
  }

  ans
}

getCovariants <- function(country) {
  countries <- c('USA', 'UK', 'Israel', 'Canada')
  ind <- match(country, countries)
  data <- fromJSON(paste0(data_path(), 'CovariantCounts.json'))
  dist <- data$regions$distributions
  df <- as.data.frame(dist) %>% filter(country %in% countries)
  df$distribution[[ind]]
}

#Source: John Hopkins https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
dailyCovidsPerCountry <- function(country) {

  nameToDate <- function(x) {
    str <- sub('.', '', x)
    as.Date(str, '%m.%d.%y')
  }

  covids <- readFile('global_covids') %>%
    filter(Country == country) %>%
    select(-c(Province, Country))
  cases <- rep(0, ncol(covids))
  for (i in seq_len(nrow(covids))) {
    cases <- cases + unlist(covids[i,])
  }
  shift <- c(0, unlist(cases[seq_along(cases) - 1]))
  daily <- cases - shift
  names(daily) <- NULL
  tibble(Cases = daily,
         Date = nameToDate(colnames(covids)))
}

comirnatyAgeEfficaciesAndrews <- function() {
  df <- readFile('Andrews_Age')
  age_group <- NULL
  for (k in seq_len(nrow(df))) {
    row <- df[k,]
    if (row$Age_End == 80) {
      group <- 'ALL'
    } else if (row$Age_End == 64) {
      group <- 'OLDER'
    } else {
      group <- 'YOUNGER'
    }
    age_group <- c(age_group, group)
  }
  df$Age_Group <- age_group

  filter(df, Age_End != 80)
}


################################ OMICRON

# subgroup: 1 two-doses, 2-6 weeks; 2: two-doses, 12-
omicronTitres <- function(vaccine = 'Comirnaty', noPriorCovid = F) {
  vac_code <- if (vaccine == 'Comirnaty') {
    'BNT162b2'
  }else { 'AZD1222' }
  df <- readFile('Crick_Omicron') %>%
    filter(dose_1 == vac_code, dose_2 == vac_code) %>%
    filter(cohort %in% c('POST-Boost', '2-6', '12-16')) %>%
    select(-c(dose_1,
              dose_2, X)) %>%
    rename_with(~c('B.1.1.7', 'B.1.617.2', 'B.1.1.529', 'Omi_ID', 'Cohort'),
                c('ic50_Alpha',
                  'ic50_Delta',
                  'ic50_Omicron', 'elig_study_id', 'cohort')) %>%                               # Replacing values
    mutate(Cohort = recode(Cohort, 'POST-Boost' = 'Boost', '2-6' = 'Early', '12-16' =
      'Late'))
  df$Vaccine <- vaccine

  if (noPriorCovid) {
    df %>% filter(symptomatic_beforeVisit == 'No')
  }
  else { df }
}

groupedOmicronTitres <- function(titres, study_divisor = 1, cohorts = c('Early',
                                                                        'Late', 'Boost')) {
  df <- tibble(
    Vaccine = character(),
    Cohort = character(),
    Variant = character(),
    Dose = integer(),
    Study_Divisor = numeric(),
    GMT = numeric(),
    Sex = list(),
    Age = list(),
    Omi_ID = list(),
    Titres = list()
  )
  #Todo: What if other cohorts?
  variants <- omiStudyVariants()
  for (cohort in cohorts) {
    dose <- if (cohort == 'Boost') { 3 } else { 2 }
    current_rows <- titres %>% filter(Cohort == cohort)
    ids <- current_rows$Omi_ID
    for (variant in variants) {
      row_titres <- current_rows[(variant)]
      good_ids <- !is.na(row_titres) #good titres
      variant_titres <- columnToRealVector(row_titres) / study_divisor
      gmt <- Gmean(variant_titres) # still censored
      ages <- as.numeric(unlist(guessAge(current_rows$age[(good_ids)])))

      df <- add_row(df,
                    Vaccine = titres$Vaccine[1],
                    Cohort = cohort,
                    Variant = variant,
                    Dose = dose,
                    Study_Divisor = study_divisor,
                    GMT = gmt,
                    Sex = list(unlist(current_rows$sex[(good_ids)])),
                    Age = list(ages),
                    Omi_ID = list(ids[(good_ids)]),
                    Titres = list(variant_titres)
      )
    }
  }
  df
}

combineOmicronEfficacyWithPredictions <- function(eff, predictions) {
  add <- NULL

  for (k in seq_len(nrow(eff))) {
    row <- eff[k,]
    predict <- filter(predictions, Variant == row$Variant, Cohort == row$Cohort, Dose == row$Dose)

    current_eff <- tibble(Vaccine = row$Vaccine,
                          Weeks = row$Weeks,
                          Efficacy = row$Efficacy,
                          Efficacy_Mean = row$Efficacy_Mean,
                          Prediction_Error = predict$Mean - row$Efficacy_Mean,
                          Is_Singleton = row$Is_Singleton,
                          LCL = row$LCL,
                          UCL = row$UCL)
    add <- rbind(add, cbind(current_eff, predict))
  }
  add

}

#Assumes time for two periods max
#normalise_by_uncertainty weighs Eff estimate for different var over intervals.
#  CI length; CI_Squared squared length
comirnatyOmicronDeltaEfficacies <- function(normalise_by_uncertainty = NULL) {
  #NB: now discards England over 50 VE, AndrewsStowe and uses 18-49 yo. But identical VE
  df <- readFile('Omicron_Delta_Efficacy') %>% filter(Cohort != "")
  studies <- unique(df$Study_ID)
  ans <- NULL

  for (id in studies) {
    study <- filter(df, Study_ID == id)
    for (variant in unique(study$Variant)) {
      for (dose in unique(study$Dose)) {
        now <- filter(study, Variant == variant, Dose == dose)
        for (cohort in unique(now$Cohort)) {
          weeks <- switch(cohort,
                          'Early' = c(2, 6),
                          'Late' = c(12, 16),
                          'Boost' = c(2, 6)
          )

          eff <- efficacyOverOneOrTwoPeriods(filter(now, Cohort == cohort), weeks,
                                             normalise_by_uncertainty = normalise_by_uncertainty)
          eff$Study <- id
          eff$Vaccine <- 'Comirnaty'
          eff$Cohort <- cohort
          eff$Weeks <- paste(weeks, collapse = '-')
          eff$Variant <- variantFromGreekNames(variant)
          eff$Dose <- dose

          ans <- rbind(ans, eff)
        }
      }
    }
  }

  ans
}

#normalise_by_uncertainty weighs Eff estimate for different var over intervals.
#  CI length; CI_Squared squared length
efficacyOverOneOrTwoPeriods <- function(rows, weeks, normalise_by_uncertainty = NULL) {

  firstRow <- function() {
    row <- rows[1,]
    tibble(Efficacy = row$Efficacy,
           LCL = row$LCL,
           UCL = row$UCL)
  }

  #accross two rows, needs combining wrt weeks. Assumes
  end_1 <- rows$End[1]

  if (nrow(rows) == 1 | end_1 >= weeks[2]) {
    return(firstRow())
  }
  count_1 <- end_1 - weeks[1] + 1
  count_2 <- weeks[2] - rows$Start[2] + 1
  w_1 <- count_1 / (count_1 + count_2)
  w_2 <- count_2 / (count_1 + count_2)

  ci_1 <- rows$UCL[1] - rows$LCL[1]
  ci_2 <- rows$UCL[2] - rows$LCL[2]
  if (!is.null(normalise_by_uncertainty))
  {
    ########## CI normalisation,
    u_1 <- 1 / ci_1
    u_2 <- 1 / ci_2
    if (normalise_by_uncertainty == 'CI_Squared') {
      u_1 <- u_1^2
      u_2 <- u_2^2
    }
    denom <- w_1 * u_1 + w_2 * u_2
    w_1 <- w_1 * u_1 / denom
    w_2 <- w_2 * u_2 / denom
  }

  eff <- w_1 * rows$Efficacy[1] + w_2 * rows$Efficacy[2]
  index <- which.min(c(ci_1, ci_2))

  left <- rows$Efficacy[index] - rows$LCL[index]
  right <- rows$UCL[index] - rows$Efficacy[index]

  tibble(Efficacy = eff,
         LCL = eff - left,
         UCL = eff + right)
}

omicronEstimationTable <- function(eff, grouped_omi, is_delta = F) {
  variant <- if (is_delta) { 'B.1.617.2' }else { 'B.1.1.529' }
  titres <- filter(grouped_omi, Variant == variant)
  eff <- filter(eff, Variant == variant)

  df <- NULL
  for (k in seq_len(nrow(eff))) {
    row <- eff[k,]
    reduced <- titres[(titres$Cohort == row$Cohort),]

    now <- tibble(Vaccine = 'Comirnaty',
                  Variant = variant,
                  Cohort = row$Cohort,
                  Dose = row$Dose,
                  Efficacy = row$Efficacy,
                  Probability = row$Efficacy / 100,
                  LCL = row$LCL,
                  UCL = row$UCL,
                  GMT = reduced$GMT,
                  Age = reduced$Age,
                  Titres = reduced$Titres,
                  Omi_ID = reduced$Omi_ID
    )
    df <- rbind(df, now)
  }

  df
}

omicronDialysisTitres <- function(vaccine = 'Comirnaty') {

  hd <- readFile('omi_HD')
  if (vaccine == 'Comirnaty') {
    df <- filter(hd, vaccine_1 == 'BNT162b2', vaccine_2 == 'BNT162b2')
  }else {
    df <- filter(hd, vaccine_1 == 'AZD1222', vaccine_2 == 'AZD1222')
  }
  # booster dose only
  boost <- filter(df, !is.na(daysSinceDose3) & daysSinceDose3 > 0, time_point == 'vax3_28d_bleed')

  boost$Vaccine <- 'Comirnaty'
  boost$Dose <- 3
  boost$Cohort <- 'Boost'
  boost <- boost %>% rename_with(~c('B.1.617.2', 'B.1.1.529', 'Omi_ID'),
                                 c('ic50_Delta',
                                   'ic50_Omicron', 'research_identifier'))
  boost
}


groupedOmicronDialysisTitres <- function(titres) {
  df <- tibble(
    Variant = character(),
    GMT = numeric(),
    Omi_ID = list(),
    Titres = list()
  )
  #NB: What if other cohorts?

  variants <- c('B.1.617.2', 'B.1.1.529')
  ids <- titres$Omi_ID
  for (variant in variants) {
    row_titres <- titres[(variant)]

    good_ids <- !is.na(row_titres) #good titres
    variant_titres <- columnToRealVector(row_titres)
    gmt <- Gmean(variant_titres) # still censored

    df <- add_row(df,
                  Variant = variant,
                  GMT = gmt,
                  Omi_ID = list(ids[(good_ids)]),
                  Titres = list(variant_titres)
    )
  }
  df$Vaccine <- titres$Vaccine[1]
  df$Dose <- titres$Dose[1]
  df$Cohort <- titres$Cohort[1]

  df
}

########## From Cromer2021 	https://doi.org/10.1016/s2666-5247(21)00267-6
getPredictionsCromer <- function() {
  variants <- c('B.1.1.7', 'B.1.351', 'B.1.617.2')

  com <- tibble(
    Vaccine = 'Comirnaty',
    Variant = variants,
    Mean = c(87.7, 57.3, 74.2),
    Mean_Low = c(73.4, 34.1, 53.6),
    Mean_Up = c(94.8, 76.1, 87.3))

  vz <- tibble(
    Vaccine = 'Vaxzevria',
    Variant = variants,
    Mean = c(62, 23.4, 41.1),
    Mean_Low = c(39.5, 9.5, 19.1),
    Mean_Up = c(79.1, 44.8, 62.5))
### NB: rename "Name" to inference type or source. As not clear.
  both <- rbind(vz, com)
  both$Dose <- 2
  both$Name <- 'CROMER'

  both
}