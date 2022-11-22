library('tidyverse')
library('stringr')
library('xtable')
library('dplyr')
library('countrycode')
library('functional')

source('selectors.R')

countryAbbrevs <- function(countries) {
  ans <- rep(NA, length(countries))

  for (k in seq_len(length(countries))) {
    ans[k] <- switch(countries[k],
                     'Scotland' = 'SCO',
                     'England' = 'ENG',
                     'England_August' = 'GBR',
                     'Israel' = 'ISR',
                     'US Trial' = 'USA',
                     'Canada' = 'CAN',
                     'France' = 'FRA',
                     'AZ Trial' = 'GBR',
                     'Qatar' = 'QAT',
                     'Sweden' = 'SWE',
                     'Spain' = 'ESP',
                     'US' = 'USA',
                     'UK' = 'GBR',
                     'SouthAfrica Trial' = 'ZAF',
                     countries[k]
    )
  }
  ans

}

greekVariantNames <- function(variants, is_title = F) {
  ans <- rep(NA, length(variants))

  for (k in seq_len(length(variants))) {
    current <- switch(variants[k],
                      'B.1.617.2' = 'DELTA',
                      'B.1.1.7' = 'ALPHA',
                      'B.1.351' = 'BETA',
                      'WildType' = 'WILD TYPE',
                      'D614G' = 'D614G',
                      'P.1' = 'GAMMA',
                      'B.1.1.529' = 'OMICRON',
                      variants[k]
    )
    if (is_title & current != 'D614G') {
      current <- str_to_title(current)
    }
    ans[k] <- current
  }
  ans
}

orderedVariants <- function() {
  c('D614G', 'B.1.1.7', 'B.1.351', 'B.1.617.2')
}

#use_brackets NULL = blank; True = Square, False = Usual
confidenceStrings <- function(lows, ups, estimates = NULL, separator = ', ',
                              use_brackets = F,
                              decimals = 1) {
  trailer <- paste0('%.', decimals, 'f')

  ans <- rep(NA, length(lows))
  for (k in seq_along(lows)) {
    a <- sprintf(trailer, round(lows[k], decimals))
    b <- sprintf(trailer, round(ups[k], decimals))
    if (is.null(use_brackets)) {
      interv <- paste0(a, separator, b) }
    else if (use_brackets) {
      interv <- paste0('[', a, separator, b, ']') }
    else {
      interv <- paste0('(', a, separator, b, ')') }

    if (!is.null(estimates)) {
      est <- sprintf(trailer, round(estimates[k], decimals))
      interv <- paste0(est, ' ', interv)
    }
    ans[k] <- interv
  }
  ans
}


latexCI <- function(lows, ups, estimates, percent_est = F, separator = ', ',
                    square_brackets = F,
                    decimals = 1, ci_bold = T) {

  trailer <- paste0('%.', decimals, 'f')

  local <- function(low, up, est) {
    a <- sprintf(trailer, round(low, decimals))
    b <- sprintf(trailer, round(up, decimals))
    est <- sprintf(trailer, round(est, decimals))

    interval <-
      if (square_brackets) {
        paste0('\\,[', a, separator, b, ']$') }
      else {
        paste0('\\,(', a, separator, b, ')$') }
    if (percent_est) {
      est <- paste0(est, '\\%')
    }
    if (ci_bold) {
      paste0('$', '\\mathbf{', est, '}', interval) }
    else {
      paste0('$', est, interval) }
  }

  ans <- rep(NA, length(lows))
  for (k in seq_along(lows)) {
    ans[k] <- local(lows[k], ups[k], estimates[k])
  }
  ans
}

figureCI <- function(est, low, up, decimals = 1, is_bold = T, start = '') {
  trail <- if (is_bold) { '<b>' } else { '' }
  ans <- paste0(trail, start, round(est, decimals), ' (95% CI: ', round(low, decimals),
                '--',
                round(up,
                      decimals), ')', trail)
}

efficacyTable <- function(df, covid_col = F) {
  #Study, Dose, Variant, Country, VE, Confidence, Ref, Comment?

  df <- df[(order(df$Study)),]

  covidReplacer <- function(x) {
    if (is.na(x)) {
      return('')
    }
    if (x == 0) { 'N' }else { 'Y' }
  }

  processStudy <- function(study) {
    ord <- order(match(study$Variant, orderedVariants()), study$Dose)
    study <- study[ord,]
    study$Ref[duplicated(study$Ref)] <- NA
    study$Study[duplicated(study$Study)] <- NA
    country_names <- countrycode(countryAbbrevs(study$Country), origin = 'iso3c',
                                 destination = 'country.name', custom_match =
                                   c('SCO' = 'Scotland', 'ENG' = 'England'))
    country_names[duplicated(country_names)] <- NA
    vocs <- greekVariantNames(study$Variant, is_title = T)
    ci <- confidenceStrings(study$LCL, study$UCL, use_brackets = T, separator = '--')
    prior_covid <- sapply(study$Prior_Covid, covidReplacer)

    tibble(Study = study$Study,
           Country = country_names,
           Dose = study$Dose,
           Covid = prior_covid,
           Variant = vocs,
           VE = sprintf('%.1f', study$Efficacy),
           CI = ci,
           Reference = study$Ref
    )
  }

  table <- NULL
  for (i in  unique(df$Study)) {
    current <- processStudy(filter(df, Study == i))
    ref_ids <- which(!is.na(current$Reference))
    current$Reference[ref_ids] <- paste0(" \\cite{", current$Reference[ref_ids], "}")
    table <- rbind(table, current)
  }
  if (!covid_col) {
    table$Covid <- NULL

  }
  starts <- list(unlist(which(!is.na(table$Study)) - 1))[[1]]
  starts <- starts[-1]
  titles <- rep("[0.2cm]", length(starts))
  addtorow <- list()
  addtorow$pos <- as.list(starts)
  addtorow$command <- as.vector(titles)

  print(xtable(table, digits = 0), include.rownames = FALSE, booktabs = TRUE,
        add.to.row = addtorow, sanitize.text.function = function(x) { x },
        tabular.environment = "longtable")
}


tableAgeEfficacies <- function() {
  #columns:  Ages, Dose, Variant, VE, CI
  df <- comirnatyAgeEfficaciesAndrews()
  df <- df[(order(df$Variant)),]

  ages <- NULL
  for (k in seq_len(nrow(df))) {
    row <- df[k,]
    ages <- c(ages, paste0(row$Age_Start, '--', row$Age_End))

  }

  ci <- confidenceStrings(lows = df$LCL, ups = df$UCL,
                          separator = '--', use_brackets = T)
  ans <- tibble(
    Age = ages,
    Dose = df$Dose,
    Variant = greekVariantNames(df$Variant, is_title = T),
    VE = sprintf('%.1f', df$Efficacy),
    CI = ci
  )

  print(xtable(ans, digits = 0), booktabs = T, sanitize.text.function = function(x) { x },
        include.rownames = FALSE)
}

tableDiseaseSubgroups <- function(df) {
  refs <- unique(df$Ref)

  diseaseStudy <- function(ref) {
    study <- filter(df, Ref == ref)
    ci <- confidenceStrings(lows = study$LCL, ups = study$UCL,
                            separator = '--', use_brackets = T)
    nas <- rep(NA, nrow(study) - 1)

    process <- function(share) {
      s <- sprintf('%.1f', round(100 * share[1], 1))
      c(s, nas)

    }

    study_ind <- match(ref, refs)
    cite <- paste0("\\cite{", ref, "}")
    tibble(
      Study = c(study_ind, nas),
      Country = c(study$Country[1], nas),
      Disease = study$Disease,
      VE = sprintf('%.1f', study$Efficacy),
      CI = ci,
      D614G = process(study$D614G),
      Alpha = process(study$B.1.1.7),
      Beta = process(study$B.1.351),
      Delta = process(study$B.1.617.2),
      Reference = c(cite, nas)
    )
  }

  studies <- NULL
  for (ref in refs)
  { studies <- rbind(studies, diseaseStudy(ref)) }

  starts <- match(refs, df$Ref) - 1
  titles <- rep("[0.2cm]", length(refs))
  addtorow <- list()
  addtorow$pos <- as.list(starts)
  addtorow$command <- as.vector(titles)

  print(xtable(studies, digits = 0), booktabs = T,
        add.to.row = addtorow,
        sanitize.text.function = function(x) { x },
        include.rownames = FALSE)
}

tableFullPrimaryModels <- function(df) {
  xtable(df, digits = 0, caption = 'Primary Fitted Models',
         label = 'T:primaryResults') %>% printTabular()
  NULL
}

################## Omicron
matchedGroupsTable <- function(rows) {
  iqr <- confidenceStrings(rows$Quart_1, rows$Quart_3, use_brackets = T, separator
    = ', ', decimals = 0)
  result <- tibble(
    Variant = rows$Variant,
    Dose = rows$Dose,
    Cohort = rows$Cohort,
    Weeks = rows$Weeks,
    N = rows$N,
    Median = rows$Median,
    IQR = iqr,
    Mean = rows$Mean,
    GMT = rows$GMT
  )
  print(xtable(result, digits = 0), include.rownames = FALSE)
}


tableEfficacyDeltaOmicron <- function(isBold = F) {
  df <- readFile('Omicron_Delta_Efficacy') %>% filter(Cohort != "")

  getAges <- function(start, end, replace49 = F) {
    if (is.na(start)) {
      return('adults')
    }
    if (replace49 & !is.na(end) & end == 49) {
      end <- NA
    }
    if (is.na(end)) {
      paste0('over ', start)
    } else {
      paste0(start, '--', end)
    }
  }

  squeezeRawStudy <- function(studies, study_id) {
    study <- filter(studies, Study == study_id)
    row <- study[1,]
    nas <- rep(NA, nrow(study) - 1)
    study$Study <- c(study_id, nas)
    study$Country <- c(row$Country, nas)
    if (study_id != 9) {
      study$Age <- c(row$Age, nas)
    } else if (nrow(study) == 4) {
      study$Age <- c(study$Age[1], NA, study$Age[3], NA)
    }
    study$Ref <- c(paste0("\\cite{", row$Ref, "}"), nas)
    study
  }

  spacer <- function(tab) {
    uniqs <- unique(filter(tab, !is.na(Study))$Study)
    starts <- match(uniqs, tab$Study) - 1
    titles <- rep("[0.2cm]", length(uniqs))
    addtorow <- list()
    addtorow$pos <- as.list(starts)
    addtorow$command <- as.vector(titles)
    addtorow
  }

  procPeriods <- function(rows) {
    first <- rows[1,]

    wk_1 <- paste(first$Start, first$End, sep = '--')
    ve_1 <- sprintf('%.1f', round(first$Efficacy, 1))
    if (nrow(rows) < 2) {
      wk_2 <- NA
      ve_2 <- NA

    } else {
      second <- rows[2,]
      wk_2 <- paste(second$Start, second$End, sep = '--')
      ve_2 <- sprintf('%.1f', round(second$Efficacy, 1))
    }

    tibble(WeeksOne = wk_1,
           VeOne = ve_1,
           WeeksTwo = wk_2,
           VeTwo = ve_2)

  }

  usedEffs <- function() {
    sum <- NULL

    for (id in unique(df$Study_ID)) {
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

            rows <- filter(now, Cohort == cohort)
            periods <- procPeriods(rows)
            eff <- efficacyOverOneOrTwoPeriods(rows, weeks)
            row <- rows[1,]

            ve <- sprintf('%.1f', round(eff$Efficacy, 1))
            if (isBold) {
              ve <- paste0('\\textbf{', ve, '}')
            }

            current <- tibble(Study = id,
                              Country = row$Country,
                              Age = getAges(row$Age_start, row$Age_end, replace49 = T),
                              Dose = dose,
                              Variant = row$Variant,
                              WeeksOne = periods$WeeksOne,
                              WeeksTwo = periods$WeeksTwo,
                              Weeks = paste(weeks[1], weeks[2], sep = '--'),
                              VeOne = periods$VeOne,
                              VeTwo = periods$VeTwo,
                              VE = ve,
                              Ref = row$Ref

            )
            sum <- rbind(sum, current)
          }
        }
      }
    }

    ans <- NULL
    for (k in unique(sum$Study)) {
      ans <- rbind(ans, squeezeRawStudy(sum, k))
    }
    ans
  }

  use_tab <- usedEffs()
  tab_sum <- xtable(use_tab, digits = 0, caption = 'Omicron and Delta VE estimates: Used', label = 'T:usedOmicronDeltaVE')
  print(tab_sum, booktabs = T, add.to.row = spacer(use_tab), sanitize.text.function = function(x) { x }, include.rownames = FALSE)
}

tableDeltaOmicronDifference <- function(df) {
  formatter0 <- Curry(confidenceStrings, separator = '--', use_brackets = T, decimals = 0)
  formatter <- Curry(confidenceStrings, separator = '--', use_brackets = T)
  omiVE <- sprintf('%.0f', round(df$Omi_Mean, 0))
  omiCI <- mapply(formatter0, df$Omi_Low, df$Omi_Up)

  deltaVE <- sprintf('%.0f', round(df$Delta_Mean, 0))
  deltaCI <- mapply(formatter0, df$Delta_Low, df$Delta_Up)

  difVE <- sprintf('%.1f', round(df$Dif_Mean, 1))
  difCI <- mapply(formatter, df$Dif_Low, df$Dif_Up)
  gmtCI <- mapply(formatter0, df$GMT_Low, df$GMT_Up)

  ans <- tibble(Cohort = df$Cohort,
                GMT = df$GMT,
                GmtCI = gmtCI,
                DeltaVE = deltaVE,
                DeltaCI = deltaCI,
                OmiVE = omiVE,
                OmiCI = omiCI,
                DifVE = difVE,
                DifCI = difCI)

  tab <- xtable(ans, digits = 0, caption = 'Difference between  Delta and Omicron', label = 'T:differenceDeltaOmi')
  print(tab, booktabs = T, include.rownames = FALSE)
}


tableDeltaOmicronModels <- function(df) {
  xtable(df, digits = 0, caption = 'Separate Comirnaty models for Delta and Omicron',
         label = 'T:separateDeltaOmicron') %>% printTabular()
}

meanSe <- function(vec, decimals_mean = 1, decimals_se = 1) {
  mn <- as.numeric(round(vec[[1]], decimals_mean))
  se <- as.numeric(round(vec[[2]], decimals_se))

  mns <- sprintf(paste0('%.', decimals_mean, 'f'), mn)
  if (mn < 0) {
    mns <- sub('.', '', mns)
    mns <- paste0('\\textminus ', mns)
  }
  se <- sprintf(paste0('%.', decimals_se, 'f'), se)
  paste0(mns, ' (', se, ')')
}


countPrimaryNabts <- function(df, vaccine = 'Comirnaty') {

  one_N <- countNabtsPerDose(df, vaccine, 1)
  two_N <- countNabtsPerDose(df, vaccine, 2)
  paste0(one_N, '/', two_N)
}

countNabtsPerDose <- function(df, vaccine = 'Comirnaty', dose) {

  df <- filter(df, Vaccine == vaccineName(vaccine), Dose == dose)
  if (nrow(df) == 0) {
    return(NA)
  }
  lens <- sapply(df$Titres, length)
  round(mean(lens), 0)
}

printTabular <- function(tab) {
  print('.............')
  print(tab, booktabs = T, sanitize.text.function = function(x) { x }, include.rownames = FALSE)
}

tablesVariantPredictions <- function(com, vax) {
  greek_variants <- c('D614G', 'Alpha', 'Beta', 'Delta')

  makeFrame <- function(df) {
    ans <- NULL
    for (i in seq_len(nrow(df))) {
      row <- df[i,]
      theta <- row$Theta[[1]]
      se <- row$SE[[1]]
      current <- tibble(Epoch = row$Epoch,
                        D614G = '',
                        Alpha = '',
                        Beta = '',
                        Delta = '',
                        alpha = meanSe(c(theta[1], se[1])),
                        beta = meanSe(c(theta[2], se[2])),
                        phi = meanSe(c(theta[3], se[3])),
                        NabtN = row$TitreCounts,
                        SampleVE = row$EffCounts
      )

      for (name in greek_variants) {
        variant <- variantFromGreekNames(name)
        if (variant %in% row$Base[[1]]) {
          current[(name)] <- '\\checkmark'
        } else if (row$Variant == variant | row$Epoch < 4) { current[(name)] <- '\\bullseye' }
      }
      ans <- rbind(ans, current)

    }
    ans <- distinct(ans, Epoch, alpha, beta, phi, .keep_all = T)
    ans <- ans[(order(ans$Epoch)),]


    ans
  }

  frame_com <- com %>%
    arrange(factor(Variant, levels = greek_variants)) %>%
    makeFrame()
  xtable(frame_com, digits = 0, caption = 'Comirnaty: Predictions per variant', label = 'T:fitsVariantComirnaty') %>%
    printTabular()

  frame_vax <- makeFrame(vax)
  frame_vax <- frame_vax[(order(frame_vax$Epoch)),]
  frame_vax$Epoch <- c('Direct', 'Direct', 'Indirect', 'Indirect')
  xtable(frame_vax, digits = 0, caption = 'Vaxzevria: Predictions per variant', label = 'fitsVariantVaxzevria') %>%
    printTabular()
}