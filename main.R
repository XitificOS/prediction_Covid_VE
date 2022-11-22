library('tidyverse')
library('dplyr')
library("plyr")
library('DescTools')
library('functional')

source('fitters.R')
source('selectors.R')
source('plotters.R')
source('presenters.R')
source('subgroups.R')



############# Helpers
#Computes Predicted VE and its bands
expectedMeanForTitres <- function(df, theta, orderGMT = T, ci = 0.95) {

  prob_df <- NULL
  for (k in seq_len(nrow(df))) {
    row <- df[k,]
    band <- efficacyPredictions(theta = theta,
                                logTitresDistribution = log10(columnToRealVector(row$Titres)),
                                probability = ci)
    current <- tibble(
      Variant = row$Variant,
      Dose = row$Dose,
      GMT = row$GMT,
      Mean = band$Mean,
      Mean_Low = band$Mean_Low,
      Mean_Up = band$Mean_Up)

    prob_df <- rbind(prob_df, current)
  }
  if (orderGMT)
  { prob_df[order(prob_df$GMT),] }
  else { prob_df }
}

#### Fit_df is used to fit the model; prediction_df is where we predict VE for the  variant; exclude_variant_fit = T means
##### we are fitting the model without data on the variant.
####### We alsways fit a single model. But predict we can for multiple variants from the model.
########## If exclude_variant_fit, than all predicted variants will be excluded from fitted model.
predictVariant <- function(fit_df, variant, prediction_df = NULL, exclude_variant_fit = T, orderGMT = T) {

  if (is.null(prediction_df)) { prediction_df <- fit_df }
  if (exclude_variant_fit)
  { fit_df <- filter(fit_df, !(Variant %in% variant)) }
  est <- fitEfficacyFullOutput(fit_df)
  theta <- c(est$Alpha[1], est$Beta[1], est$Phi[1])

  redux <- prediction_df %>%
    filter(Variant %in% variant) %>%
    distinct(Dose, Variant, .keep_all = T)   ####  NB! would  exclude same-dose cohorts, even if different

  prob_df <- expectedMeanForTitres(redux, theta = theta, orderGMT = orderGMT)
  #### meta info about predictions.
  prob_df$Name <- prediction_df$Vaccine[1]
  prob_df$Theta <- list(theta)
  prob_df$SE <- list(c(est$Alpha[2], est$Beta[2], est$Phi[2]))
  prob_df$TitreCounts <- countPrimaryNabts(fit_df, fit_df$Vaccine[1])
  prob_df$EffCounts <- nrow(fit_df)
  prob_df$Base <- list(unique(fit_df$Variant)) #predict at these variants

  prob_df
}

efficaciesWithTitres <- function(titres = NULL, vaccine = 'Comirnaty', noPriorCovid, ukOnly = F, doImpute = T) {
  if (is.null(titres))
  { titres <- getEarlyVariantTitres(vaccine = vaccine, noPriorCovid = noPriorCovid) }
  effs <- getEfficacies(vaccine = vaccine, isOnlyUK = ukOnly) %>% filter(Dose %in% unique(titres$Dose))
  ans <- combineEfficaciesTitres(efficacies = effs,
                                 titres = titres)

  if (doImpute) { imputeHelper(ans) } else { ans }
}

#first titres are used to fit distributions; 2nd, reduced titres, to evaluate the
# fitted distribution
imputeHelper <- function(titre_df, titre_redux = NULL, assignReduxOmi_ID = F) {

  if (is.null(titre_redux)) {
    titre_redux <- titre_df
  }
  eff <- titre_df
  for (k in seq_len(nrow(eff))) {
    row <- eff[k,]
    up_titre_multiplier <- if (row$Dose == 1) { 3 }else { 5 }
    augmented_x <- imputeTitres(base = unlist(row$Titres), imputed = unlist(
      (titre_redux[k,])$Titres), up_titre_multiplier = up_titre_multiplier)
    row$GMT <- Gmean(augmented_x)
    row$Titres <- list(augmented_x)
    if (assignReduxOmi_ID) {
      row$Omi_ID <- (titre_redux[k,])$Omi_ID }
    eff[k,] <- row
  }
  eff
}

resultsWithObservedMeanVE <- function(observations, results) {
  rdx <- select(observations, Vaccine, Dose, Variant, Mean) %>%
    distinct(across(-Mean), .keep_all = T) %>%
    rename('Mean_VE' = 'Mean')
  left_join(results, rdx, by = c('Vaccine', 'Dose', 'Variant'))
}

###### Computes average VE per Vaccine, variant and dose.
withMeanObservedVE <- function(df) {

  local <- function(vaccine, variant, dose) {
    delta_observed <- df %>% filter(Vaccine == vaccine, Variant == variant, Dose == dose)
    if (is.null(delta_observed)) { return(NULL) }

    eff_two <- mean(delta_observed$Efficacy)
    delta_observed$Mean <- unlist(rep(eff_two, nrow(delta_observed)))
    delta_observed
  }

  ans <- NULL
  for (vaccine in unique(df$Vaccine))
  { for (variant in unique(df$Variant)) {
    for (dose in unique(df$Dose))
    { ans <- rbind(ans, local(vaccine = vaccine, variant = variant, dose = dose)) }
  }
  }
  ans
}

################ MAIN FUNCTIONS

mainFittedModel <- function(df, predict_df = NULL, is_plotted = T, uniq_variants_doses = T, fit_only_variants = NULL, show_estimates = F) {

  fit <- if (is.null(fit_only_variants)) { df }
  else { filter(df, Variant %in% fit_only_variants) }
  theta <- fitEfficacyModel(fit)

  if (is.null(predict_df)) { predict_df <- df }
  redux <- if (uniq_variants_doses)
  { predict_df %>% distinct(Dose, Variant, .keep_all = T) }
  else { predict_df }
  prob_df <- expectedMeanForTitres(redux, theta)

  if (show_estimates) { print(theta) }
  if (is_plotted)
  { plotFittedModelLog(prob_df = prob_df, efficacies_df = df) }

  prob_df
}

mainFittedModelVaxzevria <- function(df, df_com) {

  prob_df <- mainFittedModel(df, is_plotted = F)
  theta <- fitEfficacyModel(df_com)
  redux <- df %>% distinct(Dose, Variant, .keep_all = T)
  prob_df_com <- expectedMeanForTitres(redux, theta)
  plotVaxzevriaModel(prob_df = prob_df, efficacies_df = df, extra_model_df = prob_df_com)
}


bothVaccinesModelFit <- function(noPriorCovid = F, doImpute = T) {
  az <- efficaciesWithTitres(vaccine = 'other', noPriorCovid = noPriorCovid, doImpute = doImpute)
  pf <- efficaciesWithTitres(noPriorCovid = noPriorCovid, doImpute = doImpute)
  both <- rbind(pf, az)
  prob_df <- mainFittedModel(both, is_plotted = F, uniq_variants_doses = F)
  plotBothFittedModels(prob_df = prob_df, efficacies_df = both)
  prob_df
}


kidneyDiseaseSubgroupsComirnaty <- function(noPriorCovid = T, do_impute = T, draws = 100000) {
  eff_all <- efficaciesWithTitres(noPriorCovid = noPriorCovid, doImpute = do_impute)
  theta <- fitEfficacyModel(eff_df = eff_all)

  analyseClinicalSubgroups <- function(disease) {
    titres <- comirnatyClinicalSubgroupTitres(excludeKnownCovid = noPriorCovid,
                                              disease = disease) %>%
      filter(Dose == 2) %>%
      toPredictionTitres()
    if (do_impute) { titres <- imputeHelper(titres) }
    prob_df <- expectedMeanForTitres(titres, theta, orderGMT = F)
    weights_df <- computeVartiantWeightsAllStudies(disease = disease, dose = 2)
    simulateCombinedPrediction(theta = theta, shares_df = weights_df,
                               mus = prob_df, draws = draws)
  }

  df <- rbind(analyseClinicalSubgroups(disease = 'Diabetes'), analyseClinicalSubgroups(disease = 'Kidney'))
  healthy <- comirnatyClinicalSubgroupEfficacies(disease = 'Healthy') %>% filter(Dose == 2)
  plotKidneyDiabetesSubgroups(df, healthy)
}


predictVariantsComirnaty <- function(variants = c('D614G', 'B.1.1.7', 'B.1.351', 'B.1.617.2'), noPriorCovid = F, ukOnly = F, doImpute = T,
                                     isPlot = F, with_cromer_predictions = F) {
  df <- efficaciesWithTitres(vaccine = 'Comirnaty', noPriorCovid = noPriorCovid, ukOnly = ukOnly, doImpute = doImpute)

  results <- NULL
  for (variant in variants) {
    result <- predictVariant(df, variant = variant, prediction_df = df, exclude_variant_fit = T, orderGMT = F)
    results <- rbind(results, result)
  }
  results$Name <- 'Comirnaty'
  results$Vaccine <- 'Comirnaty'
  if (with_cromer_predictions)
  { cromer <- filter(getPredictionsCromer(), Vaccine == 'Comirnaty', Variant %in% variants)
    results <- rbind.fill(results, cromer) }

  observations <- withMeanObservedVE(df %>% filter(Variant %in% variants))
  results <- resultsWithObservedMeanVE(observations = observations, results = results)

  if (isPlot) {
    plotVariantPredictionsPerDose(observed = observations, intervals = results, dose = 1, variants = variants, show_axis = T)
    plotVariantPredictionsPerDose(observed = observations, intervals = results, dose = 2, variants = variants, show_axis = T)
  }
  results
}


runOnlyUK <- function(vaccine = 'Comirnaty', noPriorCovid = T, doImpute = T,
                      isPlot = T, show_model = F) {
  if (show_model)
  { com <- efficaciesWithTitres(vaccine = vaccine, noPriorCovid = noPriorCovid,
                                ukOnly = T, doImpute = doImpute)
    mainFittedModel(df = com, is_plotted = T) }
  if (vaccineName(vaccine) == 'Comirnaty')
  { predictVariantsComirnaty(variants = c('B.1.1.7', 'B.1.617.2'), noPriorCovid = noPriorCovid,
                             ukOnly = T, doImpute = doImpute, isPlot = isPlot, with_cromer_predictions = T) }
  else { predictVaxzevriaVariants(variants = c('B.1.1.7', 'B.1.617.2'), includeCom = F, noPriorCovid = noPriorCovid,
                                  ukOnly = T, doImpute = doImpute, isPlot = isPlot, with_cromer_predictions = T) }
}


#NB! Not super-efficient algorithm. Unneeded computations
runAgePrediction <- function(noPriorCovid = T, doImpute = T) {

  titres <- getEarlyVariantTitres(noPriorCovid = noPriorCovid)
  eff <- efficaciesWithTitres(titres = titres, noPriorCovid = noPriorCovid, doImpute = doImpute)
  theta <- fitEfficacyModel(eff)
  prob_df <- expectedMeanForTitres(distinct(eff, Dose, Variant, .keep_all = T), theta)
  eff_age <- comirnatyAgeEfficaciesAndrews() %>% filter(Variant %in% c('B.1.1.7', 'B.1.617.2'))

  #NB! possibly slight inconsistency when imputing. The model is fit with all titres imputed; but predictions are made for subgroup-imputed titres; hence, slight difference
  # in rep-points

  #Ages for titres: OLDER, YOUNGER.
  predictForAge <- function(titres, theta, age_separator_inclusive, age_group = 'OLDER') {
    df_age <- titres
    if (age_group == 'OLDER') {
      df_age <- filter(titres, Age >= age_separator_inclusive & Age <= 65)
    }
    if (age_group == 'YOUNGER') {
      df_age <- filter(titres, Age <= age_separator_inclusive)
    }
    age_eff <- efficaciesWithTitres(titres = df_age, doImpute = doImpute)
    redux <- age_eff %>% distinct(Dose, Variant, .keep_all = T)
    expectedMeanForTitres(redux, theta)
  }

  predictions <- tibble(
    Age_Group = numeric(),
    Age_Start = numeric(),
    Age_End = numeric(),
    Vaccine = character(),
    Variant = character(),
    Dose = integer(),
    GMT = numeric(),
    Efficacy = numeric(),
    Mean = numeric(),
    Mean_Low = numeric(),
    Mean_Up = numeric()
  )

  for (k in seq_len(nrow(eff_age))) {
    row <- eff_age[k,]
    age_group <- row$Age_Group
    age_separator <- if (age_group == 'YOUNGER') {
      row$Age_End
    } else {
      row$Age_Start
    }

    pred <- predictForAge(titres, theta = theta, age_separator_inclusive = age_separator,
                          age_group = age_group)
    pred <- pred %>%
      filter(Dose == row$Dose, Variant == row$Variant) %>%
      select(-c(Variant, Dose))
    predictions <- rbind(predictions, cbind(row, pred))
  }

  predictions$Error <- predictions$Mean - predictions$Efficacy
  predictions <- predictions[(order(predictions$GMT)),]
  combine_prob_df <- rbind(prob_df, select(predictions, intersect(colnames(prob_df), colnames(predictions))))
  combine_prob_df <- combine_prob_df[order(combine_prob_df$GMT),] # %>% fil

  plotAgeSubgroups(prob_df = combine_prob_df, age_df = predictions)

  predictions }


predictForEpochsComirnaty <- function(noPriorCovid = T,
                                      doImpute = T, is_plot = F, with_others_predictions = T) {
  df <- efficaciesWithTitres(vaccine = 'Comirnaty', noPriorCovid = noPriorCovid, doImpute = doImpute)
  variants <- c('D614G', 'B.1.1.7', 'B.1.351', 'B.1.617.2') ### Fix all pre-Omicron variants for this investigation task.

  predictionSpecs <- function() {
    spec <- tibble(Epoch = 1,
                   Fitted = list('D614G'),
                   Predicted = list(c('B.1.1.7', 'B.1.351', 'B.1.617.2'))
    ) %>%  ###  Epoch 2
      add_row(
        Epoch = 2,
        Fitted = list(c('D614G', 'B.1.351')),
        Predicted = list(c('B.1.1.7', 'B.1.617.2'))) %>%
      add_row(
        Epoch = 2,
        Fitted = list(c('D614G', 'B.1.1.7')),
        Predicted = list(c('B.1.351', 'B.1.617.2')))

    ####### Add third epoch: "Leave one variant out", predict one variant from the othere 3
    for (predict in variants) {
      spec <- spec %>% add_row(
        Epoch = 3,
        Fitted = list(setdiff(variants, predict)),
        Predicted = list(predict))
    }
    spec
  }

  specs <- predictionSpecs()
  results <- NULL
  for (k in seq_len(nrow(specs))) {
    spec <- specs[k,]
    result <- predictVariant(filter(df, Variant %in% spec$Fitted[[1]]),
                             prediction_df = filter(df, Variant %in% spec$Predicted[[1]]),
                             variant = spec$Predicted[[1]],
                             exclude_variant_fit = T)
    result$Epoch <- spec$Epoch
    results <- rbind(results, result)
  }

  ### NB: rename "Name" to inference type or source or model type. As not clear. Also in Plots and some other methods.
  results$Vaccine <- 'Comirnaty'
  results$Name <- 'Comirnaty'
  if (with_others_predictions)
  { cromer <- getPredictionsCromer() %>%
    filter(Vaccine == 'Comirnaty', Variant %in% variants)
    cromer$Epoch <- 3
    results <- rbind.fill(results, cromer) }

  observations <- withMeanObservedVE(df %>% filter(Variant %in% variants))
  results <- resultsWithObservedMeanVE(observations = observations, results = results)
  if (is_plot)
  { plotComirnatyEpochs(observed = observations, intervals = results, variants = variants, dose = 1)
    plotComirnatyEpochs(observed = observations, intervals = results, variants = variants, dose = 2) }

  results
}


predictVaxzevriaVariants <- function(variants = c('B.1.1.7', 'B.1.617.2'),
                                     noPriorCovid = F, ukOnly = F, doImpute = F, isPlot = F, with_cromer_predictions = F) {
  if ('B.1.351' %in% variants) {
    stop('Oops!!! Cannot predicict VE viz  Beta (B.1.351) here:  No Vaxzevria VE estimate to compare with in the VE dataset.')
  }

  az <- efficaciesWithTitres(vaccine = 'Vaxzevria', noPriorCovid = noPriorCovid, ukOnly = ukOnly, doImpute = doImpute)
  com <- efficaciesWithTitres(vaccine = 'Comirnaty', noPriorCovid = noPriorCovid, ukOnly = ukOnly, doImpute = doImpute)

  results <- NULL
  for (variant in variants) {
    result <- predictVariant(fit_df = az, prediction_df = az, variant = variant, exclude_variant_fit = T)
    result$Name <- 'Vaxzevria_direct'

    result_com <- predictVariant(fit_df = com, prediction_df = az, variant = variant, exclude_variant_fit = T)
    result_com$Name <- 'Vaxzevria_from_Comirnaty'
    results <- rbind(results, result, result_com)

  }
  results$Vaccine <- 'Vaxzevria'
  if (with_cromer_predictions)
  { cromer <- getPredictionsCromer() %>% filter(Vaccine == 'Vaxzevria', Variant %in% variants)
    results <- rbind.fill(results, cromer) }

  observations <- withMeanObservedVE(az %>% filter(Variant %in% variants))
  results <- resultsWithObservedMeanVE(observations = observations, results = results)

  if (isPlot)
  { plotVariantPredictionsPerDose(observed = observations, intervals = results, dose = 1, variants = variants, show_axis = T)
    plotVariantPredictionsPerDose(observed = observations, intervals = results, dose = 2, variants = variants, show_axis = T) }

  results
}

######### vaxzevria_direct_only = F includes indirect Vaxzevria into CCC
concordanceVariantPredictions <- function(with_UK = T, is_print = T, vaxzevria_direct_only = F) {
  com <- predictVariantsComirnaty(noPriorCovid = T, doImpute = T, isPlot = F, with_cromer_predictions = F)
  com_UK <- predictVariantsComirnaty(noPriorCovid = T, doImpute = T, isPlot = F, ukOnly = T,
                                     variants = c('B.1.1.7', 'B.1.617.2'), with_cromer_predictions = F)
  vax <- predictVaxzevriaVariants(variants = c('B.1.1.7', 'B.1.617.2'), noPriorCovid = F, doImpute = F,
                                  isPlot = F, with_cromer_predictions = F)
  if (vaxzevria_direct_only) {
    vax <- filter(vax, Name == 'Vaxzevria_direct')
  }

  combined <- if (with_UK) { rbind(com, com_UK, vax) }
  else { rbind(com, vax) }
  rho <- CCC(combined$Mean_VE, combined$Mean)$rho.c

  if (is_print) {
    print('Vaxzevria')
    print(CCC(vax$Mean_VE, vax$Mean)$rho.c)
    print('Comirnaty')
    print(CCC(com$Mean_VE, com$Mean)$rho.c)
    print('Comirnaty UK')
    print(CCC(com_UK$Mean_VE, com_UK$Mean)$rho.c)
    print('Both: ')
    print(rho)
  }

  plotPredictionConcordance(combined, rho)

  combined
}


################ OMICRON
omicronStudyDivisors <- function(combined_df, theta,
                                 doImpute = F, noPriorCovid = F, union_cohorts
                                   = F, age_matched = F) {
  raw_omi_titres <- omicronTitres(noPriorCovid = noPriorCovid)
  if (age_matched) {
    raw_omi_titres <- filter(raw_omi_titres, age > 44)
  }
  omi_titres <- groupedOmicronTitres(raw_omi_titres)

  combinedMinimiser <- function(early_group, both_ids) {
    variant_other <- if (early_group$Variant == 'B.1.617.2') { 'B.1.1.7' }
    else { 'B.1.617.2' }
    other_group <- groupedOmicronTitres(raw_omi_titres %>% filter(Omi_ID %in%
                                                                    both_ids),
                                        cohorts = 'Early') %>% filter(Variant ==
                                                                        variant_other)

    if (doImpute) {
      early_other <- omi_titres %>% filter(Cohort == 'Early', Variant ==
        variant_other)
      other_group <- imputeHelper(early_other, titre_redux = other_group,
                                  assignReduxOmi_ID = T) }

    early_group <- rbind(early_group, other_group)
    initial_alpha <- columnToRealVector(filter(combined_df, Dose == 2,
                                               Variant == 'B.1.1.7')$Titres[[1]])
    initial_delta <- columnToRealVector(filter(combined_df, Dose == 2,
                                               Variant == 'B.1.617.2')$Titres[[1]])
    alpha <- columnToRealVector((early_group %>% filter(Variant ==
                                                          'B.1.1.7'))$Titres[[1]])
    delta <- columnToRealVector((early_group %>% filter(Variant ==
                                                          'B.1.617.2'))$Titres[[1]])

    minimiseLogisticDistanceTwoVariants(titresA = initial_alpha,
                                        titresBisA = alpha, titresB = initial_delta,
                                        titresBisB = delta, theta = theta)$minimum
  }

  #variant_reference: provides titres to compute the multiplier;
  #variant_predict: for which we compute the multiplier
  #NB! titres and Omi_ID do not map due to sorting in imputation
  makeRows <- function(variant_reference, variant_predict, other_cohort) {
    #From pre-omicron study. Assumed comparable weeks to omicron Early: e.g. 2 - 6
    # weeks or 2 - 9 weeks
    initial_early <- columnToRealVector(filter(combined_df, Dose == 2,
                                               Variant == variant_reference)$Titres[[1]])
    #Early Cohort from omicron study. no NAs in early cohort
    early <- omi_titres %>% filter(Cohort == 'Early',
                                   Variant == variant_reference)
    other <- omi_titres %>% filter(Cohort == other_cohort, Variant == variant_predict)
    both_ids <- if (!union_cohorts) { unique(intersect(early$Omi_ID[[1]],
                                                       other$Omi_ID[[1]])) }
    else { unique(union(early$Omi_ID[[1]], other$Omi_ID[[1]])) }

    early_group <- groupedOmicronTitres(raw_omi_titres %>% filter(Omi_ID %in% both_ids),
                                        cohorts = 'Early') %>% filter(Variant == variant_reference)
    other_predict <- groupedOmicronTitres(raw_omi_titres %>% filter(Omi_ID %in% both_ids),
                                          cohorts = other_cohort) %>% filter(Variant == variant_predict)
    if (doImpute) {
      early_group <- imputeHelper(early, titre_redux = early_group, assignReduxOmi_ID = T)
      other_predict <- imputeHelper(other, titre_redux = other_predict, assignReduxOmi_ID = T)
    }

    vector_early <- columnToRealVector(early_group$Titres[[1]])
    vector_other <- columnToRealVector(other_predict$Titres[[1]])
    minimiser_both <- if (variant_predict == 'B.1.1.529') {
      combinedMinimiser(early_group = early_group, both_ids = both_ids)
    } else { NA }

    tibble(
      Cohort = other_cohort,
      Variant = variant_predict,
      Variant_Reference = variant_reference,
      Imputed = doImpute,
      Union_Cohorts = union_cohorts,
      Age_Matched = age_matched,
      Gmt_Ratio = Gmean(vector_early) / Gmean(initial_early),
      Med_Ratio = median(vector_early) / median(initial_early),
      Minimiser_Ref = minimiseLogisticDistance(initial_early, vector_early, theta)
        $minimum,
      Minimiser_Combined = minimiser_both,
      Sample_Size = length(other_predict$Omi_ID[[1]]),
      Titres = list(vector_other),
      Omi_ID = list(other_predict$Omi_ID[[1]]))
  }

  ans <- NULL
  for (predict in c('B.1.617.2', 'B.1.1.529')) {
    for (cohort in c('Early', 'Late', 'Boost')) {
      for (reference in c('B.1.1.7', 'B.1.617.2')) {
        ans <- rbind(ans, makeRows(variant_predict = predict, variant_reference =
          reference, other_cohort = cohort))
      }
    }
  }
  ans
}

predictFromOmicronStudy <- function(theta, multiples, unique_cohorts = T) {

  updateTitres <- function() {
    cohorts <- if (unique_cohorts) {
      unique(multiples$Cohort)
    } else { multiples$Cohort }

    grouped <- NULL
    for (cohort in cohorts) {
      dose <- if (cohort == 'Boost') { 3 } else { 2 }
      now <- multiples %>% filter(Cohort == cohort)
      for (variant in unique(now$Variant)) {
        row <- now %>% filter(Variant == variant)
        kappa <- if (variant == 'B.1.1.529') { row$Minimiser_Combined }
        else { row$Minimiser_Ref }

        titres <- columnToRealVector(row$Titres) / kappa
        gmt <- Gmean(titres)
        grouped <- rbind(grouped, tibble(
          Vaccine = 'Comirnaty',
          Cohort = cohort,
          Variant = variant,
          Dose = dose,
          Divisor = kappa,
          Imputed = row$Imputed,
          GMT = gmt,
          Omi_ID = row$Omi_ID,
          Titres = list(titres))) }
    }

    grouped }

  grouped <- updateTitres()
  prob_df <- expectedMeanForTitres(grouped, theta, orderGMT = F)
  prob_df$Cohort <- grouped$Cohort

  prob_df
}


runOmicronDeltaPredictionDifferences <- function(draws = 100000, noPriorCovid = F) {
  eff <- comirnatyOmicronDeltaEfficacies()
  omi_titres <- omicronTitres(noPriorCovid = noPriorCovid)
  grouped_omi <- groupedOmicronTitres(omi_titres)

  df_delta <- omicronEstimationTable(eff = eff, grouped_omi = grouped_omi, is_delta = T)
  theta_delta <- fitEfficacyModel(df_delta)

  df_omi <- omicronEstimationTable(eff = eff, grouped_omi = grouped_omi, is_delta = F)
  theta_omi <- fitEfficacyModel(df_omi)

  cohorts <- c('Early', 'Late', 'Boost', 'BoostHD')
  difs <- NULL
  for (cohort in cohorts) {
    #Print for Delta II Early for comparison
    if (cohort == 'Early') {
      delta_titres <- filter(grouped_omi, Variant == 'B.1.617.2', Cohort == cohort)$Titres[[1]]
      delta_gmt <- bootGmtMeanCI(delta_titres, draws = draws)
      print('Early II Delta NAbTs: ')
      print(delta_gmt)
    }

    if (cohort == 'BoostHD') {
      omi_hd <- groupedOmicronDialysisTitres(omicronDialysisTitres())
      titres <- filter(omi_hd, Variant == 'B.1.1.529', Cohort == 'Boost')$Titres[[1]]
    } else {
      titres <- filter(grouped_omi, Variant == 'B.1.1.529', Cohort == cohort)$Titres[[1]]
    }
    dif <- simulatePredictionDifferenceOmicronDelta(theta_Delta = theta_delta, theta_Omicron = theta_omi,
                                                    log_titres_omi = log10(titres), draws = draws)
    difs <- rbind(difs, dif)
  }
  difs$Cohort <- c('Early II', 'Late II', 'Boost III', 'HD Boost III')
  tableDeltaOmicronDifference(difs)
}

runDeltaOmicronScaled <- function(noPriorCovid = T, doImpute = T, display_summary = F) {
  # ####################### Titres for Study One, used for fitting the model
  titres_one <- getEarlyVariantTitres(noPriorCovid = noPriorCovid) %>% filter(Dose == 1 | DaysSinceJab2 <= 42)
  combined_one <- efficaciesWithTitres(titres = titres_one, noPriorCovid = noPriorCovid, doImpute =
    doImpute)
  theta <- fitEfficacyModel(eff_df = combined_one)
  ###### divisors for cohorts
  kappas_agm <- omicronStudyDivisors(combined_one, theta = theta,
                                     doImpute = doImpute, noPriorCovid =
                                       noPriorCovid,
                                     union_cohorts = T, age_matched = F)

  if (display_summary) { print(kappas_agm) }

  meanEfficacyPerCohortVariant <- function(omi_pred) {
    ans <- NULL
    for (variant in unique(omi_pred$Variant)) {
      rows <- omi_pred %>% filter(Variant == variant)
      for (cohort in unique(rows$Cohort)) {
        current <- rows %>% filter(Cohort == cohort)
        current$Efficacy_Mean <- mean(current$Efficacy)
        current$Is_Singleton <- nrow(current) == 1
        ans <- rbind(ans, current)
      }
    }
    ans
  }

  multiples <- filter(kappas_agm, Variant_Reference == 'B.1.617.2')
  effs <- comirnatyOmicronDeltaEfficacies() %>%
    meanEfficacyPerCohortVariant()

  combined <- combineOmicronEfficacyWithPredictions(effs, predictFromOmicronStudy(multiples = multiples, theta = theta))
  redux <- combined_one %>% distinct(Dose, Variant, .keep_all = T)
  prob_df <- rbind(expectedMeanForTitres(redux, theta), select(combined, Variant, Dose, GMT,
                                                               Mean, Mean_Low, Mean_Up))
  prob_df <- prob_df[order(prob_df$GMT),]

  plotDeltaOmicronScaledBoth(prob_df = prob_df, both = combined)
}


runDeltaOmicronModel <- function(noPriorCovid = F) {
  eff <- comirnatyOmicronDeltaEfficacies()
  omi <- omicronTitres(noPriorCovid = noPriorCovid)
  grouped_omi <- groupedOmicronTitres(omi)

  df_delta <- omicronEstimationTable(eff = eff, grouped_omi = grouped_omi, is_delta = T)
  prob_delta <- mainFittedModel(df_delta, uniq_variants_doses = F, is_plotted = F)

  df_omi <- omicronEstimationTable(eff = eff, grouped_omi = grouped_omi, is_delta = F)
  prob_omi <- mainFittedModel(df_omi, uniq_variants_doses = F, is_plotted = F)

  hd_titres <- omicronDialysisTitres()
  hd <- groupedOmicronDialysisTitres(hd_titres)
  hd_delta <- mainFittedModel(df_delta, predict_df = filter(hd, Variant == 'B.1.617.2'), uniq_variants_doses = F, is_plotted = F)

  hd_omi <- mainFittedModel(df_omi, predict_df = filter(hd, Variant == 'B.1.1.529'), uniq_variants_doses = F, is_plotted = F)
  hd_predict <- rbind(hd_delta, hd_omi)

  plotDeltaOmicronFitted(prob_delta = prob_delta, eff_delta = df_delta,
                         prob_omi = prob_omi, eff_omi = df_omi, hd = hd_predict)
}

#################### Supplementary Figures. NB. Code for Figure 6  is missing but not important.
################## Also, similar code is accessed from  the next funcion
runCensoredDistributions <- function(vaccine = 'Comirnaty', dose = 2, showVariant = 'B.1.351', noPriorCovid = F, is_print = F) {
  titres <- getEarlyVariantTitres(vaccine = vaccine, noPriorCovid = noPriorCovid)
  x <- (titres %>% filter(Dose == dose))[(as.symbol(showVariant))]

  cens <- titresToCensoredDF(x)
  best <- fitCensoredDistribution(cens, show_summary = F)

  qfun <- curryCDF(best)
  reps <- censoredRepPoints(titres = x, qfun = qfun)
  reps10 <- 10^reps
  middle <- columnToRealVector(x)
  sum_1 <- summary(middle)
  gmt_before <- Gmean(middle)
  middle <- middle[40 < middle & middle < 2560]
  middle <- c(reps10, middle)
  gmt_after <- Gmean(middle)
  sum_2 <- summary(middle)

  if (is_print) {
    print('------')
    print(best)
    print(sum_1)
    print(gmt_before)
    print(sum_2)
    print(gmt_after)
  }

  dfun <- curryPDF(best)
  xs <- seq(from = 0.01, to = log10(5760), by = 0.01)
  pdfs <- tibble(
    Titres = list(xs),
    Variant = showVariant,
    Distribution = best$distname,
    Density = list(dfun(xs) / 3),
    Reps = list(reps)
  )
  plotEarlyTitresHistogram(titres %>% filter(Dose == dose), pdfs =
    pdfs, singleVariant = showVariant)
}


simplifiedVersusFullModels <- function(df) {
  betafit <- fitBetaRegression(df, estimation_type = 'BC')
  theta <- optim(par = c(betafit$coefficients$mean, betafit$coefficients$precision), fn
    = betaLikelihoodData, df = df)$par

  prob_df <- expectedMeanForTitres(df %>% distinct(Dose, Variant, .keep_all = T), theta)
  beta <- betafit$coefficients$mean
  simpleFun <- Curry(logisticFunction, alpha = beta[1], beta = beta[2])

  plotSimpleFullModelsLog(prob_df = prob_df, efficacies_df = df, simple_fun = simpleFun)
}


################ Estimates all the models and prints them in LaTeX tables:

fittedPrimaryModelsTables <- function()
{

  countEffs <- function(df, vaccine = 'Comirnaty') {
    nrow(filter(df, Vaccine == vaccine))
  }

  basics <- function(df, noPriorCovid, doImpute, head = '') {
    out <- fitEfficacyFullOutput(df)

    tibble(Model = head,
           Alpha = meanSe(out$Alpha),
           Beta = meanSe(out$Beta),
           Phi = meanSe(out$Phi),
           ComNab = countPrimaryNabts(df),
           VaxNab = countPrimaryNabts(df, 'Vaxzevria'),
           Dummy = '',
           VE = countEffs(df) + countEffs(df, vaccine = 'Vaxzevria'),
           Covid = if (noPriorCovid) { 'N' } else { 'Y' },
           Imputed = if (doImpute) { 'Y' } else { 'N' })
  }

  com_df <- efficaciesWithTitres(noPriorCovid = T, doImpute = T)
  com_basic <- basics(head = 'Comirnaty', df = com_df, doImpute = T, noPriorCovid = T)

  vax_indirect <- basics(head = 'Vaxzevria Indirect', df = efficaciesWithTitres(noPriorCovid = F, doImpute = F), doImpute = F, noPriorCovid = F)
  vax_direct <- basics(head = 'Vaxzevria Direct', df = efficaciesWithTitres(vaccine = 'Vaxzevria', noPriorCovid = F, doImpute = F), doImpute = F, noPriorCovid = F)
  az <- efficaciesWithTitres(vaccine = 'Vax', noPriorCovid = F, doImpute = T)
  pf <- efficaciesWithTitres(noPriorCovid = F, doImpute = T)
  both <- basics(head = 'Combined', df = rbind(pf, az), doImpute = T, noPriorCovid = F)

  titres_one <- getEarlyVariantTitres(noPriorCovid = T) %>% filter(Dose == 1 | DaysSinceJab2 <= 42)
  combined_one <- efficaciesWithTitres(titres = titres_one, noPriorCovid = T, doImpute = T)

  com_rsc <- basics(head = 'Comirnaty\\textdagger', df = combined_one, doImpute = T, noPriorCovid = T)
  res <- rbind(com_basic, vax_indirect, vax_direct, both, com_rsc)
  tableFullPrimaryModels(res)
}

deltaOmicronModelTables <- function(noPriorCovid = F) {
  #Columns: Head, alpha, beta, phi, EarlyTwo, LateTwo, BoosterThree, VE, Covid/Imputation,  Note

  countNabts <- function(df) {
    df <- distinct(df, Cohort, .keep_all = T)
    early <- sapply(filter(df, Cohort == 'Early')$Titres, length)
    late <- sapply(filter(df, Cohort == 'Late')$Titres, length)
    boost <- sapply(filter(df, Cohort == 'Boost')$Titres, length)

    tibble(Early = early,
           Late = late,
           Boost = boost,
           Total = early + late + boost)
  }

  basics <- function(df, noPriorCovid, doImpute, head = '') {
    nab <- countNabts(df)
    out <- fitEfficacyFullOutput(df)

    tibble(Model = head,
           Alpha = meanSe(out$Alpha),
           Beta = meanSe(out$Beta),
           Phi = meanSe(out$Phi),
           EarlyTwo = nab$Early,
           LateTwo = nab$Late,
           BoostThree = nab$Boost,
           Total = nab$Total,
           Eff = nrow(df),
           Covid = if (noPriorCovid) { 'N' } else { 'Y' },
           Impute = if (doImpute) { 'Y' } else { 'N' })
  }

  eff <- comirnatyOmicronDeltaEfficacies()
  grouped_omi <- groupedOmicronTitres(omicronTitres(noPriorCovid = noPriorCovid))

  df_delta <- omicronEstimationTable(eff = eff, grouped_omi = grouped_omi, is_delta = T)
  delta <- basics(df_delta, head = 'Delta', noPriorCovid = noPriorCovid, doImpute = F)

  df_omi <- omicronEstimationTable(eff = eff, grouped_omi = grouped_omi, is_delta = F)
  omicron <- basics(df_omi, head = 'Omicron', noPriorCovid = noPriorCovid, doImpute = F)
  df <- rbind(delta, omicron)

  tableDeltaOmicronModels(df)
}

latexPrimaryVariantPredictions <- function() {
  # Comirnaty
  full <- predictForEpochsComirnaty(noPriorCovid = T, doImpute = T, is_plot = F, with_others_predictions = F) %>%
    distinct(Variant, Base, .keep_all = T)
  uk <- predictVariantsComirnaty(variants = c('B.1.1.7', 'B.1.617.2'), noPriorCovid = T,
                                 ukOnly = T, doImpute = T, isPlot = F, with_cromer_predictions = F)
  uk$Epoch <- 4

  vax <- predictVaxzevriaVariants(noPriorCovid = F, doImpute = F, isPlot = F, with_cromer_predictions = F) %>%
    filter(Dose == 1)
  vax$Epoch <- c(5, 6, 5, 6)
  tablesVariantPredictions(com = rbind(full, uk), vax = vax)
}


############ Runs all figures in the paper  (and tables, if required.)
paperFiguresTables <- function(show_figures = T, show_tables = T) {
  ########## PAPER FIGURES
  if (show_figures)
  { #### Figure 1: Comirnaty, Vaxzevria, Age Comirnaty, Both Models
    mainFittedModel(efficaciesWithTitres(noPriorCovid = T, doImpute = T))
    mainFittedModelVaxzevria(efficaciesWithTitres(vaccine = 'Vaxzevria', noPriorCovid =
      F, doImpute = F), df_com = efficaciesWithTitres(vaccine = 'Comirnaty', noPriorCovid = F, doImpute
      = F))
    runAgePrediction(noPriorCovid = T, doImpute = T)
    bothVaccinesModelFit(noPriorCovid = F, doImpute = T)

    ######### Fig 2: Disease subgroups
    kidneyDiseaseSubgroupsComirnaty()

    #############  Variant Predictions
    predictVaxzevriaVariants(variants = c('B.1.1.7', 'B.1.617.2'), noPriorCovid = F, doImpute = F, isPlot = T, with_cromer_predictions = T) #  Fig 3: Dose 2, Sup Fig 9: Dose 1
    predictForEpochsComirnaty(noPriorCovid = T, doImpute = T, is_plot = T) # Fig 3: Dose 2; Sup Fig 9 Dose 1
    concordanceVariantPredictions(with_UK = T) # Fig  4: Concordance incl. UK only for Comirnaty

    ############### Omicron and Delta
    runDeltaOmicronScaled() # Fig 5: Top
    runDeltaOmicronModel(noPriorCovid = F)  # Fig 5: Bottom

    ############ SUPPLEMENT FIGURES
    runCensoredDistributions() # Fig 7: Illustrate Imputation for Comirnaty 2 doses vs Beta
    simplifiedVersusFullModels(efficaciesWithTitres(noPriorCovid = T, doImpute = T)) # Fig 8: Comirnaty

    ################## Comirnaty: Variant Predictions UK only
    runOnlyUK(noPriorCovid = T, doImpute = T) # Fig 10.
    ########### Vaxzevria suggested:
    # runOnlyUK(vaccine = 'Vaxzevria', noPriorCovid = T, doImpute = F) ####
  }

  ############## PAPER TABLES:  VE estimates and fitted models (LaTeX)
  if (show_tables)
  {
    runOmicronDeltaPredictionDifferences(noPriorCovid = F) # Table 1: Predicted VE differences: Omicron vs Delta models
    efficacyTable(df = getEfficacies()) # Table 2
    efficacyTable(df = getEfficacies(vaccine = 'Vaxzevria')) # Table 3
    tableAgeEfficacies() # Table 4: age-based VE
    tableDiseaseSubgroups(df = frameDiseaseVE()) # Table 5: VE for kidney diseas subgroups
    tableEfficacyDeltaOmicron() #Table 6: Delta and Omicron VE
    fittedPrimaryModelsTables() # Table 7;  pre-Omicron fitted models
    deltaOmicronModelTables() # Table 8;  Delta and Omicron fitted models
    latexPrimaryVariantPredictions() #### Tables 9 and 10 for variant predictions
  }
}

#### Runs all figures and tables (optional) in the paper.
paperFiguresTables()






















