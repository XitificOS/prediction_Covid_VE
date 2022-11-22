library('tidyverse')
library('stringr')
library('DescTools')
library('dplyr')
library("plyr")
library('functional')
library('purrr')

source('fitters.R')
source('selectors.R')

casesForPeriods <- function(variant_dates, daily_covids) {

  toDateIntervals <- function(dates) {
    len <- length(dates)
    if (len <= 1) {
      return(NULL)
    }
    ends <- dates[2:len] - 1
    ends[len - 1] <- dates[len]
    tibble(Start = dates[1:len - 1], End = ends)
  }

  periods <- toDateIntervals(variant_dates)
  totals <- NULL
  for (i in seq_len(nrow(periods))) {
    row <- periods[i,]
    covids <- filter(daily_covids, row$Start <= Date & Date <= row$End)
    totals <- c(totals, sum(covids$Cases))
  }

  periods$Cases <- totals
  periods
}

computeVariantShares <- function(variants_table, periods) {

  decodeCovariantName <- function(covariant) {
    switch(covariant,
           "20H (Beta, V2)" = 'B.1.351',
           "20I (Alpha, V1)" = 'B.1.1.7',
           "21A (Delta)" = 'B.1.617.2',
           "21I (Delta)" = 'B.1.617.2',
           "21J (Delta)" = 'B.1.617.2',
           'D614G'
    )
  }

  zeroVariantsTable <- function() {
    tibble(
      'B.1.617.2' = 0,
      'B.1.1.7' = 0,
      'B.1.351' = 0,
      'D614G' = 0
    )
  }

  variant_frequences <- NULL
  for (i in seq_len(nrow(variants_table))) {
    row <- zeroVariantsTable()
    for (j in seq_len(ncol(variants_table))) {
      strain <- colnames(variants_table)[j]
      count <- variants_table[i, j]
      legacy_variant <- decodeCovariantName(strain)
      row[(legacy_variant)] <- row[(legacy_variant)] + count
    }
    row <- row / sum(row) * periods$Cases[i]
    variant_frequences <- rbind(variant_frequences, row)
  }
  colSums(variant_frequences) / sum(variant_frequences)

}

variantWeightPerStudy <- function(study) {

  variantDates <- function(variants) {

    nearestDateIndex <- function(dates, to_date) {
      which(abs(dates - to_date) == min(abs(dates - to_date), na.rm = TRUE))[1]
    }

    study_start <- as.Date(study$Date_Start, format = '%d %b %Y')
    study_end <- as.Date(study$Date_End, format = '%d %b %Y')
    variant_dates <- as.Date(variants$week)
    variant_dates[nearestDateIndex(variant_dates, study_start)
                    :nearestDateIndex(variant_dates, study_end)]
  }

  variants <- getCovariants(country = study$Country)
  variant_dates <- variantDates(variants = variants)
  daily_covids <- dailyCovidsPerCountry(study$Country) %>%
    filter(variant_dates[1] <= Date & Date <= tail(variant_dates, n = 1))

  variants <- filter(variants, variant_dates[1] <= week & week < tail(variant_dates, n
    = 1))

  periods <- casesForPeriods(variant_dates = variant_dates, daily_covids = daily_covids)
  computeVariantShares(variants_table =
                         variants$stand_estimated_cases,
                       periods = periods)
}

#'Diabetes', 'Kidney'
computeVartiantWeightsAllStudies <- function(
  disease = 'Diabetes', dose = 2, exclude_healthy_Whittaker = T) {

  eff <- comirnatyClinicalSubgroupEfficacies(disease = disease, exclude_healthy_Whittaker = exclude_healthy_Whittaker) %>% filter(Dose == dose)
  variant_weights <- NULL
  for (i in seq_len(nrow(eff))) {
    study <- eff[i,]
    variant_weights <- rbind(variant_weights, variantWeightPerStudy(study))
  }
  cbind(eff, variant_weights)
}

frameDiseaseVE <- function() {
  shares <- rbind(computeVartiantWeightsAllStudies(disease = 'Diabetes'),
                  computeVartiantWeightsAllStudies(disease = 'Kidney'),
                  computeVartiantWeightsAllStudies(disease = 'Healthy', exclude_healthy_Whittaker = F))
  shares[(order(shares$Ref)),]
}

