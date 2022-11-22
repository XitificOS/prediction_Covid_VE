library('gtools')
library('DescTools')
library('betareg')
library('functional')
library('transport')
library('boot')
library('fitdistrplus')
library('actuar')
library('HelpersMG')

# File to fit model severe model, including step
# x is  logtitres. Epoch is at which step the variant appeared?
#e.g., 0 for D614G; 1 for Alpha and Beta; 2 for Delta; 3 for Omicron
logisticFunction <- function(logTitres, alpha, beta) {
  term <- exp(alpha + beta * logTitres)
  term / (1.0 + term)
}

# estimation_type: 'ML', 'BC' or 'BR', ie, standard max likelihood,
# or Bias-correrted, or Bias-reduced. Using
# Intercept only for precision
fitBetaRegression <- function(efficacies_df, estimation_type =
  'BC') {
  betareg(Probability ~ log10(GMT) | 1, link.phi = 'log',
          data = efficacies_df, type = estimation_type)
}

# integrating wrt empirical distribution of logtitres per variant
# for a given dose regimen
populationProtectionPerPoint <- function(logTitres, alpha,
                                         beta) {
  probabilities <- logisticFunction(logTitres = logTitres,
                                    alpha = alpha, beta = beta)
  mean(probabilities, na.rm = T)
}

#LOG LIKELIHOOD FUNCTION!!!
# input logTitres distribution per variant/odose;
#y is the efficacy (0,1) for variant/dose;
#alpha, beta are parameters for the mean, ie, logisticTwo
#a, b are precision parameters, where phi is modelled via explicit
# model
# if theta has 3 parameters, then use constant phi
#Theta has alpha, beta, a and maybe b.

# computes log-likelihood for logTitres distributions of variants.
# Theta includes mean and precision parameters
betaLikelihoodData <- function(theta, df) {

  betaLikelihoodPerPoint <- function(logTitres, y) {
    mu <- populationProtectionPerPoint(logTitres,
                                       alpha = theta[1], beta =
                                         theta[2])
    # the last parameter is for phi! from estimation
    phi <- theta[length(theta)]
    # ref gnlr.R line 1106
    m <- mu * phi
    s <- (1 - mu) * phi
    dbeta(y, m, s, ncp = 0, log = TRUE)
  }

  wt <- 1 / nrow(df)
  lik <- 0
  for (k in seq_len(nrow(df))) {
    row <- df[k,]
    distribution <- log10(columnToRealVector(row$Titres))
    lik <- lik - wt * betaLikelihoodPerPoint(distribution, y =
      row$Probability)
  }
  lik
}

efficacyPredictions <- function(theta,
                                logTitresDistribution, probability = 0.95) {
  mu <- populationProtectionPerPoint(logTitres = logTitresDistribution,
                                     alpha = theta[1],
                                     beta = theta[2])

  phi <- theta[length(theta)]
  p <- mu * phi
  q <- phi - p

  probBound <- (1 - probability) / 2
  probs <- c(probBound, 1 - probBound)
  bands <- qbeta(p = probs, shape1 = p, shape2 = q)

  tibble(
    Mean = 100 * mu,
    Mean_Low = 100 * bands[1],
    Mean_Up = 100 * bands[2]
  )
}


simulatePredictionDifferenceOmicronDelta <- function(theta_Delta, theta_Omicron, log_titres_omi,
                                                     probability = 0.95, draws = 10000) {

  simulate <- function(mu, theta) {
    phi <- theta[length(theta)]
    p <- mu * phi
    q <- phi - p
    rbeta(shape1 = p, shape2 = q, n = draws)
  }

  mu_omi <- efficacyPredictions(logTitresDistribution = log_titres_omi,
                                theta = theta_Omicron, probability = probability)
  mu_delta <- efficacyPredictions(logTitresDistribution = log_titres_omi,
                                  theta = theta_Delta, probability = probability)
  dif <- simulate(mu_delta$Mean / 100, theta = theta_Delta) - simulate(mu_omi$Mean / 100, theta = theta_Omicron)

  probBound <- (1 - probability) / 2
  probs <- c(probBound, 1 - probBound)
  quint <- quantile(dif, probs = probs)
  gmt <- bootGmtMeanCI(vec =10^log_titres_omi, draws = draws,
                       conf= probability)

  tibble(
    GMT = gmt$GMT,
    GMT_Low = gmt$Low,
    GMT_Up = gmt$Up,
    Omi_Mean = mu_omi$Mean,
    Omi_Low = mu_omi$Mean_Low,
    Omi_Up = mu_omi$Mean_Up,
    Delta_Mean = mu_delta$Mean,
    Delta_Low = mu_delta$Mean_Low,
    Delta_Up = mu_delta$Mean_Up,
    TrueDiff = mu_delta$Mean - mu_omi$Mean,
    Dif_Mean = 100 * mean(dif),
    Dif_Median = 100 * median(dif),
    Dif_Low = 100 * quint[1],
    Dif_Up = 100 * quint[2]
  )
}


simulateCombinedPrediction <- function(theta, shares_df, mus, draws = 100000,
                                       probability = 0.95, seed = 42) {
  set.seed(seed)
  variants <- intersect(legacyVariantsAll(), colnames(shares_df))
  shares <- select(shares_df, all_of(variants))
  randoms <- tibble(Dummy = 1:draws)
  phi <- theta[length(theta)]

  for (variant in variants) {
    mu <- filter(mus, Variant == variant)$Mean / 100
    p <- mu * phi
    q <- phi - p

    draw <- rbeta(draws, shape1 = p, shape2 = q)
    randoms <- randoms %>% add_column(!!(variant) := as.numeric(draw))
  }

  result <- NULL
  percentil <- NULL
  alpha <- (1 - probability) / 2
  for (i in seq_len(nrow(shares))) {
    sim_effs <- NULL
    row <- ceiling(shares[i,] * draws)
    row.names(row) <- NULL

    for (variant in variants) {
      count <- as.integer(row[(variant)])
      names(count) <- NULL
      index <- sample(draws, size = count[1], replace = T)
      take <- unlist(randoms[(variant)])
      sim_effs <- c(sim_effs, take[index])
    }

    result <- c(result, mean(sim_effs))
    percentil <- rbind(percentil, quantile(sim_effs, c(alpha, 1 - alpha), na.rm = T))
  }

  ans <- shares_df
  ans$Predicted_Mean <- result * 100
  ans$Predicted_LCL <- percentil[, 1] * 100
  ans$Predicted_UCL <- percentil[, 2] * 100

  ans
}

fitEfficacyModel <- function(eff_df, estimation_type = 'BC') {
  betafit <- fitBetaRegression(efficacies_df = eff_df, estimation_type = estimation_type)
  optim(par = c(betafit$coefficients$mean, betafit$coefficients$precision), fn
    = betaLikelihoodData, df = eff_df)$par

}

fitEfficacyFullOutput <- function(eff_df) {
  betafit <- fitBetaRegression(efficacies_df = eff_df, estimation_type = 'BC')
  out <- optim(par = c(betafit$coefficients$mean, betafit$coefficients$precision), fn
    = betaLikelihoodData, df = eff_df, hessian = T)

  se <- SEfromHessian(out$hessian, hessian = FALSE, silent = FALSE)
  tibble(
    Type = c('Est', 'SE'),
    Alpha = as.numeric(c(out$par[[1]], se[[1]])),
    Beta = as.numeric(c(out$par[[2]], se[[2]])),
    Phi = as.numeric(c(out$par[[3]], se[[3]]))
  )
}


fitCensoredDistribution <- function(cens, show_summary = F) {
  distributions <- c('gamma', 'weibull', 'norm', 'lnorm', 'logis',
                     'exp')
  fits <- NULL
  for (dist in distributions) {
    fit <- fitdistcens(cens, distr = dist)
    if (show_summary) {
      print(summary(fit))
    }
    fits <- rbind(fits, fit)
  }

  fits[which.min(fits[, 'aic']),]
}

curryCDF <- function(best) {

  est <- best$estimate
  switch(
    best$distname,
    "gamma" = Curry(qgamma, shape = est[1], rate = est[2]),
    'gumbel' = Curry(qgumbel, alpha = est[1], scale = est[2]),
    "lgamma" = Curry(qlgamma, shapelog = est[1], ratelog = est[2]),
    'weibull' = Curry(qweibull, shape = est[1], scale = est[2]),
    "norm" = Curry(qnorm, mean = est[1], sd = est[2]),
    "lnorm" = Curry(qlnorm, meanlog = est[1], sdlog = est[2]),
    'logis' = Curry(qlogis, location = est[1], scale = est[2]),
    'exp' = Curry(qexp, rate = est[1])
  )
}

curryPDF <- function(best) {

  est <- best$estimate
  switch(
    best$distname,
    "gamma" = Curry(dgamma, shape = est[1], rate = est[2]),
    'gumbel' = Curry(dgumbel, alpha = est[1], scale = est[2]),
    "lgamma" = Curry(dlgamma, shapelog = est[1], ratelog = est[2]),
    'weibull' = Curry(dweibull, shape = est[1], scale = est[2]),
    "norm" = Curry(dnorm, mean = est[1], sd = est[2]),
    "lnorm" = Curry(dlnorm, meanlog = est[1], sdlog = est[2]),
    'logis' = Curry(dlogis, location = est[1], scale = est[2]),
    'exp' = Curry(dexp, rate = est[1])
  )
}

# qfun inverse distribution function
# todo: is the upper bound ok? maybe higher for two doses?
censoredRepPoints <- function(titres, qfun, low = 40, up = 2560, treshold_min = 1,
                              max_cut_times = 4,
                              show_probs = F) {

  x <- columnToRealVector(titres)
  n <- length(x)
  below_count <- length(which(x < low))
  low_reps <- c()
  counter <- n
  while (length(low_reps) != below_count | counter <= 0) {
    js <- (seq(counter) - 0.5) / counter
    low_reps <- 10^qfun(js[js <= below_count / n])
    low_reps <- low_reps[low_reps <= low]
    if (length(low_reps) < below_count) {
      counter <- counter + 1
    }
    if (length(low_reps) > below_count) {
      counter <- counter - 1
    }
  }

  above_count <- length(which(x > up))
  up_reps <- c()
  counter <- n
  while (length(up_reps) != above_count | counter <= 0) {
    js <- (seq(counter) - 0.5) / counter
    up_reps <- 10^qfun(js[js >= 1 - above_count / n])
    up_reps <- up_reps[up_reps >= up]
    if (length(up_reps) < above_count) {
      counter <- counter + 1
    }
    if (length(up_reps) > above_count) {
      counter <- counter - 1
    }
  }

  reps <- c(low_reps, up_reps)
  if (is.null(reps)) {
    return(NULL)
  }

  if (show_probs) {
    below <- length(which(reps < low))
    above <- length(which(reps > up))
    print(paste0('true below:', below_count, ' ', below))
    print(paste0('true above:', above_count, ' ', above))
  }

  treshold_max <- up * max_cut_times
  replaceThresholds <- function(x) {
    if (x > treshold_max) {
      treshold_max
    }
    else if (x < treshold_min) { treshold_min }
    else { x }
  }

  log10(sapply(reps, replaceThresholds))
}

minimiseLogisticDistance <- function(titresA, titresB, thetaA) {
  alpha <- thetaA[1]
  beta <- thetaA[2]

  a <- logisticFunction(logTitres = log10(titresA), alpha = alpha,
                        beta = beta)
  mean_a <- mean(a)
  distance <- function(divisor) {
    b <- logisticFunction(logTitres = log10(titresB / divisor), alpha = alpha,
                          beta = beta)
    abs(mean_a - mean(b))
  }
  optimize(interval = c(0, 100), f = distance)
}

minimiseLogisticDistanceTwoVariants <- function(titresA, titresBisA, titresB, titresBisB, theta) {
  alpha <- theta[1]
  beta <- theta[2]

  a <- logisticFunction(logTitres = log10(titresA), alpha = alpha,
                        beta = beta)
  b <- logisticFunction(logTitres = log10(titresB), alpha = alpha,
                        beta = beta)

  distance <- function(divisor) {
    aBis <- logisticFunction(logTitres = log10(titresBisA / divisor), alpha = alpha,
                             beta = beta)
    a_dif <- abs(mean(a) - mean(aBis))
    bBis <- logisticFunction(logTitres = log10(titresBisB / divisor), alpha = alpha,
                             beta = beta)
    b_dif <- abs(mean(b) - mean(bBis))
    max(c(a_dif, b_dif))
  }

  optimize(interval = c(0, 100), f = distance)
}

#NB! always orders as sorted. Hence need to track beforehand, eg for subject ID
imputeTitres <- function(base, imputed = NULL, up_titre_multiplier = 5) {
  if (is.null(imputed)) {
    imputed <- base
  }
  middle <- imputed[40 < imputed & imputed < 2560]
  if (length(middle) == length(imputed)) {
    return(sort(middle))

  }
  cens <- titresToCensoredDF(base)
  best <- fitCensoredDistribution(cens, show_summary = F)
  reps <- censoredRepPoints(titres = imputed, qfun = curryCDF(best), max_cut_times =
    up_titre_multiplier)
  sort(c(10^reps, middle))
}

bootGmtMeanCI <- function(vec, draws = 100000, conf = 0.95) {
  b <- boot(data = columnToRealVector(vec), statistic = function(x, i) exp(mean(log(x[i]))), R = draws)
  ans <- boot.ci(b, conf = conf, type = 'bca')
  tibble(GMT = ans$t0,
         Low = ans$bca[[4]],
         Up = ans$bca[[5]])
}