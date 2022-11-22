library('tidyverse')
library('plotly')
library('DescTools')
library('dplyr')
library('plyr')
# library('data.table')
library('functional')

source('selectors.R')
source('fitters.R')
source('presenters.R')


######### Colours used
faintBlue <- function() { '#B0E0E6' }
blue <- function() { '#1f77b4' }
red <- function() { '#d62728' }
green <- function() { '#2ca02c' }

black <- function() { '#000000' }
blueLight <- function() { '#56B4E9' }
grey <- function() { '#BCBCBC' }
greyMedium <- function() { '#808080' }
orange <- function() { "#E69F00" }

# HELPER FUNCTIONS
# works for MacOSs with a chrome browser. Tested in Catalina + zsh. Works for Linux (Ubuntu) with a chrome browser.
# Shows inconsequential error messages.
# Not tested on Windows.
showInChrome <- function(widget, showFigure = T) {
  if (!showFigure) {
    return()
  }
  # Generate random file name
  temp <- paste(tempfile('plotly'), 'html', sep = '.')
  # Save.
  htmlwidgets::saveWidget(widget, temp, selfcontained = FALSE)
  # Show in chrome
  # system(sprintf("open -a 'google chrome'  /%s", temp))
  system(sprintf("google-chrome-stable  /%s", temp))
}

errorBar <- function(low, mid, up, colour = black()) {
  list(
    symmetric = FALSE,
    arrayminus = mid - low,
    array = up - mid,
    color = colour)
}

addVariantAnnotations <- function(fig, df, showDose = F, isLog = F, isBold = F, y_shift
  = 0, distinct_only = T, text_colour = black(), labels = NULL, font_size = 20) {
  redux <- df
  if (distinct_only)
  { redux <- df %>% distinct(Dose, Variant, Vaccine, .keep_all = T) }
  if (is.null(labels))
  { labels <- greekVariantNames(redux$Variant) }
  if (showDose) {
    for (k in seq_len(nrow(redux))) {
      labels[k] <- paste0(labels[k], '  ', switch(redux$Dose[k], 'I', 'II', 'III'))
    }
  }

  if (isBold) {
    labels <- sprintf("<b>%s</b>", labels)
  }

  x <- if (isLog) {
    log10(redux$GMT)
  } else {
    redux$GMT
  }

  fig %>% add_annotations(x = x,
                          y = 1 + y_shift,
                          xref = "x",
                          yref = "y",
                          yanchor = 'bottom',
                          text = labels,
                          align = 'center',
                          showarrow = FALSE,
                          textposition = "top right",
                          textangle = '-90',
                          font = list(color = text_colour, size = font_size
                          )
  )
}

addAnnotations <- function(fig, x, labels, isLog = F, isBold = F, y_shift
  = 0, text_colour = black(),  yanchor='bottom',font_size = 16) {

  if (isBold) {
    labels <- sprintf("<b>%s</b>", labels)
  }
  if (isLog) {
    x <- log10(x)
  }

  fig %>% add_annotations(x = x,
                          y = 1 + y_shift,
                          xref = "x",
                          yref = "y",
                          yanchor = yanchor,
                          text = labels,
                          align = 'center',
                          showarrow = FALSE,
                          textposition = "top right",
                          textangle = '-90',
                          font = list(color = text_colour,
                                      size = font_size)
  )
}

addFullRegression <- function(fig = NULL, prob_df, meanColour = black(), meanName =
  'VE Model', showBands = T, meanOnLegend = F, lowOnLegend = F, lowName
                                = '95% Bands', showPredictionDots = F, showMean =
                                T, predictionColour = blue(), predictionSize = 12) {
  if (is.null(fig))
  { fig <- plot_ly() }
  if (showMean)
  { fig <- fig %>%
    add_trace(x = prob_df$GMT, y = prob_df$Mean, type = 'scatter', name = meanName,
              mode = 'lines', showlegend = meanOnLegend,
              line = list(color = meanColour, width = 3.5)) }

  if (showBands) {
    fig <- fig %>%
      add_trace(x = prob_df$GMT, y = prob_df$Mean_Low, type = 'scatter', name =
        lowName, mode = 'lines', showlegend = lowOnLegend,
                line = list(color = grey(), width = 2.5)) %>%
      add_trace(x = prob_df$GMT, y = prob_df$Mean_Up, type = 'scatter', name = '97.5%',
                mode = 'lines', showlegend = F,
                line = list(color = grey(), width = 2.5))
  }

  if (showPredictionDots) {
    fig <- fig %>%
      add_trace(x = prob_df$GMT, y = prob_df$Mean, type = 'scatter', name = 'Prediction',
                mode = 'markers', showlegend = meanOnLegend,
                marker = list(color = predictionColour, size = predictionSize, symbol = 'square'))
  }

  fig
}

addObservations <- function(fig, eff_df, colour = red(), symbol = 'circle', name =
  'VE Estimates', showlegend = T, legendgroup = NULL) {

  fig %>% add_trace(data = eff_df, x = ~GMT, y = ~Efficacy, type = 'scatter', name =
    name, mode = 'markers', showlegend = showlegend, legendgroup = legendgroup,
                    marker = list(size = 10, color = colour, symbol = symbol))

}

addVariantPredictionLabel <- function(fig, variant, y = 0.185, x = 10, font_size = 20, extra_y = 0.0015) {
  label <- paste0('<b>', toupper(greekVariantNames(variant)), '<b>')
  fig %>% add_text(x = x, y = y + extra_y, textfont = list(size = font_size), text = label, showlegend = F)
}

############ MAIN FUNCTIONS

plotEarlyTitresHistogram <- function(titres_df, pdfs = NULL, singleVariant = NULL) {

  trans <- function(x) { log10(x) }
  local <- function(variant, colour) {
    data <- columnToRealVector(titres_df[(variant)])
    gmean <- Gmean(data, method = 'classic')

    ymax <- 0.35
    xmid <- 1

    fig <- plot_ly() %>%
      add_trace(x = trans(data), type = "histogram", name = variant,
                histnorm = "probability", marker = list(color = colour), xbins = list
        (size = 0.1)) %>%
      add_trace(x = trans(gmean), y = 0, name = 'GMT', type = 'scatter', mode = 'markers',
                marker = list(color = black(), size = 15, symbol = 'square'),
                showlegend = F) %>%
      add_trace(x = trans(40), y = c(0, 0.15), name = 'GMT', type = 'scatter', mode = 'line',
                line = list(color = black(), width = 2, dash = 'dash'),
                showlegend = F) %>%
      add_text(x = trans(40), y = 0.16, name = 'GMT', text = '<b>LOWER 40<b>', type = 'text',
               showlegend = F) %>%
      add_text(x = trans(2560), y = 0.16, name = 'GMT', text = '<b>UPPER 2560<b>', type = 'text',
               showlegend = F) %>%
      add_trace(x = trans(2560), y = c(0, 0.15), name = 'GMT', type = 'scatter', mode = 'line',
                line = list(color = black(), width = 2, dash = 'dash'),
                showlegend = F) %>%

      layout(xaxis = list(showgrid = F, showticklabels = F),
             yaxis = list(range = c(-0.02, ymax), showticklabels = T, showgrid = T,
                          title = ''),
             showlegend = F) %>%
      add_annotations(x = xmid,
                      y = ymax,
                      text = greekVariantNames(variant),
                      showarrow = FALSE,
                      textposition = "top left")

    if (variant == 'WildType') {
      fig <- fig %>%
        add_annotations(x = trans(c(5, 10)),
                        y = 0.01,
                        text = c('None', 'Weak'),
                        showarrow = FALSE,
                        textposition = "top center") %>%
        add_annotations(x = trans(c(40, 5760)),
                        y = -0.01,
                        text = c('40*', '5760*'),
                        showarrow = F,
                        textposition = "top center")

    }

    if (is.null(pdfs) | !variant %in% pdfs$Variant) {
      return(fig)
    }

    row <- filter(pdfs, Variant == variant)
    fig <- fig %>%
      add_trace(
        x = unlist(row$Titres), y = unlist(row$Density), type = 'scatter', name =
          'pdf', mode =
          'lines',
        showlegend = F,
        line = list(color = black(), width = 2)
      ) %>%
      add_trace(
        x = unlist(row$Reps), y = 0, type = 'scatter', name =
          'rep', mode =
          'markers',
        showlegend = F,
        marker = list(color = red(), size = 10)
      )
    fig
  }

  fig <- if (!is.null(singleVariant)) {
    local(variant = singleVariant, colour = "#009E73") %>%
      layout(
        xaxis = list(
          title = '<b>LOG NABT<b>'),
        showlegend = F)
  }
  else {
    subplot(
      local(variant = 'D614G', colour = "#009E73"),
      local(variant = 'B.1.1.7', colour = "#E69F00"), #orange
      local(variant = 'B.1.351', colour = blue()),
      local(variant = 'B.1.617.2', colour = red()),
      nrows = 2, shareX = T) %>%
      layout(
        xaxis = list(
          title = '<b>LOG NABT<b>'),
        showlegend = F) }

  showInChrome(fig)
  # fig
}

plotKidneyDiabetesSubgroups <- function(df, healthy) {

  ve_font <- 16
  local <- function(fig, df, observed_ys) {
    shift_y <- 0.025
    disease <- df$Disease[1]
    label <- sprintf("<b>%s</b>", toupper(disease))
    fig <- fig %>% add_text(x = 20, y = max(observed_ys), text = label, size = 16,
                            textposition = 'middle right', showlegend = F)

    show_all_legend <- disease == 'Diabetes'
    for (i in seq_len(nrow(df))) {
      show_legend <- (disease == 'Diabetes' & i == 1)
      y <- observed_ys[i]
      row <- df[i,]
      all_eff <- filter(healthy, Comment.source == row$Comment.source)
      if (!is.null(all_eff) & nrow(all_eff) == 1) {

        fig <- add_trace(fig, x = all_eff$Efficacy, y = y + shift_y / 2, type =
          'scatter',
                         mode =
                           'markers',
                         marker = list(color = blueLight(), size = 15), showlegend =
                           show_all_legend, name
                           = '<b>VE HEALTHY<b>',
                         error_x = errorBar(
                           all_eff$LCL,
                           all_eff$Efficacy,
                           all_eff$UCL,
                           colour = blueLight()
                         )) %>%
          add_text(x = all_eff$Efficacy, y = y + shift_y / 2 + 0.01, text = round
          (all_eff$Efficacy, 0), textfont = list(size = ve_font),
                   showlegend = F)
        show_all_legend <- F
      }

      fig <- add_trace(fig, x = row$Efficacy, y = y, type = 'scatter', mode =
        'markers',
                       marker = list(color = red(), size = 15), showlegend = show_legend, name
                         = '<b>VE DISEASE<b>',
                       error_x = errorBar(
                         row$LCL,
                         row$Efficacy,
                         row$UCL,
                         colour = red()
                       )) %>%
        add_text(x = row$Efficacy, y = y + 0.01, text = round(row$Efficacy, 0), textfont = list(size = ve_font),
                 showlegend = F) %>%
        add_text(x = row$LCL - 3, y = y - 0.0005, text =
          countryAbbrevs(row$Country), textfont = list(size = ve_font),
                 showlegend = F) %>%
        add_trace(x = row$Predicted_Mean, y = y - shift_y, type = 'scatter', mode =
          'markers',
                  marker = list(color = black(), size = 15, symbol = 'square'),
                  showlegend = (disease == 'Diabetes' & i == 2), name
                    = '<b>PREDICT. DISEASE<b>',
                  error_x = errorBar(
                    row$Predicted_LCL,
                    row$Predicted_Mean,
                    row$Predicted_UCL)) %>%
        add_text(x = row$Predicted_Mean, y = y + 0.011 - shift_y, text = round
        (row$Predicted_Mean), textfont = list(size = ve_font), showlegend = F)
    }

    fig
  }

  df <- df[order(df$Efficacy),]
  diab <- filter(df, Disease == 'Diabetes')
  kidney <- filter(df, Disease == 'Kidney')
  observed_ys_diab <- seq(from = 1, by = -0.09, length.out = nrow(diab))
  observed_ys_kidney <- seq(from = min(observed_ys_diab) - 0.12, by =
    -0.09, length.out = nrow(kidney))

  fig <- plot_ly() %>%
    local(diab, observed_ys = observed_ys_diab) %>%
    local(kidney, observed_ys = observed_ys_kidney) %>%
    layout(
      yaxis = list(showticklabels = FALSE, showgrid = F, zeroline = T),
      xaxis = list(title = '<b>VE<b>',
                   titlefont = list(size = ve_font + 2),
                   tickfont = list(size = ve_font),
                   range = c(0, 100),
                   showgrid = F, zeroline = F))

  showInChrome(fig)
}

plotFittedModelLog <- function(prob_df, efficacies_df, title = '') {

  fig <- addFullRegression(prob_df = prob_df,
                           meanName = 'Comirnaty Model', meanOnLegend = T, lowOnLegend = T,
                           showPredictionDots = T) %>%
    addObservations(efficacies_df) %>%
    addVariantAnnotations(efficacies_df, isLog = T, showDose = T, isBold = T) %>%
    layout(xaxis = list(title = '<b>LOG GMT<b>', titlefont = list(size = 22), tickfont = list(size = 18),
                        showline = T, range = c(1, log10
      (420)),
                        showgrid = F, type = 'log',
                        tickvals = round(prob_df$GMT, 0)),
           yaxis = list(title = '<b>VE %<b>', tickfont = list(size = 18),
                        titlefont = list(size = 22), range
                          = c(0, 100)), title = toupper
      (title),
           showlegend = T, legend = list(x = 0.78, y = 0.3, font = list(size = 22)))

  showInChrome(fig)
}

#Plots vaxzevria
plotVaxzevriaModel <- function(prob_df, efficacies_df, extra_model_df, title = '') {

  diffs <- extra_model_df$Mean - prob_df$Mean
  fig <- addFullRegression(prob_df = extra_model_df, showBands = F, meanColour = faintBlue(),
                           meanName = 'Comirnaty Model', meanOnLegend = T) %>%
    addFullRegression(prob_df = prob_df, meanName = 'Vaxzevria Model', meanOnLegend = T,
                      lowName = 'Vaxzevria 95% Bands', lowOnLegend = T) %>%
    add_text(x = prob_df$GMT, y = prob_df$Mean_Low - 2, text = sprintf
    ("<b>%s</b>", round(diffs, 1)), textfont = list(size = 15), showlegend = F,
             name =
               'Vaxzevria -
              Comirnaty') %>%
    addObservations(efficacies_df, name = 'Vaxzevria VE Estimates') %>%
    addVariantAnnotations(efficacies_df, isLog = T, showDose = T, isBold = T) %>%

    layout(xaxis = list(title = '<b>LOG GMT<b>', titlefont = list(size = 22), tickfont = list(size = 18), showline = T, range = log10(c(18,
                                                                                                                                        120)),
                        showgrid = F, type = 'log',
                        tickvals = round(prob_df$GMT, 0)),
           yaxis = list(title = '<b>VE %<b>', tickfont = list(size = 18),
                        titlefont = list(size = 22), range
                          = c(0, 100)), title = toupper
      (title),
           showlegend = T, legend = list(x = 0.78, y = 0.3, font = list(size = 22)))

  showInChrome(fig)
}

plotBothFittedModels <- function(prob_df, efficacies_df, title = '') {
  eff_az <- filter(efficacies_df, Vaccine == 'Vaxzevria')
  eff_pf <- filter(efficacies_df, Vaccine != 'Vaxzevria')

  fig <- addFullRegression(prob_df = prob_df, meanName = 'Combined Model', meanOnLegend = T, lowOnLegend = T) %>%
    addObservations(eff_df = eff_az, colour = black(),
                    name = 'Vaxzevria VE Estimates') %>%
    addObservations(eff_df = eff_pf, colour = blue(),
                    name = 'Comirnaty VE Estimates') %>%
    addVariantAnnotations(eff_az, isLog = T, showDose = T, isBold = T) %>%
    addVariantAnnotations(eff_pf, isLog = T, showDose = T, isBold = T, text_colour = blue()) %>%
    layout(xaxis = list(title = '<b>LOG GMT<b>', titlefont = list(size = 22), tickfont = list(size = 18),
                        showline = T, range = log10(c(15, 420)),
                        showgrid = F, type = 'log',
                        tickvals = round(prob_df$GMT, 0)),
           yaxis = list(title = '<b>VE %<b>', tickfont = list(size = 18),
                        titlefont = list(size = 22), range
                          = c(0, 100)), title = toupper
      (title),
           showlegend = T, legend = list(x = 0.78, y = 0.3, font = list(size = 22)))
  showInChrome(fig)
}

plotPredictionConcordance <- function(both, rho) {

  label <- figureCI(rho$est, rho$lwr.ci, rho$upr.ci, 2, start = '\U03C1: ')
  fig <- plot_ly() %>%
    add_trace(x = c(30, 95), y = c(30, 95), mode = 'lines',  type = 'scatter', showlegend
      = F, line = list(color = greyMedium(), width = 3)) %>%
    add_trace(data = filter(both, Vaccine == 'Vaxzevria', Dose == 1), x =~Mean, y = ~Mean_VE, name = 'Vaxzevria I',
              type = 'scatter', mode = 'markers',
              marker = list(size = 20, color = red())) %>%
    add_trace(data = filter(both, Vaccine == 'Vaxzevria', Dose == 2), x =~Mean, y = ~Mean_VE, name = 'Vaxzevria II',
              type = 'scatter', mode = 'markers',
              marker = list(size = 20, color = red(), symbol = 'square')) %>%
    add_trace(data = filter(both, Vaccine == 'Comirnaty', Dose == 1), x =~Mean, y = ~Mean_VE, name = 'Comirnaty I',
              type = 'scatter',
              mode = 'markers', marker = list(size = 20, color = blueLight())) %>%
    add_trace(data = filter(both, Vaccine == 'Comirnaty', Dose == 2), x =~Mean, y = ~Mean_VE, name = 'Comirnaty II',
              type = 'scatter',
              mode = 'markers', marker = list(size = 20, color = blueLight(),
                                              symbol = 'square')) %>%
    add_text(x = 40, y = 95, text = label, textfont = list(size = 22),
             showlegend = F) %>%
    layout(xaxis = list(tickfont = list(size = 18),
                        titlefont = list(size = 22), title =
                          '<b>PREDICTED VE %<b>', showline = T,
                        showgrid = F),
           yaxis = list(tickfont = list(size = 18),
                        titlefont = list(size = 22), title = '<b>MEAN ESTIMATED VE %<b>', showline = T,
                        showgrid = F),
           legend = list(x = 0.78, y = 0.3, font = list(size = 22))
    )

  showInChrome(fig)
}


plotSimpleFullModelsLog <- function(prob_df, efficacies_df, simple_fun) {

  log_titres <- log10(seq(10, 420))
  simple_y <- 100 * simple_fun(log_titres)

  fig <- addFullRegression(prob_df = prob_df,
                           lowOnLegend = T, lowName = 'Full bands',
                           showMean = F) %>%
    add_trace(x = 10^log_titres, y = simple_y, type = 'scatter',
              name = 'Simple Model',
              mode = 'lines', showlegend = T,
              line = list(color = black(), width = 3.5)) %>%
    addObservations(efficacies_df) %>%
    add_trace(x = prob_df$GMT, y = prob_df$Mean, type = 'scatter', name =
      'Full Predictions',
              mode = 'markers', showlegend = T,
              marker = list(color = blue(), size = 12, symbol = 'square')) %>%
    addVariantAnnotations(efficacies_df, isLog = T, showDose = T, isBold = T) %>%
    layout(xaxis = list(title = '<b>LOG GMT<b>', titlefont = list(size = 22), tickfont = list(size = 18),
                        showline = T, range = c(1, log10
      (420)),
                        showgrid = F, type = 'log',
                        tickvals = round(prob_df$GMT, 0)),
           yaxis = list(title = '<b>VE %<b>', tickfont = list(size = 18),
                        titlefont = list(size = 22), range
                          = c(0, 100)),
           showlegend = T, legend = list(x = 0.78, y = 0.3, font = list(size = 22)))

  showInChrome(fig)
}

plotAgeSubgroups <- function(prob_df, age_df) {

  labels <- NULL
  for (k in seq_len(nrow(age_df))) {
    row <- age_df[k,]
    labels <- c(labels, paste0(row$Age_Start, '-', row$Age_End))
  }

  y_shift_pred <- c(-2.5, 2, 2, 2, -3.5, -2.5)
  show_legend <- T
  fig <- plot_ly(showlegend = show_legend) %>%
    addFullRegression(prob_df = prob_df,
                      meanName = 'Comirnaty Model',
                      meanOnLegend = show_legend, lowOnLegend = show_legend) %>%
    add_trace(x = age_df$GMT, y = age_df$Mean, type = 'scatter', name =
      'Prediction', mode = 'markers',
              marker = list(size = 15, color = blueLight(), symbol = 'square')) %>%
    add_text(x = age_df$GMT, y = age_df$Mean + y_shift_pred, text = sprintf
    ("<b>%s</b>",
     round(age_df$Mean, 1)), textfont = list(size = 15), showlegend = F) %>%
    addObservations(age_df, name = 'VE Estimates', showlegend =
      show_legend) %>%
    addAnnotations(x = age_df$GMT, labels = labels, isLog = T, isBold = T, y_shift = age_df$Mean_Low) %>%
    add_text(x = age_df$GMT, y = age_df$Mean_Up + 1.15, textfont = list(size = 15), text = sprintf("<b>%s</b>",
                                                                                                   round(age_df$Error, 1)),
             textfont = list(size = 15), showlegend = F, textfont = list(color = black())) %>%
    addVariantAnnotations(age_df, isLog = T, distinct_only = F, showDose = T, isBold = T) %>%
    layout(xaxis = list(title = '<b>LOG GMT<b>', titlefont = list(size = 22), tickfont = list(size = 18),
                        showline = T, range = log10(c(20, 275)),
                        showgrid = F, type = 'log',
                        tickvals = round(age_df$GMT, 0)),
           yaxis = list(title = '<b>VE %<b>', tickfont = list(size = 18),
                        titlefont = list(size = 22), range
                          = c(0, 102)),
           showlegend = T, legend = list(x = 0.78, y = 0.3, font = list(size = 22)))

  showInChrome(fig)
}

plotComirnatyEpochs <- function(observed, intervals, variants, dose) {

  observed <- filter(observed,
                     Dose == dose)
  intervals <- filter(intervals, Dose == dose)

  local <- function(variant, full_show = T) {
    observed_one <- observed %>% filter(Variant == variant)
    observed_one <- observed_one[order(observed_one$Efficacy),]
    observed_one$Country <-  countryAbbrevs(observed_one$Country)

    interval_one <- intervals %>%
      filter(Variant == variant) %>%
      distinct(Epoch, Mean, .keep_all = T)
    top_y <- 0.185

    predict_ys <- seq(from = top_y - 0.03, by = -0.015, length.out = nrow(interval_one))
    observed_shift_y <- if (variant != 'B.1.617.2') { min(predict_ys) - 0.045 } else {
      min(predict_ys) - 0.03
    }

    observed_ys <- seq(from = observed_shift_y, by = -0.01, length.out = nrow(observed_one))

    fig <- plot_ly(showlegend = F) %>%
      add_trace(x = observed_one$Mean[1], y = c(0, top_y), name = 'Mean VE Estimate', type = 'scatter', mode = 'lines',
                line = list(color = greyMedium(), width = 2, dash = 'dot'), showlegend = F) %>%
      add_text(x = observed_one$Mean[1], y = top_y + 0.003, textfont = list(size = 15),
               text = sprintf("<b>%s</b>", round(observed_one$Mean[1], 1)), showlegend = F) %>%
      addVariantPredictionLabel(variant = variant)

    delta_count <- 0
    for (k in seq_len(length(predict_ys))) {
      row <- interval_one[k,]
      mark_size <- 15
      colour <- black()
      symbol <- 'square'
      ep <- 'I'

      if (row$Epoch == 2) {
        colour <- green()
        ep <- 'II'
      }
      if (row$Epoch == 3) {
        colour <- blue()
        ep <- 'III'
      }
      if (row$Name == 'CROMER') {
        colour <- black()
        symbol <- 'triangle-up'
        mark_size <- 22
        ep <- 'CRO'
      }

      fig <- add_trace(fig, x = row$Mean, y = predict_ys[k], type = 'scatter', mode =
        'markers',
                       marker = list(color = colour, size = mark_size, symbol = symbol),
                       showlegend = F,
                       error_x = errorBar(
                         row$Mean_Low,
                         row$Mean,
                         row$Mean_Up
                       )) %>%
        add_text(x = row$Mean_Up + 3, y = predict_ys[k], textfont = list(size = 15),
                 text = sprintf("<b>%s</b>", round(row$Mean, 1)),
                 showlegend = F) %>%
        add_text(x = row$Mean_Low - 1, y = predict_ys[k], textfont = list(size = 15),
                 text = paste0('<b>', ep, '<b>'), textposition = 'middle left',
                 showlegend = F)

      if (row$Variant == 'B.1.617.2' & row$Epoch == 2) {
        delta_count <- delta_count + 1 }
    }

    for (k in seq_len(length(observed_ys))) {
      row <- observed_one[k,]
      country_x <- if (row$LCL - 3 > 0) {
        row$LCL - 3
      } else { row$UCL + 3 }

      fig <- add_trace(fig, x = row$Efficacy, y = observed_ys[k], type = 'scatter', mode = 'markers',
                       marker = list(color = red(), size = 15), showlegend = F, name = 'VE Estimates',
                       error_x = errorBar(
                         row$LCL,
                         row$Efficacy,
                         row$UCL
                       )) %>%
        add_text(x = row$UCL + 3, y = observed_ys[k], textfont = list(size = 15),
                 text = sprintf("<b>%s</b>", round(row$Efficacy, 1)), showlegend = F) %>%
        add_text(x = country_x, y = observed_ys[k] - 0.0005, textfont = list(size = 15), text = paste0('<b>', row$Country, '<b>'), showlegend = F)
    }

    x_title <- if (full_show) {
      '<b>VE %<b>'
    }else { '' }

    #range=c(-0.05,0.18)
    fig %>% layout(
      yaxis = list(showticklabels = FALSE, showgrid = F, zeroline = F),
      xaxis = list(title = x_title,
                   range = c(0, 110), tickfont = list(size = 18),
                   titlefont = list(size = 22),
                   showgrid = F, zeroline = T),
      showlegend = F)
  }

  count <- length(variants)
  if (count == 1) {
    fig <- local(variant = variants[1])
  }
  if (count == 2) {
    fig <- subplot(local(variant = variants[1]),
                   local(variant = variants[2]),
                   nrows = 1, shareX = T, shareY = T)
  }
  if (count == 3) {
    fig <- subplot(local(variant = variants[1]),
                   local(variant = variants[2]),
                   local(variant = variants[3]),
                   nrows = 1, shareX = T, shareY = T)
  }
  if (count == 4) {
       variants <-   c('B.1.1.7', 'B.1.617.2','D614G',  'B.1.351') ######## re-arrange the order to suit the paper.
    fig <- subplot(local(variant = variants[1]),
                   local(variant = variants[2], full_show = F),
                   local(variant = variants[3]),
                   local(variant = variants[4], full_show = F),
                   nrows = 2, shareX = T, shareY = T)
  }

  showInChrome(fig %>% layout(showlegend = F))
}

plotVariantPredictionsPerDose <- function(observed, intervals, variants, dose, show_axis = T) {

  local <- function(variant, full_show = T) {
    observed_one <- observed %>% filter( Variant == variant, Dose == dose)
    observed_one <- observed_one[order(observed_one$Efficacy),]
    observed_one$Country <-  countryAbbrevs(observed_one$Country)
    interval_one <- intervals %>% filter( Variant == variant, Dose == dose)

    top_y <- 0.185
    predict_ys <- seq(from = top_y - 0.03, by = -0.015, length.out = nrow(interval_one))
    observed_ys <- seq(from = min(predict_ys) - 0.03, by = -0.01, length.out = nrow(observed_one))

    fig <- plot_ly(showlegend = F) %>%
      add_trace(x = observed_one$Mean[1], y = c(0, top_y), name = 'Mean VE Estimate', type = 'scatter', mode = 'lines',
                line = list(color = greyMedium(), width = 2, dash = 'dot'), showlegend = F) %>%
      add_text(x = observed_one$Mean[1], y = top_y + 0.003, textfont = list(size = 15),
               text = sprintf("<b>%s</b>", round(observed_one$Mean[1], 1)), showlegend = F) %>%
      addVariantPredictionLabel(variant = variant)

    for (k in seq_len(length(predict_ys))) {
      row <- interval_one[k,]
      prediction_label <- 'PRE'
      colour <- blue()
      symbol <- 'square'
      mark_size <- 15

      if (row$Name == 'CROMER') {
        prediction_label <- 'CRO'
        colour <- black()
        mark_size <- 22
        symbol <- 'triangle-up'
      }

      if (row$Name %in% c('Vaxzevria_from_Comirnaty', 'Comirnaty_from_Vaxzevria')) {
        prediction_label <- 'IND'
        colour <- green()
      }

      fig <- add_trace(fig, x = row$Mean, y = predict_ys[k], type = 'scatter', mode =
        'markers',
                       marker = list(color = colour, size = mark_size, symbol = symbol), showlegend = F,
                       error_x = errorBar(
                         row$Mean_Low,
                         row$Mean,
                         row$Mean_Up
                       )) %>%
        add_text(x = row$Mean_Up + 3, y = predict_ys[k], textfont = list(size = 15),
                 text = sprintf("<b>%s</b>", round(row$Mean, 1)),
                 showlegend = F) %>%
        add_text(x = row$Mean_Low - 1, y = predict_ys[k], textfont = list(size = 15),
                 text = paste0('<b>', prediction_label, '<b>'), textposition = 'middle left',
                 showlegend = F)
    }


    for (k in seq_len(length(observed_ys))) {
      row <- observed_one[k,]
      country_x <- if (row$LCL - 3 > 0) {
        row$LCL - 3
      } else { row$UCL + 3 }

      fig <- add_trace(fig, x = row$Efficacy, y = observed_ys[k], type = 'scatter', mode = 'markers',
                       marker = list(color = red(), size = 15), showlegend = F, name = 'VE Estimates',
                       error_x = errorBar(
                         row$LCL,
                         row$Efficacy,
                         row$UCL
                       )) %>%
        add_text(x = row$UCL + 3, y = observed_ys[k], textfont = list(size = 15),
                 text = sprintf("<b>%s</b>", round(row$Efficacy, 1)), showlegend = F) %>%
        add_text(x = country_x, y = observed_ys[k] - 0.0005, textfont = list(size = 15), text = paste0('<b>', row$Country, '<b>'),
                 showlegend = F)
    }

    x_title <- if (full_show & show_axis) {
      '<b>VE %<b>'
    }else { '' }

    fig %>% layout(
      yaxis = list(showticklabels = FALSE, showgrid = F, zeroline = F),
      xaxis = list(title = x_title,
                   range = c(0, 110), tickfont = list(size = 18),
                   titlefont = list(size = 22), showticklabels = full_show & show_axis,
                   showgrid = F, zeroline = T),
      showlegend = F)
  }

  count <- length(variants)
  if (count == 1) {
    fig <- local(variant = variants[1])
  }
  if (count == 2) {
    fig <- subplot(local(variant = variants[1]),
                   local(variant = variants[2]),
                   nrows = 1, shareX = T, shareY = T)
  }
  if (count == 3) {
    fig <- subplot(local(variant = variants[1]),
                   local(variant = variants[2]),
                   local(variant = variants[3]),
                   nrows = 1, shareX = T, shareY = T)
  }
  if (count == 4) {
     variants <-   c('B.1.1.7', 'B.1.617.2','D614G',  'B.1.351') ######## re-arrange the order to suit the paper.
    fig <- subplot(local(variant = variants[1]),
                   local(variant = variants[2]),
                   local(variant = variants[3]),
                   local(variant = variants[4]),
                   nrows = 1, shareX = T, shareY = T)
  }

  fig <- fig %>%
    layout(legend = list(x = 0.002, y = 0.9))

  showInChrome(fig)
}


plotDeltaOmicronFitted <- function(prob_delta, eff_delta, prob_omi, eff_omi, hd) {

  omi_cohorts <- distinct(eff_omi, Cohort, .keep_all = T)
  omi_cohorts <- omi_cohorts[(order(omi_cohorts$Cohort)),]
  delta_cohorts <- distinct(eff_delta, Cohort, .keep_all = T)
  delta_cohorts <- delta_cohorts[(order(delta_cohorts$Cohort)),]

  hd_omi <- filter(hd, Variant == 'B.1.1.529')
  hd_delta <- filter(hd, Variant == 'B.1.617.2')

  omi_x <- c(omi_cohorts$GMT, hd_omi$GMT)
  delta_x <- c(delta_cohorts$GMT, hd_delta$GMT)
  labels <- c('BOOST III', 'EARLY II', 'LATE II', 'HD BOOST III')

  delta_tix <- c(prob_delta$GMT, hd_delta$GMT)
  omi_tix <- c(prob_omi$GMT, hd_omi$GMT)
  x_range <- log10(c(22, 730))
  x2 <- list(
    tickfont = list(color = blue(), size = 18),
    titlefont = list(color = blue(), size = 22),
    tickvals = delta_tix,
    ticktext = round(delta_tix, 0),
    overlaying = "x",
    side = "top",
    anchor = "y",
    showline = T,
    showgrid = F,
    range = x_range,
    type = 'log',
    title = "<b>DELTA: LOG GMT <b>")

  delta <- rbind(prob_delta, hd_delta)
  delta <- delta[(order(delta$GMT)),]
  delta_preds <- unique(delta$Mean)
  omi <- rbind(prob_omi, hd_omi)
  omi <- omi[(order(omi$GMT)),]
  omi_preds <- unique(omi$Mean)

  fig <- addFullRegression(prob_df = delta,  meanName = 'Delta Model', meanOnLegend = T, lowOnLegend = T,
                           showPredictionDots = T, meanColour = blue(), predictionColour = blue(), predictionSize = 15) %>%
    addObservations(eff_delta) %>%
    add_text(x = unique(delta$GMT) + c(6, 12, 20, 28), y = delta_preds - 1.5, xaxis = 'x2', text = sprintf("<b>%s</b>",
                                                                                                           round(delta_preds, 1)), textfont = list(size = 17)) %>%
    add_text(x = 139, y = 85, text = '<b>DELTA  FIT<b>', textfont = list(size = 19, color = blue())) %>%
    addAnnotations(x = delta_x, labels = labels, isLog = T, isBold = T, text_colour = blue(), y_shift = 107,  font_size = 18, yanchor = 'top') %>%
    addFullRegression(prob_df = omi,  meanName = 'Omicron Model', meanOnLegend = T, lowOnLegend = T,
                      showPredictionDots = T, predictionColour = black, predictionSize = 15) %>%
    addObservations(eff_omi, colour = green()) %>%
    add_text(x = unique(omi$GMT), y = omi_preds - 2, text = sprintf("<b>%s</b>",
                                                                    round(omi_preds, 1)), textfont = list(size = 17)) %>%
    add_text(x = 26, y = 47, text = '<b>OMICRON  FIT<b>', textfont = list(size = 19, color = black)) %>%
    addAnnotations(x = omi_x, labels = labels, isLog = T,  font_size = 18, isBold = T) %>%
    layout(xaxis = list(title = '<b>OMICRON: LOG GMT<b>', tickfont = list(size = 18),
                        titlefont = list(size = 22), showline = T,
                        showgrid = F, type = 'log',
                        tickvals = omi_tix,
                        range = x_range,
                        ticktext = round(omi_tix, 0)),
           xaxis2 = x2,
           yaxis = list(title = '<b>VE %<b>', range = c(0, 109), tickfont = list(size = 18),
                        titlefont = list(size = 22), tickvals = seq(from = 0, to = 100, by = 20)),
           margin = list(t = 55),
           showlegend = F)

  showInChrome(fig)
}

plotDeltaOmicronScaledBoth <- function(prob_df, both) {
  both <- both[(order(both$GMT)),]
  omi <- filter(both, Variant == 'B.1.1.529')
  delta <- filter(both, Variant == 'B.1.617.2')

  expected <- omi %>% distinct(Cohort, .keep_all = T)
  expected_delta <- delta %>% distinct(Cohort, .keep_all = T)

  pred_shift <- c(1.3, 3.5, 12)
  pred_shift_delta <- c(4, 8.5, 20) * 1.3
  delta_tix <- expected_delta$GMT
  x_range <- log10(c(20, 430))
  x2 <- list(
    tickfont = list(color = blue(), size = 18),
    titlefont = list(color = blue(), size = 22),
    tickvals = delta_tix,
    ticktext = round(delta_tix, 0),
    overlaying = "x",
    side = "top",
    anchor = "y",
    showline = T,
    showgrid = F,
    range = x_range,
    type = 'log',
    title = "<b>DELTA: LOG GMT (Primary Scale)<b>")

  labels <- c('LATE II', 'EARLY II', 'BOOST III')

  fig <- addFullRegression(prob_df = prob_df,  meanOnLegend = T, meanName = 'Primary Model', predictionSize = 17) %>%
    addAnnotations(x = expected$GMT, labels = labels,  font_size = 18,
                   isBold = T, isLog = T, y_shift = 0) %>%
    addAnnotations(x = delta_tix, labels = labels, isLog = T, isBold = T,  font_size = 18,text_colour = blue(), y_shift = 107, yanchor = 'top') %>%
    addObservations(omi, colour = green(), name =
      'VE Estimates', legendgroup = 'omi') %>%
    addObservations(delta, colour = red(), name =
      'VE Estimates', showlegend = T, legendgroup = 'delta') %>%
    add_text(x = expected$GMT - pred_shift, y = expected$Efficacy_Mean,
             text = sprintf("<b>%s</b>", round(expected$Efficacy_Mean, 1)), textfont = list(size = 17), showlegend = F) %>%
    add_text(x = expected_delta$GMT - pred_shift_delta, y = expected_delta$Efficacy_Mean, xaxis = 'x2',
             text = sprintf("<b>%s</b>", round(expected_delta$Efficacy_Mean, 1)), textfont = list(size = 17), showlegend = F) %>%
    add_trace(x = expected$GMT, y = expected$Mean, type = 'scatter', name =
      'Prediction', mode = 'markers', legendgroup = "omi",
              legendgrouptitle = list(text = "Omicron", font = list(size = 22)),
              marker = list(size = 17, color = black, symbol = 'square')) %>%
    add_trace(x = expected_delta$GMT, y = expected_delta$Mean, type = 'scatter', name =
      'Prediction', mode = 'markers', legendgroup = "delta",
              legendgrouptitle = list(text = "Delta", font = list(size = 22)),
              marker = list(size = 17, color = blue(), symbol = 'square')) %>%
    add_text(x = expected_delta$GMT, y = expected_delta$Mean + c(1.6, -2, -2),
             text = sprintf("<b>%s</b>", round(expected_delta$Mean, 1)), textfont = list(size = 17), showlegend = F) %>%
    add_text(x = expected$GMT, y = expected$Mean - 2,
             text = sprintf("<b>%s</b>", round(expected$Mean, 1)), textfont = list(size = 17), showlegend = F) %>%
    addObservations(delta, colour = red(), name =
      'VE Delta', symbol = 'circle', showlegend = F) %>%
    add_trace(x = expected$GMT, y = expected$Efficacy_Mean, type =
      'scatter', name =
                'Mean VE', mode = 'markers', showlegend = F,
              marker = list(size = 20, color = black, symbol = '141', line = list
              (width = 3))) %>%
    add_trace(x = expected_delta$GMT, y = expected_delta$Efficacy_Mean, type =
      'scatter', name =
                'Mean VE', mode = 'markers', showlegend = F,
              marker = list(size = 20, color = black, symbol = '141', line = list
              (width = 3))) %>%
    layout(xaxis = list(title = '<b>OMICRON: LOG GMT (Primary Scale)<b>', tickfont = list(size = 18),
                        titlefont = list(size = 22), showline = T,
                        showgrid = F, type = 'log', range = x_range,
                        ticktext = round(expected$GMT, 0),
                        tickvals = expected$GMT),
           xaxis2 = x2,
           yaxis = list(title = '<b>VE %<b>', range = c(0, 109), tickfont = list(size = 18),
                        titlefont = list(size = 22)),
           margin = list(t = 55),
           legend = list(x = 0.83, y = 0.27, font = list(size = 22)))

  showInChrome(fig)
}

##### NB?: add optional plot titles here and in Main.
