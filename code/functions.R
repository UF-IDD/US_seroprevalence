###############################################################################
# Functions used throughout the scripts.
# ~~~
# 
###############################################################################


#' Get root-mean-square error.
#'
#' @param dat :numeric vector: Data to which the model was fit to.
#' @param pred :numeric vector: Model predictions.
#'
#' @return the RMSE.

RMSE = function(dat, pred)
{
  if (length(dat) != length(pred))
  {
    stop('Error in function RMSE : length of vectors not the same')
  }

  return(sqrt(mean((dat - pred)^2, na.rm = TRUE)))
}



#' Adjust number of cases to account for waning using a step function.
#'
#' @description Function used by `GetAdjustedCases` below to produce time series
#' of cases adjusted for waning.
#'
#' @param x :numeric vector: number of weekly cases (non-cumulative) for a
#'        state.
#' @param detect :numeric: number of weeks lead or lag between a case being
#'        reported and becoming seropositive.
#' @param ab_dur :named numeric vector: number of weeks after testing positive
#'        that a case turns negative again (seroreversion), per assay. The items
#'        of the vector are assumed to have names 'abbott', 'ortho', and
#'        'roche'.
#' @param props_* :numeric vectors: what proportion of tests used each assay in
#'        turn, per round. E.g., `props_abbott` gives the proportion of tests
#'        that used the Abbott assay, for each round. For each round, the
#'        `props_*` add to one.
#'
#' @return a matrix, with a number of time series equivalent to the number of
#' rounds. The rows are the time points, and there is a column per survey round,
#' with what the time series would be with the proportions of assays used in
#' that specific round. The "real" adjusted time series is then the diagonal of
#' the subset of the matrix only including time points for which there is a
#' serosurvey.

GetWanedCases = function(x, detect, ab_dur, 
                         props_abbott = 1, 
                         props_ortho = 1, 
                         props_roche = 1)
{
  # Add the lead or lag between reporting and seroconversion.
  x = data.table::shift(x = x, n = detect, fill = 0)
  cases_decay = matrix(data = 0, nrow = length(x), ncol = length(props_abbott))

  for (k in seq_along(props_abbott))
  {
    for (i in seq_along(x))
    {
      # Deal with situations when the times to seroreversion are longer than the
      # remaining time points.

      # Abbott assay.
      if ((i + ab_dur['abbott'] - 1) < length(x))
      {
        times = i:(i + ab_dur['abbott'] - 1)
        cases_decay[times, k] = cases_decay[times, k] + rep(x[i] * props_abbott[k], length(times))
      } else
      {
        times = i:length(x)
        cases_decay[times, k] = cases_decay[times, k] + rep(x[i] * props_abbott[k], length(times))
      }

      # Ortho assay.
      if ((i + ab_dur['ortho'] - 1) < length(x))
      {
        times = i:(i + ab_dur['ortho'] - 1)
        cases_decay[times, k] = cases_decay[times, k] + rep(x[i] * props_ortho[k], length(times))
      } else
      {
        times = i:length(x)
        cases_decay[times, k] = cases_decay[times, k] + rep(x[i] * props_ortho[k], length(times))
      }

      # Roche assay.
      if ((i + ab_dur['roche'] - 1) < length(x))
      {
        times = i:(i + ab_dur['roche'] - 1)
        cases_decay[times, k] = cases_decay[times, k] + rep(x[i] * props_roche[k], length(times))
      } else
      {
        times = i:length(x)
        cases_decay[times, k] = cases_decay[times, k] + rep(x[i] * props_roche[k], length(times))
      }
    }
  }

  cases_decay = data.frame(cases_decay)
  colnames(cases_decay) = paste0('cases_waning_', seq_len(ncol(cases_decay)))

  return(cases_decay)
}



#' Get the number of cases adjusted to account for waning.
#'
#' See `GetWanedCases` above.
#'
#' @param dat :data.table: data with the numbers of weekly cases.
#' @param durs :named numeric vector: vector with the numbers of weeks to
#'        seroreversion per assay, and the lead/lag between a case being
#'        reported and it seroconverting. Names for the elements are assumed to
#'        be: "abbott", "ortho", "roche", "lag".
#'
#' @return A version of `dat` with an additional column with the adjusted number
#' of cases.

#' @description The number of cases are adjusted to account for waning by using
#' step functions (one for each assay used) which indicate the probability of a
#' case being reported (being either zero or one). The rationale is that once a
#' case seroconverts (accounting for lead/lag between seroconversion and
#' reporting), it contributes to the cumulative reported cases only for a number
#' of weeks. 

GetAdjustedCases = function(dat, durs)
{
  dat = dat[, cases_w_cum_d := data.table::shift(cases_cumulative, n = durs[['lag']], 
                                                 fill = 0), 
            by = 'state']

  s = copy(dat)
  s = s[!is.na(seroprevalence), ]
  s$cases_waning = as.numeric(NA)

  st = unique(s$state)

  for (i in seq_along(st))
  {
    w_st = which(s$state == st[i])
    state = subset(dat, state == st[i])

    # This produces a matrix (see function above).
    c_waning = state[, GetWanedCases(cases, detect = durs['lag'], 
                                     ab_dur = durs[c('abbott', 'ortho', 'roche')],
                                     props_abbott = Abbott_Architect_perc[!is.na(Abbott_Architect_perc)] / 100,
                                     props_ortho = Ortho_VITROS_perc[!is.na(Ortho_VITROS_perc)] / 100,
                                     props_roche = Roche_Elecsys_perc[!is.na(Roche_Elecsys_perc)] / 100)]
    w_r = which(!is.na(state$seroprevalence))
    # Get time series with the adjusted number of cases as per the proportions
    # of assays used in each survey round.
    c_waning = diag(as.matrix(c_waning[w_r, ]))

    s = s[w_st, cases_waning := c_waning]
  }

  s = s[, `:=`(cases_perc_waning = cases_waning / state_population * 100)]

  return(s)
}



#' Simplified version to get waned cases, for plotting only.
#'
#' @param x :numeric vector: the number of daily cases (non-cumulative).
#' @param ab_dur :numeric: number of days after which cases test negative.
#'
#' @return cumulative number of cases, accounting for waning after `wane_days`
#' number of days.

GetWanedCasesSimple = function(x, detect, ab_dur)
{
  x = data.table::shift(x, detect, fill = 0)
  cases_decay = vector('numeric', length = length(x))

  for (i in seq_along(x))
  {
    if ((i + ab_dur - 1) < length(x))
    {
      times = i:(i + ab_dur - 1)
      cases_decay[times] = cases_decay[times] + rep(x[i], length(times))
    } else
    {
      times = i:length(x)
      cases_decay[times] = cases_decay[times] + rep(x[i], length(times))
    }
  }

  return(cases_decay)
}



#' Pad survey:::predict.svyglm output with NAs.
#'
#' This is necessary if you want survey:::predict.svyglm to behave as if `na.act
#' = na.exclude` had been passed. Currently, `na.act` is not passed on to
#' the predict function, so the output needs to be manually padded with any NAs
#' present in the dataset. This is an inelegant way round the problem, but
#' leaves the original predict function as is.
#'
#' @param object :svyglm model object: the model to predict from.
#' @param newdata :data.frame or data.table: the data used for predictions.
#' @param type :character: predictions on link or response scale?
#' @param se.fit :boolean: return standard errors too? If so, these will be
#'        returned in a list, like predict.(g)lm.
#'
#' @return either a vector of length equal to the number of rows of newdata,
#' padded with NAs if necessary, with the predicted values from the model
#' object, or a list with the fit and standard errors, each vectors of length
#' equal to the number of rows of newdata.

PredictPaddedNAs = function(object, newdata = NULL, type = 'link', se.fit = FALSE)
{
  if (class(object)[1] != 'svyglm')
  {
    stop('Error in PredictPaddedNAs : object must be a svyglm model.')
  }

  cols    = all.vars(formula(object)[-2])
  w_comp  = complete.cases(subset(newdata, , cols))
  out_vec = rep(as.numeric(NA), length = nrow(newdata))

  if (isFALSE(se.fit))
  {
    pred = predict(object, newdata = newdata, type = type)
    out_vec[w_comp] = pred
    return(out_vec)
  } else
  {
    pred = predict(object, newdata = newdata, se.fit = se.fit, type = type)
    out  = list(fit = out_vec, se.fit = out_vec)
    out[['fit']][w_comp]    = as.numeric(pred)
    out[['se.fit']][w_comp] = sqrt(attr(pred, 'var'))
    return(out)
  }
}



#' Get predicted seroprevalences and incidences (the latter assuming specific
#' times to seroreversion for the three assays.
#'
#' @param d :data.table: data used in model fitting.
#' @param ex :list: list of length 4 numeric vectors, with the elements
#'        being the times to seroreversion for Abbott, Ortho, and Roche assays,
#'        and the lead/lag.
#' @param f :character: formula object to use (defined in functions.R).
#' @param get_cis :boolean: get uncertainty ranges (across different times to
#'        seroreversion and lead/lag times). This is the output from having run
#'        `seroprevalence_get_overall_uncertainty.R`.
#' 
#' @return Augmented data.table including predicted seroprevalences and
#' incidences, with and without various metrics of vaccination, for reference
#' and "best" (adjusted cases) models.

GetPredictions = function(d, ex, f, get_cis = TRUE)
{
  # Get the adjusted cases for the times to seroreversion in `ex`.
  s_ex = GetAdjustedCases(dat = d, durs = ex[[1]])
  s_ex = merge(d, s_ex[, c('state', 'week', 'cases_perc_waning')], all.x = TRUE)

  s_ex_dat = subset(s_ex, !is.na(seroprevalence))
  s_ex_dat = s_ex_dat[, wghts := state_population / n_total]

  # Fit GLM for both reference and waning model.
  sdesign = svydesign(id = ~ 1,
                      weights = ~ wghts,
                      data = s_ex_dat)
  m_ex = svyglm(get(f), 
                design = sdesign,
                na.action = na.exclude,
                family = 'binomial')
  m_ref = svyglm(get(gsub('_3w', '', f)), 
                 design = sdesign,
                 na.action = na.exclude,
                 family = 'binomial')

  if (isTRUE(get_cis))
  {
    cis = readRDS('../output/cis_glm_3-wanings.rds')

    s_ex = merge(s_ex, cis, by = c('week', 'state'), all.x = TRUE)
  }

  # Get predicted seroprevalences from the model, to include their 95%
  # confidence envelopes.
  link = family(m_ex)$linkinv

  p_seroprev = PredictPaddedNAs(m_ex, se.fit = TRUE, newdata = s_ex)

  # Get predicted seroprevalences, and the upper and lower confidence intervals
  # (from the specific model).
  s_ex = s_ex[, `:=`(pred_seroprev = link(p_seroprev$fit),
                     pred_seroprev_modlo = link(p_seroprev$fit - 
                                                (2 * p_seroprev$se.fit)),
                     pred_seroprev_modhi = link(p_seroprev$fit + 
                                                (2 * p_seroprev$se.fit)))]

  # Get predicted seroprevalences, in wich each of the assays are exclusively
  # used, in turn (i.e., predicted seroprevalence if Abbott assays are used
  # exclusively throughout the US).
  s_ex = s_ex[, `:=`(Abbott_Architect_p_real = Abbott_Architect_perc,
                     Ortho_VITROS_p_real = Ortho_VITROS_perc,
                     Roche_Elecsys_p_real = Roche_Elecsys_perc)]

  s_ex = s_ex[, `:=`(Abbott_Architect_perc = 100,
                     Ortho_VITROS_perc = 0,
                     Roche_Elecsys_perc = 0)]

  p_seroprev = PredictPaddedNAs(m_ex, se.fit = FALSE, newdata = s_ex) 
  s_ex = s_ex[, `:=`(pred_seroprev_abbott = link(p_seroprev))]

  s_ex = s_ex[, `:=`(Abbott_Architect_perc = 0,
                     Ortho_VITROS_perc = 100,
                     Roche_Elecsys_perc = 0)]

  p_seroprev = PredictPaddedNAs(m_ex, se.fit = FALSE, newdata = s_ex) 
  s_ex = s_ex[, `:=`(pred_seroprev_ortho = link(p_seroprev))]

  s_ex = s_ex[, `:=`(Abbott_Architect_perc = 0,
                     Ortho_VITROS_perc = 0,
                     Roche_Elecsys_perc = 100)]

  p_seroprev = PredictPaddedNAs(m_ex, se.fit = FALSE, newdata = s_ex) 
  s_ex = s_ex[, `:=`(pred_seroprev_roche = link(p_seroprev))]


  # Get predicted incidences, by using the model fit above to predict, but using
  # the real cumulative number of reported cases in place of the waned cases.
  s_ex = s_ex[, cases_perc_waning := cases_cumulative_perc]

  p_incidence = PredictPaddedNAs(m_ex, se.fit = TRUE, newdata = s_ex)

  # Get estimated proportions infected, and the upper and lower confidence
  # intervals (from the specific model).
  s_ex = s_ex[, `:=`(pred_incidence = link(p_incidence$fit),
                     pred_incidence_modlo = link(p_incidence$fit - 
                                                 (2 * p_incidence$se.fit)),
                     pred_incidence_modhi = link(p_incidence$fit + 
                                                 (2 * p_incidence$se.fit)))]


  #' Estimated EPIV, assuming no correlation between vaccination and incidence.
  GetNoCor = function(x, y)
  {
    return(x + y - x * y)
  }

  #' Estimated EPIV, assuming a negative correlation between vaccination and
  #' incidence.
  GetNegCor = function(x, y)
  {
    return(pmin(x + y, 1))
  }

  s_ex = s_ex[, `:=`(prevalence = NA,
                     prevalence_no_vacc = NA,
                     seroprev_v = GetNoCor(seroprevalence / 100, 
                                           vaccinated_cumulative_perc / 100), 
                     seroprev_v_neg = GetNegCor(seroprevalence / 100, 
                                                vaccinated_cumulative_perc / 100), 
                     pred_epiv = GetNoCor(pred_incidence, 
                                          vaccinated_cumulative_perc / 100),
                     pred_epiv_neg = GetNegCor(pred_incidence,
                                               vaccinated_cumulative_perc / 100),
                     pred_epiv_d1 = GetNoCor(pred_incidence,
                                                    vaccinated_dose1_cumulative_perc / 100),
                     pred_epiv_d1_neg = GetNegCor(pred_incidence,
                                                  vaccinated_dose1_cumulative_perc / 100))]

  if (isTRUE(get_cis))
  {
    s_ex = s_ex[, `:=`(pred_epiv_lo = GetNoCor(pred_incidence_lo,
                                               vaccinated_cumulative_perc / 100),
                       pred_epiv_hi = GetNoCor(pred_incidence_hi, 
                                               vaccinated_cumulative_perc / 100),
                       pred_epiv_neg = GetNegCor(pred_incidence,
                                                 vaccinated_cumulative_perc / 100),
                       pred_epiv_d1_lo = GetNoCor(pred_incidence_lo,
                                                  vaccinated_dose1_cumulative_perc / 100), 
                       pred_epiv_d1_hi = GetNoCor(pred_incidence_hi,
                                                  vaccinated_dose1_cumulative_perc / 100),
                       pred_epiv_d1_neg_lo = GetNegCor(pred_incidence_lo, 
                                                       vaccinated_dose1_cumulative_perc / 100), 
                       pred_epiv_d1_neg_hi = GetNegCor(pred_incidence_hi, 
                                                       vaccinated_dose1_cumulative_perc / 100))]
  }

  return(s_ex)
}




#' Get predicted seroprevalences when a single assay is used.
#'
#' @param d :data.table: the data.
#' @param ex :numeric vector: times to seroreversion for each of the assays and
#'        the lead/lag between a case being reported and seroconverting. The
#'        vector is assumed to be named with "abbott", "ortho", "roche", "lag".
#'
#' @return An augmented data.table with predicted seroprevalences for each of
#' the assays.

GetPredsPerAssay = function(d, ex, f)
{
  d = copy(d)
  s_ex = GetAdjustedCases(dat = d, durs = ex)
  s_ex = s_ex[, wghts := state_population / n_total]

  sdesign = svydesign(id = ~ 1,
                      weights = ~ wghts,
                      data = s_ex)
  m_ex = svyglm(get(f), 
                design = sdesign,
                na.action = na.exclude,
                family = 'binomial')

  link = family(m_ex)$linkinv

  s_ex = s_ex[, `:=`(Abbott_Architect_perc = 0,
                     Ortho_VITROS_perc = 0,
                     Roche_Elecsys_perc = 100)]

  p_seroprev = PredictPaddedNAs(m_ex, se.fit = FALSE, newdata = s_ex) 
  s_ex = s_ex[, `:=`(pred_seroprev_roche = link(p_seroprev))]

  s_ex = s_ex[, `:=`(Abbott_Architect_perc = 0,
                     Ortho_VITROS_perc = 100,
                     Roche_Elecsys_perc = 0)]

  p_seroprev = PredictPaddedNAs(m_ex, se.fit = FALSE, newdata = s_ex)
  s_ex = s_ex[, `:=`(pred_seroprev_ortho = link(p_seroprev))]

  s_ex = s_ex[, `:=`(Abbott_Architect_perc = 100,
                     Ortho_VITROS_perc = 0,
                     Roche_Elecsys_perc = 0)]

  p_seroprev = PredictPaddedNAs(m_ex, se.fit = FALSE, newdata = s_ex)
  s_ex = s_ex[, `:=`(pred_seroprev_abbott = link(p_seroprev))]

  # Get difference between the estimated seroprevalences, and predicted
  # seroprevalences when a single assay is used.
  s_ex$pred_seroprev_ortho_d  = s_ex$pred_seroprev_ortho - s_ex$seroprevalence / 100
  s_ex$pred_seroprev_abbott_d = s_ex$pred_seroprev_abbott - s_ex$seroprevalence / 100
  s_ex$pred_seroprev_roche_d  = s_ex$pred_seroprev_roche - s_ex$seroprevalence / 100

  s_ex = pivot_longer(s_ex, c('pred_seroprev_ortho_d',
                              'pred_seroprev_abbott_d',
                              'pred_seroprev_roche_d'), 
                      values_to = 'preds_diff', names_to = 'assay_preds')
  s_ex$round_pl = paste('Round', s_ex$survey_round)
  s_ex$assay_preds[s_ex$assay_preds == 'pred_seroprev_abbott_d'] = 'Abbott'
  s_ex$assay_preds[s_ex$assay_preds == 'pred_seroprev_ortho_d']  = 'Ortho'
  s_ex$assay_preds[s_ex$assay_preds == 'pred_seroprev_roche_d']  = 'Roche'

  return(as.data.table(s_ex))
}





#' Get predicted US-wide time series of seroprevalence per assay.
#'
#' @param dat :data.table: main data (with seroprevalence etc. per state.
#'
#' @return A data.table with overall seroprevalence across the US, and with
#' prediced seroprevalences for each of the assays (if each assay in turn had
#' been exclusively used by all states).

GetPredTS = function(dat)
{
  s_all = dat[, .(seroprev_us = weighted.mean(seroprevalence / 100, 
                                              state_population), 
                  pred_seroprev_ortho = weighted.mean(pred_seroprev_ortho, 
                                                      state_population),
                  pred_seroprev_abbott = weighted.mean(pred_seroprev_abbott, 
                                                       state_population),
                  pred_seroprev_roche = weighted.mean(pred_seroprev_roche, 
                                                      state_population)),
              by = 'survey_round']
  rs = dat[, .(week = median(week)), by = 'survey_round']
  s_all = merge(s_all, rs, by = 'survey_round')

  s_all = s_all[order(week), ]

  s_all = pivot_longer(data = s_all, 
                       cols = c('seroprev_us', 'pred_seroprev_abbott',
                                'pred_seroprev_ortho', 'pred_seroprev_roche'),
                       values_to = 'values', 
                       names_to = 'assay')
  s_all$assay_lab = NA
  s_all$assay_lab[s_all$assay == 'seroprev_us']          = 'Surveys'
  s_all$assay_lab[s_all$assay == 'pred_seroprev_abbott'] = 'Abbott'
  s_all$assay_lab[s_all$assay == 'pred_seroprev_ortho']  = 'Ortho'
  s_all$assay_lab[s_all$assay == 'pred_seroprev_roche']  = 'Roche'

  s_all$assay_lab = factor(s_all$assay_lab, levels = c('Surveys', 'Abbott',
                                                       'Ortho', 'Roche'))

  return(as.data.table(s_all))
}





###############################################################################
# Model formulae.
###############################################################################

#' Replace a variable in a given formula.
#'
#' @param f :model formula:
#' @param variable :character: the variable to be replaced.
#' @param replacement :character: the replacement.
#'
#' @return a formula with one of its variables replaced.

ReplaceVarFormula = function(f, variable, replacement)
{
  f = gsub(variable, replacement, f)
  f = as.formula(paste(f[2], f[1], f[3]))
  return(f)
}



# Main reference GLM model.
f_glm = as.formula(cbind(n_positive, n_negative) ~ 
                   week_decimal * factor(state) +
                   sqrt(cases_cumulative_perc) +
                   sqrt(deaths_cumulative_perc) +
                   excess_deaths_difference_perc +
                   sqrt(hospitalisation_cumulative_perc) +
                   log(tests_total_cumulative_perc) +
                   vaccinated_cumulative_perc +
                   Abbott_Architect_perc +
                   Roche_Elecsys_perc +
                   cases_00_19_perc +
                   cases_20_49_perc +
                   cases_50_69_perc)

# Waning models (only difference is that the cumulative reported cases from the
# reference model is replaced with cumulative number of "waned" cases).
f_glm_3w = ReplaceVarFormula(f = f_glm, 
                             variable = 'cases_cumulative_perc', 
                             replacement = 'cases_perc_waning')

# Main reference GLM model with splines.
f_spline = as.formula(cbind(n_positive, n_negative) ~ 
                      week_decimal * state +
                      ns(sqrt(cases_cumulative_perc), df = 5) +
                      ns(sqrt(deaths_cumulative_perc), df = 5) + 
                      ns(excess_deaths_difference_perc, df = 5) +
                      ns(sqrt(hospitalisation_cumulative_perc), df = 5) +
                      ns(log(tests_total_cumulative_perc), df = 5) +
                      vaccinated_cumulative_perc +
                      Abbott_Architect_perc +
                      Roche_Elecsys_perc +
                      ns(cases_00_19_perc, df = 5) +
                      ns(cases_20_49_perc, df = 5) +
                      ns(cases_50_69_perc, df = 5))

