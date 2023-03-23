###############################################################################
# GLM for the main reference model.
# ~~~
#
# This script fits the main reference model, which assumes no waning in cases.
###############################################################################

library('survey')
library('dplyr')
library('data.table')

source('functions.R')

# Load complete data. This includes time points for which seroprevalence
# estimates are not available (this will be necessary to produce time series of
# cases adjusted for waning).
dat = fread('../data/supporting_data.csv')
dat = dat[, week := as.Date(week)]

# Subset of the data with seroprevalence estimates.
dat_rounds = subset(dat, !is.na(seroprevalence))

# Weights are the inverse of the sampling proportion.
dat_rounds = dat_rounds[, wghts := state_population / n_total]

# Main model.
sdesign = svydesign(id = ~ 1,
                    weights = ~ wghts,
                    data = dat_rounds)
m = svyglm(f_glm, 
           design = sdesign,
           na.action = na.exclude,
           family = 'binomial')


# Get LOO RMSE values - leaving one round out at a time.
rs    = sort(unique(dat_rounds$survey_round))
rmses = rep(NA, length(rs))

for (j in seq_along(rs))
{
  # Split out current round j from the data.
  this_s = dat_rounds[abs(dat_rounds$survey_round - rs[j]) > 1e-2, ]
  this_r = dat_rounds[abs(dat_rounds$survey_round - rs[j]) < 1e-2, ]

  # Fit model to data excluding current round j.
  this_sdesign = svydesign(id = ~ 1,
                           weights = ~ wghts,
                           data = this_s)
  this_m = svyglm(f_glm, 
                  design = this_sdesign,
                  na.action = na.exclude,
                  family = 'binomial')

  if (abs(length(this_m$residuals) - nrow(this_s)) > 1e-2)
    stop('Error in seroprevalence_glm_3wanings.R : model dim not the same as data')

  # Predict round j using model.
  p_r      = predict(this_m, newdata = this_r, type = 'response')
  rmses[j] = RMSE(dat = this_r$seroprevalence / 100, pred = p_r)
}

# Model RMSE and LOO RMSE
m$rmse      = RMSE(dat = dat_rounds$seroprevalence / 100, 
                   pred = predict(m, type = 'response'))
m$rmse_loov = median(rmses, na.rm = TRUE)

saveRDS(m, file = '../output/main_glm.rds')

