###############################################################################
# Run the main GAM reference model (and its leave-one-out RMSE).
# ~~~
#
# This script fits the main reference model, which assumes no waning in cases.
#
# This code will run in parallel on a personal computer.
###############################################################################

library('mgcv')
library('data.table')
library('doParallel')

source('functions.R')


# Load complete data. This includes time points for which seroprevalence
# estimates are not available (this will be necessary to produce time series of
# cases adjusted for waning).
dat = fread('../data/supporting_data.csv')
dat = dat[, week := as.Date(week)]

# Subset of the data with seroprevalence estimates.
dat_rounds = subset(dat, !is.na(seroprevalence))

# Main model.
# Weights account for difference sampling proportions across states. Note that
# here, weights are normalised by their mean.
m = gam(f_gam, data = dat_rounds, 
        weights = (state_population / n_total) / mean(state_population / n_total),
        method = 'REML', select = TRUE, family = 'binomial')

m$rmse = RMSE(dat = dat_rounds$seroprevalence / 100, 
              pred = predict(m, type = 'response'))


rs = sort(unique(dat$survey_round))

cl = makeCluster(10)
registerDoParallel(cl)  # Register the cluster.

mods = foreach(i = seq_along(rs), .packages = c('mgcv', 'data.table')) %dopar%
{
  this_s = dat_rounds[abs(survey_round - rs[i]) > 1e-2, ]
  this_r = dat_rounds[abs(survey_round - rs[i]) < 1e-2, ]

  this_m = gam(f_gam, data = this_s, 
               weights = (state_population / n_total) / mean(state_population / n_total),
               method = 'REML', select = TRUE, 
               family = 'binomial')
  p_r = predict.gam(this_m, newdata = this_r, type = 'response')
  return(RMSE(dat = this_r$seroprevalence / 100, pred = p_r))
}

stopCluster(cl)

m$rmse_loov = median(unlist(mods))

saveRDS(m, file = '../output/main_gam.rds')

