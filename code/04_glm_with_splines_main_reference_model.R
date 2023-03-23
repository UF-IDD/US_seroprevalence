###############################################################################
# Run the main reference model with splines.
# ~~~
#
# This script fits the main reference model, which assumes no waning in cases,
# and allows for non-parametric nonlinear functions between variables and
# seroprevalence.
###############################################################################

library('data.table')
library('survey')
library('splines')

source('functions.R')


# Load complete data. This includes time points for which seroprevalence
# estimates are not available (this will be necessary to produce time series of
# cases adjusted for waning).
dat = fread('../data/supporting_data.csv')
dat = dat[, week := as.Date(week)]

# Subset of the data with seroprevalence estimates.
dat_rounds = subset(dat, !is.na(seroprevalence))
dat_rounds = dat_rounds[, wghts := state_population / n_total]

# Main model.
# Weights account for different sampling proportions across states. 
sdesign = svydesign(id = ~ 1,
                    weights = ~ wghts,
                    data = dat_rounds)

m = svyglm(f_spline, 
           design = sdesign,
           family = 'binomial')

m$rmse = RMSE(dat = dat_rounds$seroprevalence / 100, 
              pred = predict(m, type = 'response'))

saveRDS(m, file = '../output/main_glm_splines.rds')
