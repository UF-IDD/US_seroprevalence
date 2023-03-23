###############################################################################
# Get overall uncertainty around fits and estimates.
# ~~~
#
# The uncertainty estimated here is across the waning GLM models, to account for
# uncertainty in the selection of the times to seroreversion of the three
# assays, and the lead/lag time. These estimates should include uncertainty in
# each model fit. US-wide estimates only include uncertainty around times to
# seroreversion and lead/lag time, but not uncertainty around each individual
# model fit (which in any case is far smaller).
#
# This code was intended to run on an HPC; it can take some time if run
# serially or in parallel on a personal computer.
###############################################################################

library('survey')
library('data.table')


source('functions.R')

# Get array number from the HPC.
hpc = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load complete data. This includes time points for which seroprevalence
# estimates are not available (this will be necessary to produce time series of
# cases adjusted for waning).
states_w = fread('../data/supporting_data.csv')
states_w = states_w[, week := as.Date(week)]

# Models with cases adjusted for 3 different waning rates. 
# This assumes that `seroprevalence_glm_3wanings.R` has been run, and that its
# outputs have been combined into one big data.table.
out_3w_glm = readRDS('../output/glm_3-wanings.rds')

# Select the best 5 percentile models by LOO median RMSE.
w_bests = which(out_3w_glm$rmse_loov <= quantile(out_3w_glm$rmse_loov, prob = 0.05))
bests = out_3w_glm[w_bests, c('waning_wks_abbott', 'waning_wks_ortho',
                              'waning_wks_roche', 'detects')]
bests = as.matrix(bests)
bests = split(bests, seq_len(nrow(bests)))
bests = lapply(X = bests, FUN = function(x) 
               {
                 names(x) = c('abbott', 'ortho', 'roche', 'lag')
                 return(x)
               })

n = length(bests)
n = n / 324

# Run only a subset on this CPU.
r = (n * (hpc - 1) + 1):(hpc * n)
bests = bests[r]

# Get estimated proportion infected (and associated uncertainty envelopes) for
# each of the models in the best 5 percentile by LOO median RMSE (so for all
# the different assumed times to seroreversion and lead/lag times included in
# that best 5 percentile). This yields a number of time series (with associated
# envelopes due to uncertainty in the predicted values) for each of those
# models. From these, for each point in time, get the range across the
# predictions (including the uncertainty envelopes). Do this for each state, and
# also get a US-wide estimate (which excludes each model-specific uncertainty
# envelopes).
preds_out = GetPredictions(d = states_w, ex = bests[1], f = 'f_glm_3w', 
                           get_cis = FALSE)

preds_out_us = preds_out[, .(state = 'US-wide', 
                             week = mean(week),
                             state_population = as.character(NA),
                             pred_seroprev_modlo = weighted.mean(pred_seroprev, 
                                                                 w = state_population),
                             pred_seroprev_modhi = weighted.mean(pred_seroprev, 
                                                                 w = state_population),
                             pred_incidence_modlo = weighted.mean(pred_incidence, 
                                                                  w = state_population),
                             pred_incidence_modhi = weighted.mean(pred_incidence, 
                                                                  w = state_population)),
                         by = 'survey_round']

preds_out = rbind(preds_out[, c('state', 'week', 'state_population',
                                'pred_seroprev_modlo', 'pred_seroprev_modhi',
                                'pred_incidence_modlo', 'pred_incidence_modhi')],
                  preds_out_us[, survey_round := NULL])

for (i in 2:length(bests))
{
  preds = GetPredictions(d = states_w, ex = bests[i], f = 'f_glm_3w', 
                         get_cis = FALSE)

  preds_us = preds[, .(state = 'US-wide', 
                       week = mean(week),
                       state_population = NA,
                       pred_seroprev_modlo = weighted.mean(pred_seroprev, 
                                                           w = state_population),
                       pred_seroprev_modhi = weighted.mean(pred_seroprev, 
                                                           w = state_population),
                       pred_incidence_modlo = weighted.mean(pred_incidence, 
                                                            w = state_population),
                       pred_incidence_modhi = weighted.mean(pred_incidence, 
                                                            w = state_population)),
                   by = 'survey_round']

  preds = rbind(preds[, c('state', 'week', 'state_population',
                          'pred_seroprev_modlo', 'pred_seroprev_modhi',
                          'pred_incidence_modlo', 'pred_incidence_modhi')],
                preds_us[, survey_round := NULL])

  preds_out = rbind(preds_out, preds)
}

out = preds_out[, .(pred_seroprev_lo = min(pred_seroprev_modlo),
                    pred_seroprev_hi = max(pred_seroprev_modhi),
                    pred_incidence_lo = min(pred_incidence_modlo),
                    pred_incidence_hi = max(pred_incidence_modhi)), 
                by = c('state', 'week')]

saveRDS(out, file = paste0('../output/cis/cis_glm_3-wanings_', hpc, '.rds'))

