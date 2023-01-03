###############################################################################
# GLMs for waning models.
# ~~~
#
# This script fits the waning models, which makes many different assumptions on
# waning.
#
# This code was intended to run on an HPC; it can take some time if run
# seriously or in parallel on a personal computer.
###############################################################################

library('dplyr')
library('data.table')

source('functions.R')

# Get array number from the HPC.
hpc = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load complete data. This includes time points for which seroprevalence
# estimates are not available (this will be necessary to produce time series of
# cases adjusted for waning).
dat = fread('supporting_data.csv')
dat = dat[, week := as.Date(week)]

# Subset of the data with seroprevalence estimates.
dat_rounds = subset(dat, !is.na(seroprevalence))



###############################################################################
# Use cumulative cases adjusted for waning instead. 
############################################################################### 
detects = -1:2  # Lead/lag between a case being reported and seroconversion.

max_date = max(dat_rounds$week)

# Total number of weeks in the study.
m_wks = as.numeric(floor((max_date - as.Date('2020-03-01')) / 7)) 

# Number of weeks to seroreversion to explore in the analysis.
waning_wks = seq(8, m_wks, 1)

# Maximum number of weeks we can explore for the Abbott and Ortho assays
# (because both were phased out at different times). The Roche assay was used
# throughout the study, so its maximum is m_wks.
abbott_max_weeks = 70
ortho_max_weeks  = 49

# All combinations of lead/lags and times to seroreversion across the three
# assays that will be explored.
pars = expand.grid(detects = detects, 
                   waning_wks_abbott = waning_wks[waning_wks <= abbott_max_weeks], 
                   waning_wks_ortho = waning_wks[waning_wks <= ortho_max_weeks],
                   waning_wks_roche = waning_wks,
                   aic = NA,
                   rmse = NA,
                   rmse_loov = NA)

# We run all combinations across 2940 different processes on the HPC. 
# Get the combinations that will be run in this process.
n = nrow(pars)
n = n / 2940

r = (n * (hpc - 1) + 1):(hpc * n)
pars = pars[r, ]

for (i in seq_len(nrow(pars))) 
{
  # Get waned cases for current assumptions.
  s = GetAdjustedCases(dat = dat, 
                       durs = c('abbott' = pars$waning_wks_abbott[i],
                                'ortho' = pars$waning_wks_ortho[i], 
                                'roche' = pars$waning_wks_roche[i], 
                                'lag' = pars$detects[i]))
  s = s[, wghts := mean(n_total / state_population) * state_population / n_total]

  m = glm(f_glm_3w, data = s, 
          weights = wghts,
          family = 'binomial')

  pars$aic[i]  = AIC(m)
  pars$rmse[i] = RMSE(dat = s$seroprevalence / 100, pred = predict(m, type = 'response'))

  rs    = sort(unique(s$survey_round))
  rmses = rep(NA, length(rs))

  for (j in seq_along(rs))
  {
    this_s = s[abs(s$survey_round - rs[j]) > 1e-2, ]
    this_r = s[abs(s$survey_round - rs[j]) < 1e-2, ]

    this_m = glm(f_glm_3w, data = this_s, 
                 weights = wghts,
                 family = 'binomial')

    p_r      = predict(this_m, newdata = this_r, type = 'response')
    rmses[j] = RMSE(dat = this_r$seroprevalence / 100, pred = p_r)
  }

  pars$rmse_loov[i] = median(rmses, na.rm = TRUE)
}

saveRDS(pars, file = paste0('../output/glms/glm_3-wanings_', hpc, '.rds'))

