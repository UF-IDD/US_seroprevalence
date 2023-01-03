###############################################################################
# Produce figures for manuscript and supporting information.
###############################################################################

library('data.table')
library('mgcv')
library('colorspace')
library('dplyr')
library('tidyr')
library('fields')
library('RColorBrewer')
library('tidycensus')
library('doParallel')
library('ggplot2')
library('ggpubr')
library('corrplot')
library('Hmisc')

source('functions.R')




###############################################################################
# Script.
###############################################################################

# Data.
states_w = fread('../data/supporting_data.csv')

# Keep a subset which is only rows with seroprevalence.
s = states_w[!is.na(seroprevalence), ]


# Reference model (no adjustment for waning).
m_ref_glm = readRDS('../output/main_glm.rds')
m_ref_gam = readRDS('../output/main_gam.rds')

# Models with cases adjusted for 3 different waning rates. 
# This assumes that `seroprevalence_glm_3wanings.R` has been run, and that its
# outputs have been combined into one big data.table.
out_3w_glm = readRDS('../output/glm_3-wanings.rds')

# Best model by metric.
out_3w_glm[which.min(out_3w_glm$aic), ]
out_3w_glm[which.min(out_3w_glm$rmse), ]
(b = out_3w_glm[which.min(out_3w_glm$rmse_loov), ])

# We'll be using the best model below.
best = list(c('abbott' = b$waning_wks_abbott, 
              'ortho' = b$waning_wks_ortho, 
              'roche' = b$waning_wks_roche, 
              'lag' = b$detects))



###############################################################################
# Fig 1: Data
###############################################################################

s_sub = subset(s, survey_round == 24)
ord   = order(s_sub$cases_cumulative_perc, decreasing = TRUE)
sero  = s_sub$seroprevalence[ord] * 100
case  = s_sub$cases_cumulative_perc[ord]
death = s_sub$deaths_cumulative_perc[ord]
vacc  = s_sub$vaccinated_cumulative_perc[ord]

w_abb = which(s_sub$Abbott_Architect_p[ord] > 99)

pdf('../plots/figure1.pdf', width = 8.7/2.54, height = 11/2.54, pointsize = 6)
par(mar = c(2.5, 1.8, 1.5, 0.8), mgp = c(1.3, 0.3, 0), tcl = 0.20, las = 1)

split.screen(c(1, 4))
par(oma = c(0, 0, 0, 0.9))

screen(1)
plot(case, seq_along(ord), type = 'o', pch = 16, yaxt = 'n', ylim = c(1.7, 48.3),
     xlab = '% cumul. cases', ylab = '')
abline(h = w_abb, col = adjustcolor('black', alpha = 0.07), lwd = 7)
axis(side = 2, at = seq_along(ord), labels = s_sub$state[ord])
axis(side = 4, at = seq_along(ord), labels = FALSE)
axis(side = 3, labels = FALSE)
mtext('A', side = 3, adj = -0.02, line = 0.1, cex = 1.4)

screen(2)
plot(sero, seq_along(ord), type = 'o', pch = 16, yaxt = 'n', ylim = c(1.7, 48.3),
     xlab = 'Seroprevalence', ylab = '')
abline(h = w_abb, col = adjustcolor('black', alpha = 0.07), lwd = 7)
lines(sero, seq_along(ord), type = 'o', pch = 16)
mgp.axis(side = 2, at = seq_along(ord), labels = s_sub$state[ord], hadj = 0.5,
         mgp = c(1.3, 1.3, 0))
axis(side = 4, at = seq_along(ord), labels = FALSE)
axis(side = 3, labels = FALSE)
mtext('B', side = 3, adj = -0.02, line = 0.1, cex = 1.4)

screen(3)
plot(death,seq_along(ord),  type = 'o', pch = 16, yaxt = 'n', ylim = c(1.7, 48.3),
     xlab = '% cumul. deaths', ylab = '')
mgp.axis(side = 2, at = seq_along(ord), labels = s_sub$state[ord], hadj = 0.5,
         mgp = c(1.3, 1.3, 0))
axis(side = 4, at = seq_along(ord), labels = FALSE)
axis(side = 3, labels = FALSE)
mtext('C', side = 3, adj = -0.02, line = 0.1, cex = 1.4)

screen(4)
plot(vacc, seq_along(ord), type = 'o', pch = 16, yaxt = 'n', ylim = c(1.7, 48.3),
     xlab = '% vaccinated', ylab = '')
mgp.axis(side = 2, at = seq_along(ord), labels = s_sub$state[ord], hadj = 0.5,
         mgp = c(1.3, 1.3, 0))
axis(side = 4, at = seq_along(ord), labels = s_sub$state[ord])
axis(side = 3, labels = FALSE)
mtext('D', side = 3, adj = -0.02, line = 0.1, cex = 1.4)

close.screen(all.screens = TRUE)
dev.off()




################################################################################
# Fig 2: Model performance when cases are adjusted using three different waning
# rates (one per assay).
################################################################################

b_aic       = out_3w_glm[which.min(out_3w_glm$aic), ]
b_loov_rmse = out_3w_glm[which.min(out_3w_glm$rmse_loov), ]

diag_df = data.frame(x = c(0, 2, 2, 2, 48, 48, 54, 0, 2, 2, 2, 24, 24, 54), 
                     y = c(0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0),
                     w = factor(c(rep('Slower waning', 7), rep('Faster waning', 7)),
                                levels = c('Slower waning', 'Faster waning')))

cols = brewer.pal(n = 3, name = 'Dark2')

p1 = ggplot(diag_df, aes(x = x, y = y, colour = w)) + 
  geom_line(alpha = 0.5, size = 0.8) + 
  labs(x = 'Weeks', y = 'Prob. seropos.') +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.5, 1)) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = 'none')

# Compare three time series: cumulative cases (no waning), and cumulative cases,
# adjusted for two different times to seroreversion, faster and slower.

# No waning.
fl = subset(states_w, state == 'FL')$cases_cumulative
fl_0 = subset(states_w, state == 'FL')$cases

# Assuming 24 weeks of time to seroreversion.
fl_24 = GetWanedCases(x = fl_0, detect = 0, 
                      ab_dur = c(abbott = 24, ortho = 60, roche = 97),
                      props_abbott = rep(1, length(fl)),
                      props_ortho = rep(0, length(fl)),
                      props_roche = rep(0, length(fl)))

# Assuming 60 weeks of time to seroreversion.
fl_60 = GetWanedCases(x = fl_0, detect = 0, 
                      ab_dur = c(abbott = 24, ortho = 60, roche = 97),
                      props_abbott = rep(0, length(fl)),
                      props_ortho = rep(1, length(fl)),
                      props_roche = rep(0, length(fl)))

df = data.frame(x = rep(seq_along(fl), 3),
                fl = c(fl, diag(as.matrix(fl_60)), diag(as.matrix(fl_24))),
                w = factor(c(rep('No waning', length(fl)), 
                             rep('Slower waning', nrow(fl_24)), 
                             rep('Faster waning', nrow(fl_24))),
                           levels = c('No waning', 'Slower waning', 'Faster waning')))

p2 = ggplot(df, aes(x = x, y = fl / 1e6, colour = w)) + 
  geom_line(size = 0.8, alpha = 0.5) +
  labs(x = 'Week', y = 'Cases (millions)', colour = NULL) +
  scale_colour_manual(values = c('black', cols[1:2])) +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.spacing.x = unit(0, 'lines'),
        legend.key.height = unit(0.3, 'cm'),
        legend.justification = c(-0.4, 1),
        legend.margin = margin(-23, 0, 0, 0),
        legend.text = element_text(margin = margin(r = -1.2, unit = 'cm'), size = 7))

w_b_aic  = quantile(out_3w_glm$aic, prob = 0.05) - m_ref_glm$aic
w_b_loov = quantile(out_3w_glm$rmse_loov, prob = 0.05) - m_ref_glm$rmse_loov

p3 = ggplot(subset(out_3w_glm, detects < 0 &
                   waning_wks_ortho == 49), 
            aes(x = waning_wks_roche, y = waning_wks_abbott)) + 
  geom_tile(aes(fill = aic - m_ref_glm$aic)) +
  geom_line(data = data.frame(x = c(7.5, 70.5), y = c(7.5, 70.5)), 
            aes(x = x, y = y), colour = 'white', linetype = '33') +
  geom_contour(aes(z = aic - m_ref_glm$aic), breaks = w_b_aic, 
               col = 'black', alpha = 0.6) +
  geom_point(data = b_aic, size = 1, colour = 'green',
             aes(x = waning_wks_roche, y = waning_wks_abbott, fill = NULL)) + 
  geom_text(data = data.frame(x = 95, y = 65), aes(x = x, y = y), 
            label = expression(Delta*'AIC'), 
            hjust = 1, size = 3) +
  coord_fixed() + 
  labs(x = NULL,
       y = NULL,
       fill = NULL) +
  scale_fill_continuous_diverging(palette = 'Blue-Red 3',
                                  breaks = c(-1e3, 0, 1e3, 2e3)) +
  theme_bw() +
  theme(legend.position = 'top', 
        legend.text = element_text(size = 7),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-0, -10, -5, -10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))
    

p4 = ggplot(subset(out_3w_glm, detects < 0 &
                   waning_wks_ortho == 10), 
            aes(x = waning_wks_roche, y = waning_wks_abbott)) + 
  geom_tile(aes(fill = rmse_loov - m_ref_glm$rmse_loov)) +
  geom_line(data = data.frame(x = c(7.5, 70.5), y = c(7.5, 70.5)), 
            aes(x = x, y = y), colour = 'white', linetype = '33') +
  geom_contour(aes(z = rmse_loov - m_ref_glm$rmse_loov), breaks = w_b_loov, 
               col = 'black', alpha = 0.6) +
  geom_point(data = b_loov_rmse, size = 1, colour = 'green',
             aes(x = waning_wks_roche, y = waning_wks_abbott, fill = NULL)) + 
  geom_text(data = data.frame(x = 95, y = 65), aes(x = x, y = y), 
            label = expression(Delta*'LOO RMSE'), 
            hjust = 1, size = 3) +
  coord_fixed() + 
  labs(x = NULL, 
       y = NULL,
       fill = NULL) +
  scale_fill_continuous_diverging(palette = 'Blue-Red 3',
                                  breaks = c(-0.002, -0.001, 0, 0.001, 0.002, 0.003)) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.text = element_text(size = 7),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-0, -10, -5, -10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

p34 = ggarrange(p3, p4, nrow = 1)
p34 = annotate_figure(p34, 
                      bottom = text_grob('Roche time to seroreversion (weeks)',
                                         size = 9),
                      left = text_grob('Abbott time to seroreversion (weeks)',
                                       rot = 90, size = 9))

p = ggarrange(ggarrange(NULL, p1, NULL, p2, ncol = 1, 
                        heights = c(0.03, 0.32, -0.03, 0.5), align = 'v'), 
              NULL, p34, nrow = 1, widths = c(0.26, 0.03, 0.7))

ggsave('../plots/figure2.pdf', 
       width = 17.8, height = 18 / 3, unit = 'cm')





################################################################################
# Fig 3: Example time series of seroprevalence and predicted incidence.
################################################################################

preds_glm_all = GetPredictions(d = states_w, ex = best, f = 'f_glm_3w')

ex_states = c('ID', 'NE', 'NJ', 'NY', 'OK', 'LA', 'GA')

vacc_start = s[, .(start = max(week[vaccinated_cumulative_perc < 1e-2])), by = 'state']
vacc_start = subset(vacc_start, state %in% ex_states)
vacc_start$state = factor(vacc_start$state, levels = ex_states)


preds_glm_all_ex = subset(preds_glm_all, state %in% ex_states & 
                          week >= as.Date('2020-08-01') &
                          week <= as.Date('2022-01-15'))
preds_glm_all_ex = preds_glm_all_ex[, c('state', 'week', 'survey_round',
                                        'seroprevalence',
                                        'pred_seroprev', 'pred_incidence',
                                        'pred_seroprev_lo', 'pred_seroprev_hi',
                                        'pred_incidence_lo', 'pred_incidence_hi')]

# Uncertainty across models (to account for selection of times to seroreversion
# and lead/lag times).
#
# This assumes that `seroprevalence_get_overall_uncertainty.R` has been run, and
# that its outputs have been combined into one big data.table.
preds_glm_cis = readRDS('../output/cis_glm_3-wanings.rds')
preds_glm_cis_us = subset(preds_glm_cis, state == 'US-wide')

preds_glm_us = preds_glm_all[, .(state = 'US-wide', 
                                 week = mean(week),
                                 seroprevalence = weighted.mean(seroprevalence, 
                                                                w = state_population),
                                 pred_seroprev = weighted.mean(pred_seroprev, 
                                                               w = state_population),
                                 pred_incidence = weighted.mean(pred_incidence, 
                                                                w = state_population)),
                             by = 'survey_round']
preds_glm_us = merge(preds_glm_us, preds_glm_cis_us, all.x = TRUE)

preds_glm_all_ex = rbind(preds_glm_all_ex, preds_glm_us)
preds_glm_all_ex$state = factor(preds_glm_all_ex$state, 
                                levels = c(ex_states, 'US-wide'))




p = ggplot(preds_glm_all_ex) +
  geom_vline(data = vacc_start,
             aes(xintercept = start), linetype = '11', col = 'black', 
             alpha = 0.2, size = 1) +
  geom_ribbon(data = subset(preds_glm_all_ex, !is.na(pred_seroprev)),
              aes(x = week, ymin = pred_seroprev_lo, ymax = pred_seroprev_hi), 
              alpha = 0.1) +
  geom_ribbon(data = subset(preds_glm_all_ex, !is.na(pred_incidence_lo)),
              aes(x = week, ymin = pred_incidence_lo, ymax = pred_incidence_hi), 
              alpha = 0.2, fill = 'royalblue2') +
  geom_point(aes(x = week, y = seroprevalence / 100, colour = 'Seroprevalence surveys'),
             size = 0.8, alpha = 0.6) +
  geom_line(data = subset(preds_glm_all_ex, !is.na(pred_seroprev)),
            aes(x = week, y = pred_seroprev, colour = 'Fitted seroprevalence'), 
            size = 0.5, alpha = 1) + 
  geom_line(data = subset(preds_glm_all_ex, !is.na(pred_incidence)), 
            aes(x = week, y = pred_incidence, 
                colour = 'Estimated proportion infected'), 
            size = 0.5, alpha = 0.6) + 
  scale_x_date(date_breaks = '6 months', date_labels = '%b %Y') +
  facet_wrap(~ state, ncol = 4) + 
  labs(x = NULL, y = 'Seroprevalence and proportions infected') +
  scale_colour_manual('', 
                      breaks = c('Seroprevalence surveys',
                                 'Fitted seroprevalence',
                                 'Estimated proportion infected'),
                      values = c('Seroprevalence surveys' = 'black',
                                 'Fitted seroprevalence' = 'darkgray',
                                 'Estimated proportion infected' = 'royalblue2')) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, NA, NA),
                                                   linetype = c(NA, 1, 1),
                                                   alpha = c(0.6, 0.6, 0.6)))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = 'bottom',
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -15, -5, -10),
        strip.text = element_text(margin = margin(0.05, 0, 0.05, 0, 'cm')),
        panel.spacing = unit(0.1, 'cm'),
        plot.margin = margin(0, 0.34, 0.1, 0.1, 'cm'),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))

ggsave('../plots/figure3.pdf', width = 17.8, height = 7.0, unit = 'cm')




################################################################################
# Fig 4: Maps of datastreams and predicted incidence (with and without
# vaccination), negative correlations between incidence and vaccination, and
# example time series for two states.
################################################################################

#' Plot maps of various metrics, for a specific round.
#'
#' @param preds_all :data.table: augmented data.table with predicted values
#' (output from `GetPredictions` above).
#' @param r :numeric: survey round for which to plot the maps.
#' @param zl :numeric: the maximum value to use in the plots for seroprevalence
#' and incidence maps (all of which have common ranges in values).
#'
#' @return Four maps,
#' (i) anti-N seroprevalence (as per the surveys); 
#' (ii) estimated cumulative infections;
#' (iii) vaccination coverage;
#' (iv) estimated proportion infected and/or vaccinated.
#'
#' Plot of Pearson correlations between infections and vaccinations, per round.
#'
#' Example time series of anti-N seroprevalence, estimated cumulative
#' infections, and vaccination coverage, for two contrasting states.

PlotIncidenceVacc = function(preds_all, r = 26, zl = 0.7)
{
  preds = subset(preds_all, !is.na(seroprevalence))
  map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                  by = c('GEOID' = 'state_code'))
  map = left_join(map, subset(preds, survey_round %in% r), 
                  by = 'state')

  pars = sort(unique(map$par))

  pl1 = ggplot(map, aes(fill = seroprevalence / 100)) +
    geom_sf(colour = NA) +
    scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                     lim = c(0, zl)) +
    labs(fill = 'CDC\nseroprev.  ') +
    theme_void() +
    theme(legend.position = 'top', 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.width = unit(0.5, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          plot.title = element_text(hjust = 0.5))

  pl2 = ggplot(map, aes(fill = pred_incidence)) +
    geom_sf(colour = NA) +
    scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                     lim = c(0, zl)) +
    labs(fill = 'Estimated\nproportion  \ninfected') +
    theme_void() +
    theme(legend.position = 'top', 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.width = unit(0.5, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          plot.title = element_text(hjust = 0.5))

  pl3 = ggplot(map, aes(fill = vaccinated_cumulative_perc / 100)) +
    geom_sf(colour = NA) +
    scale_fill_continuous_sequential(palette = 'Blues', rev = TRUE,
                                     lim = c(0, zl)) +
    labs(fill = 'Vaccinated ') +
    theme_void() +
    theme(legend.position = 'top', 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.width = unit(0.5, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          plot.title = element_text(hjust = 0.5))

  pl4 = ggplot(map, aes(fill = pred_epiv)) +
    geom_sf(colour = NA) +
    scale_fill_continuous_sequential(palette = 'Blue-Yellow', rev = TRUE,
                                     lim = c(0, zl)) +
    labs(fill = 'EPIV ') +
    theme_void() +
    theme(legend.position = 'top', 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.width = unit(0.5, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          plot.title = element_text(hjust = 0.5))


  cors_i_v = preds[, .(week = mean(week),
                       r = cor.test(vaccinated_cumulative_perc, 
                                    pred_incidence)$estimate, 
                       r_raw = cor.test(vaccinated_cumulative_perc, 
                                        seroprevalence)$estimate, 
                       ci_lo = cor.test(vaccinated_cumulative_perc, 
                                        pred_incidence)$conf.int[1],
                       ci_hi = cor.test(vaccinated_cumulative_perc, 
                                        pred_incidence)$conf.int[2]), 
                   by = 'survey_round']
  cors_i_v = na.omit(cors_i_v)

  pl6 = ggplot(cors_i_v, aes(x = week, y = r)) + geom_point(size = 1) + geom_line() +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.1) +
    geom_hline(yintercept = 0, linetype = '22', size = 0.3) +
    scale_x_date(date_breaks = '3 months', date_labels = '%b %Y') +
    labs(x = NULL, y = 'Pearson corr.') +
    theme_bw() +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8))

  w = max(preds$week[preds$survey_round == r])

  pl7 = ggplot(subset(preds, state %in% c('WA', 'AL')), 
               aes(x = week, y = seroprevalence / 100)) +
    geom_ribbon(aes(ymin = pred_incidence_lo, ymax = pred_incidence_hi), 
                alpha = 0.2, fill = 'royalblue2') +
    geom_point(size = 0.5) +
    geom_vline(xintercept = w, linetype = '22', alpha = 0.3) +
    geom_line(aes(y = vaccinated_cumulative_perc / 100, 
                  colour = 'Proportion vaccinated')) +
    geom_line(data = subset(preds_all, state %in% c('WA', 'AL') & 
                            week >= as.Date('2020-08-01') &
                            week <= as.Date('2022-01-15')),
              aes(x = week, y = pred_incidence, 
                  colour = 'Estimated proportion infected')) +
    scale_x_date(date_breaks = '6 months', date_labels = '%b %Y') +
    facet_wrap(~ state, nrow = 2, ncol = 1) +
    theme_bw() +
    labs(x = NULL, y = 'Seroprev., vacc.') +
    scale_colour_manual('', 
                        breaks = c('Estimated proportion infected',
                                   'Proportion vaccinated'),
                        values = c('Estimated proportion infected' = 'royalblue2',
                                   'Proportion vaccinated' = 'red')) +
    theme(legend.position = 'none', 
          legend.spacing.y = unit(-0.2, 'cm'),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.key.width = unit(1, 'cm'),
          strip.background = element_blank(),
          strip.text = element_text(size = 7, margin = margin(0.15, 0, 0.05, 0, 'cm')),
          panel.spacing = unit(0, 'cm'),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          plot.margin = margin(0, 0.3, 0, 0, 'cm'))

  pl = ggarrange(ggarrange(pl1, pl2, nrow = 1, ncol = 2,
                           font.label = list(size = 8, face = 'plain'),
                           common.legend = FALSE, labels = c('A', 'B')), 
                 ggarrange(pl3, pl4, ncol = 2, 
                           font.label = list(size = 8, face = 'plain'),
                           common.legend = FALSE, labels = c('C', 'D')), 
                 ggarrange(pl6, pl7, ncol = 2,
                           font.label = list(size = 8, face = 'plain'),
                           common.legend = FALSE, widths = c(0.55, 0.45), 
                           labels = c('E', 'F')), 
                 nrow = 3, heights = c(0.35, 0.35, 0.3))
}

pl_glm = PlotIncidenceVacc(preds_all = preds_glm_all, r = 24, zl = 0.72)
ggsave('../plots/figure4.pdf', width = 11.4, height = 11.4, unit = 'cm')




###############################################################################
# Fig 5: Compare predicted incidences with blood donors dataset.
###############################################################################

don = fread('../data/blood_donors_surveys_processed.csv')
don$week = as.Date(don$week)


# Merge blood donors serosurvey with our own dataset. Need to do this manually
# to be able to have some flexibility in what dates are matched.
preds_glm_all = GetPredictions(d = states_w, ex = best, f = 'f_glm_3w')
preds_glm = subset(preds_glm_all, !is.na(seroprevalence))
preds_glm$blood_prevalence         = as.numeric(NA)
preds_glm$blood_prevalence_no_vacc = as.numeric(NA)

for (i in seq_len(nrow(preds_glm))) 
{
  w_st = which(don$state == preds_glm$state[i])
  w_wk = which(don$week[w_st] == preds_glm$week[i])

  if (length(w_wk) > 0) 
  {
    preds_glm$blood_prevalence[i] = don$blood_prevalence[w_st][w_wk]
    preds_glm$blood_prevalence_no_vacc[i] = don$blood_prevalence_no_vacc[w_st][w_wk]
  } else 
  {
    # If dates for rounds do not match, merge rounds that are within a week of
    # each other.
    d = abs(don$week[w_st] - as.Date(preds_glm$week[i]))

    if (any(d <= 7)) 
    {
      w_min = which(d == min(d))

      if (length(w_min) > 1) 
      {
        preds_glm$blood_prevalence[i] = mean(don$blood_prevalence[w_st][w_min])
        preds_glm$blood_prevalence_no_vacc[i] = mean(don$blood_prevalence_no_vacc[w_st][w_min])
      } else 
      {
        preds_glm$blood_prevalence[i] = don$blood_prevalence[w_st][w_min]
        preds_glm$blood_prevalence_no_vacc[i] = don$blood_prevalence_no_vacc[w_st][w_min]
      }
    }
  }
}

# Get LOESS predictions to add nonlinear trends to plots.
plot_preds = subset(preds_glm, !is.na(blood_prevalence) & !is.na(pred_epiv))

l = loess(pred_epiv ~ blood_prevalence, data = plot_preds)
plot_preds$pred_epiv_l = predict(l)

l = loess(pred_epiv_d1 ~ blood_prevalence, data = plot_preds)
plot_preds$pred_epiv_d1_l = predict(l)

l = loess(pred_epiv_d1_neg ~ blood_prevalence, data = plot_preds)
plot_preds$pred_epiv_d1_neg_l = predict(l)


b = brewer.pal(3, 'Set1')


# Top-left panel:
# - No correlation between vaccination and incidence.
# - Only include individuals with a complete series of vaccination.
plot_preds$lab = 'Ind. prob., complete series'

p1 = ggplot(plot_preds, aes(x = blood_prevalence, y = pred_epiv)) +
  geom_abline(intercept = 0, slope = 1, lty = '33') +
  geom_linerange(aes(ymin = pred_epiv_lo,
                     ymax = pred_epiv_hi), 
                 colour = b[2], size = 0.2, alpha = 0.3) +
  geom_point(pch = 16, size = 0.6, alpha = 0.3, colour = b[2]) +
  geom_line(aes(x = blood_prevalence, y = pred_epiv_l), 
            colour = b[2], size = 1) + 
  coord_fixed() + 
  lims(x = c(0, 1), y = c(0, 1)) + 
  facet_wrap(~ lab) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# Top-right panel:
# - No correlation between vaccination and incidence.
# - Include individuals with at least one dose of a vaccine.
plot_preds$lab = 'Independent prob., dose 1'

p2 = ggplot(plot_preds, aes(x = blood_prevalence, y = pred_epiv_d1)) +
  geom_abline(intercept = 0, slope = 1, lty = '33') +
  geom_linerange(aes(ymin = pred_epiv_d1_lo,
                     ymax = pred_epiv_d1_hi), 
                 colour = b[1], size = 0.2, alpha = 0.3) +
  geom_point(pch = 16, size = 0.8, alpha = 0.3, colour = b[1]) +
  geom_line(aes(x = blood_prevalence, y = pred_epiv_d1_l), 
            colour = b[1], size = 1) + 
  coord_fixed() + 
  lims(x = c(0, 1), y = c(0, 1)) + 
  facet_wrap(~ lab) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# Bottom-left panel:
# - Negative correlation between vaccination and incidence.
# - Only include individuals with a complete series of vaccination.
plot_preds$lab = 'Negative correlation, dose 1'

p3 = ggplot(plot_preds, aes(x = blood_prevalence, y = pred_epiv_d1_neg)) +
  geom_abline(intercept = 0, slope = 1, lty = '33') +
  geom_linerange(aes(ymin = pmin(pred_epiv_d1_neg_lo),
                     ymax = pmin(pred_epiv_d1_neg_hi)), 
                 colour = b[3], size = 0.2, alpha = 0.3) +
  geom_point(pch = 16, size = 0.8, alpha = 0.3, colour = b[3]) +
  geom_line(aes(x = blood_prevalence, y = pred_epiv_d1_neg_l), 
            colour = b[3], size = 1) + 
  coord_fixed() + 
  lims(x = c(0, 1), y = c(0, 1)) + 
  facet_wrap(~ lab) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# Bottom-right panel:
# Compare LOESS curves for the other three panels.
plot_preds$lab = 'Comparison'

p4 = ggplot(plot_preds, aes(x = blood_prevalence, y = pred_epiv_l)) +
  geom_abline(intercept = 0, slope = 1, lty = '33') +
  geom_line(colour = b[2], size = 1, show.legend = TRUE) + 
  geom_line(aes(x = blood_prevalence, y = pred_epiv_d1_l), 
            colour = b[1], size = 1) + 
  geom_line(aes(x = blood_prevalence, y = pred_epiv_d1_neg_l), 
            colour = b[3], size = 1) + 
  coord_fixed() + 
  lims(x = c(0, 1), y = c(0, 1)) + 
  facet_wrap(~ lab) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

pl = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, 
               font.label = list(face = 'plain'),
               common.legend = FALSE)
pl = annotate_figure(pl, 
                     bottom = text_grob('Anti-S seroprevalence from the blood donors dataset', 
                                        size = 9),
                     left = text_grob('Estimated proportion infected and/or vaccinated', 
                                      rot = 90, size = 9))

ggsave('../plots/figure5.pdf', width = 11.4, height = 11.4, units = 'cm')




################################################################################
# Supplementary information plots
################################################################################


###############################################################################
# Fig S1: Data
###############################################################################

s_ex = subset(s, survey_round %in% c(4, 13, 24))
s_ex = rbind(s_ex, s_ex[(nrow(s_ex) + 1):(nrow(s_ex) + 3), ])
s_ex$state[(nrow(s_ex) - 2):nrow(s_ex)] = 'ND'
s_ex$survey_round[(nrow(s_ex) - 2):nrow(s_ex)] = c(4, 13, 24)

s_ex = s_ex[, `:=`(assay_cat = ifelse(Abbott_Architect_perc == 100, 'Abbott',
                                      ifelse(Roche_Elecsys_perc == 100, 'Roche', 
                                             ifelse(Ortho_VITROS_perc == 100, 
                                                    'Ortho', 'Abbott & Ortho'))))]

map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                by = c('GEOID' = 'state_code'))
map = left_join(map, s_ex, by = 'state')
map$round = paste('Round', map$survey_round)
map$round = factor(map$survey_round, levels = c('Round 4', 'Round 13', 'Round 24'))
map = subset(map, state != 'DC')

cols = qualitative_hcl(4, palette = 'Dark')

pl1 = ggplot(map, aes(fill = assay_cat)) +
  geom_sf(colour = 'black', size = 0.01) +
  scale_fill_manual(values = cols, na.value = 'gray50', 
                    limits = c('Abbott', 'Abbott & Ortho', 'Ortho', 'Roche')) +
  facet_wrap(~ survey_round, nrow = 1) +
  labs(fill = '') +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.8, 'cm'),
        plot.title = element_text(hjust = 0.5))

pl2 = ggplot(map, aes(fill = seroprevalence / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE) +
  facet_wrap(~ survey_round, nrow = 1) +
  labs(fill = 'CDC seroprevalence') +
  theme_void() +
  theme(legend.position = 'top',
        legend.title = element_text(size = 9))


s_sum = s[, .(week = mean(week),
              Abbott = mean(Abbott_Architect_perc), 
              Ortho = mean(Ortho_VITROS_perc),
              Roche = mean(Roche_Elecsys_perc)), 
          by = 'survey_round']
s_sum = melt(s_sum, id.vars = c('week', 'survey_round'), variable.name = 'Assay')

pl3 = ggplot(s_sum, aes(x = week, y = value, colour = Assay)) +
  geom_line(show.legend = FALSE) + geom_point(size = 0.6, show.legend = FALSE) +
  scale_colour_manual(values = cols[-2]) +
  geom_vline(data = subset(s_sum, survey_round %in% c(4, 13, 24)), 
             aes(xintercept = week), alpha = 0.05, size = 3) +
  scale_x_date(date_breaks = '6 months', date_labels = '%b %Y') +
  labs(x = NULL, y = 'Mean % use of assay') +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0.3, 0, 0), 'cm'))

s_glm = GetPredsPerAssay(d = states_w, ex = best[[1]], f = 'f_glm')
s_all_glm = GetPredTS(dat = s_glm)

pl4 = ggplot(subset(s_all_glm, assay_lab == 'Surveys'), 
             aes(x = week, y = values)) + 
  geom_line() +
  geom_point(size = 0.7) +
  geom_vline(data = subset(s_all_glm, survey_round %in% c(4, 13, 24)), 
             aes(xintercept = week), alpha = 0.05, size = 3) +
  scale_x_date(date_breaks = '6 months', date_labels = '%b %Y') +
  xlab(NULL) + ylab('Seroprevalence') +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0.3, 0, 0), 'cm'))

pl = ggarrange(pl1, pl2, ggarrange(pl3, pl4, nrow = 1), ncol = 1, heights = c(1, 1.12, 0.8))

ggsave('../plots/figures1.pdf', width = 6, height = 5, unit = 'in')




###############################################################################
# SI Fig S2: %assay vs. seroprevalence.
###############################################################################

#' Jitter points on plot.

Jit = function(x)
{
  w_x = which(x == 100 | x == 0)
  x[w_x] = x[w_x] + runif(length(w_x), -2, 2)
  return(x)
}


pdf(file = '../plots/figures2.pdf', width = 5, height = 1.5, pointsize = 8)

par(mar = c(2.2, 1.6, 0.5, 0.5), mgp = c(1.2, 0.3, 0), tcl = 0.15, las = 1)

split.screen(c(1, 3))

par(oma = c(0, 1.4, 0, 0))

screen(1)
plot(Jit(s$Abbott_Architect_perc), s$seroprevalence / 100, pch = 16, 
     cex = 0.6, col = adjustcolor('black', alpha = 0.2),
     xlab = '% Abbott', ylab = '',
     xlim = c(-3, 103))
mtext('(a)', side = 3, adj = 0.90, line = -1.2)

screen(2)
plot(Jit(s$Ortho_VITROS_perc), s$seroprevalence / 100, pch = 16, 
     cex = 0.6, col = adjustcolor('black', alpha = 0.2),
     xlab = '% Ortho', ylab = '',
     xlim = c(-3, 103))
mtext('(b)', side = 3, adj = 0.90, line = -1.2)

screen(3)
plot(Jit(s$Roche_Elecsys_perc), s$seroprevalence / 100, pch = 16, 
     cex = 0.6, col = adjustcolor('black', alpha = 0.2),
     xlab = '% Roche', ylab = '',
     xlim = c(-3, 103))
mtext('(c)', side = 3, adj = 0.90, line = -1.2)

close.screen(all.screens = TRUE)

par(las = 0)

mtext(text = 'Seroprevalence', side = 2, outer = TRUE, line = -1)

dev.off()



###############################################################################
# SI Fig S3: Map of assays used.
###############################################################################

round_max = 25
s_ex = subset(s, survey_round < (round_max + 1))

s_ex = s_ex[, `:=`(assay_cat = ifelse(Abbott_Architect_perc == 100, 'Abbott',
                                      ifelse(Roche_Elecsys_perc == 100, 'Roche', 
                                             ifelse(Ortho_VITROS_perc == 100, 
                                                    'Ortho', 'Abbott & Ortho'))))]
s_ex = subset(s_ex, survey_round <= 15 | survey_round == 24)
s_ex$survey_round[which(s_ex$survey_round == 15)] = '15-23'

map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                by = c('GEOID' = 'state_code'))
map = left_join(map, s_ex, by = 'state')
map = subset(map, !state %in% c('DC', 'ND'))

map$survey_round = factor(map$survey_round, levels = unique(map$survey_round))

cols = qualitative_hcl(4, palette = 'Dark')

pl = ggplot(map, aes(fill = assay_cat)) +
  geom_sf(colour = 'black', size = 0.1) +
  scale_fill_manual(values = cols, na.value = 'gray50', 
                    limits = c('Abbott', 'Abbott & Ortho', 'Ortho', 'Roche')) +
  facet_wrap(~ survey_round, ncol = 4) +
  labs(fill = 'Assay used') +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.8, 'cm'),
        plot.title = element_text(hjust = 0.5))

ggsave('../plots/figures3.pdf', width = 8.5, height = 6.5)



###############################################################################
# SI Figs S4-S6: Comparing seroprevalence, incidence, vaccination for all
# rounds.
###############################################################################

map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                by = c('GEOID' = 'state_code'))
map = left_join(map, subset(preds_glm, survey_round < 11), by = 'state')
map = map[!map$state %in% c('ND', 'DC'), ]

pl1 = ggplot(map, aes(fill = seroprevalence / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                   lim = c(0, 0.7)) +
  labs(fill = 'Observed\nseroprev.') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl2 = ggplot(map, aes(fill = pred_incidence)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                   lim = c(0, 0.7)) +
  labs(fill = 'Estimated\nproportion\ninfected') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl3 = ggplot(map, aes(fill = pred_incidence - seroprevalence / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_diverging(palette = 'Blue-Red 2', rev = TRUE) +
  labs(fill = 'Infections -\nseroprev.') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl4 = ggplot(map, aes(fill = vaccinated_cumulative_perc / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Blues', rev = TRUE, lim = c(0, 0.79)) +
  labs(fill = 'Proportion\nvaccinated') +
  facet_wrap(~ survey_round, ncol = 1, strip.position = 'right') +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15))

pl = ggarrange(pl1, pl2, pl3, pl4, ncol = 4,
               font.label = list(face = 'plain'),
               common.legend = FALSE) 

ggsave('../plots/figures4.png', width = 10, height = 14, unit = 'in')



map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                by = c('GEOID' = 'state_code'))
map = left_join(map, subset(preds_glm, survey_round > 10 & survey_round < 21), 
                by = 'state')
map = map[!map$state %in% c('ND', 'DC'), ]

pl1 = ggplot(map, aes(fill = seroprevalence / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                   lim = c(0, 0.7)) +
  labs(fill = 'Observed\nseroprev.') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl2 = ggplot(map, aes(fill = pred_incidence)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                   lim = c(0, 0.7)) +
  labs(fill = 'Estimated\nproportion\ninfected') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl3 = ggplot(map, aes(fill = pred_incidence - seroprevalence / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_diverging(palette = 'Blue-Red 2', rev = TRUE) +
  labs(fill = 'Infections -\nseroprev.') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl4 = ggplot(map, aes(fill = vaccinated_cumulative_perc / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Blues', rev = TRUE, lim = c(0, 0.79)) +
  labs(fill = 'Proportion\nvaccinated') +
  facet_wrap(~ survey_round, ncol = 1, strip.position = 'right') +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15))

pl = ggarrange(pl1, pl2, pl3, pl4, ncol = 4,
               font.label = list(face = 'plain'),
               common.legend = FALSE) 

ggsave('../plots/figures5.png', width = 10, height = 14, unit = 'in')



map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                by = c('GEOID' = 'state_code'))
map = left_join(map, subset(preds_glm, survey_round > 20), by = 'state')
map = map[!map$state %in% c('ND', 'DC'), ]

pl1 = ggplot(map, aes(fill = seroprevalence / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                   lim = c(0, 0.7)) +
  labs(fill = 'Observed\nseroprev.') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl2 = ggplot(map, aes(fill = pred_incidence)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Greens', rev = TRUE, 
                                   lim = c(0, 0.7)) +
  labs(fill = 'Estimated\nproportion\ninfected') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl3 = ggplot(map, aes(fill = pred_incidence - seroprevalence / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_diverging(palette = 'Blue-Red 2', rev = TRUE) +
  labs(fill = 'Infections -\nseroprev.') +
  facet_wrap(~ survey_round, ncol = 1) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())

pl4 = ggplot(map, aes(fill = vaccinated_cumulative_perc / 100)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Blues', rev = TRUE, lim = c(0, 0.79)) +
  labs(fill = 'Proportion\nvaccinated') +
  facet_wrap(~ survey_round, ncol = 1, strip.position = 'right') +
  theme_void() +
  theme(legend.position = 'top', 
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15))

pl = ggarrange(pl1, pl2, pl3, pl4, ncol = 4,
               font.label = list(face = 'plain'),
               common.legend = FALSE) 

ggsave('../plots/figures6.png', width = 10, height = 13.5, unit = 'in')




###############################################################################
# SI Fig S7: Time series of overall seroprevalence adjusting by assay.
###############################################################################

#' Plot maps and time series of predicted seroprevalences per assay.
#'
#' @param data_map :data.table: data to use in maps (result of using
#' GetPredsPerAssay above).
#' @param data_ts :data.table: data to use in time series (result of using
#' GetPredTS above).
#'
#' @return A plot with the change in seroprevalence had each of the three assays
#' in turn been used exclusively, at four different points in time (3 x 4 maps),
#' and a timeseries with overall seroprevalence in the US (including what
#' seroprevalence would have been had a single assay been used exclusively).

PlotMapAndTs = function(data_map, data_ts)
{
  data(state_laea)
  data(fips_codes)

  rounds = seq(4, 28, 8)

  prev_map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                       by = c('GEOID' = 'state_code'))
  prev_map = left_join(prev_map, subset(data_map, survey_round %in% rounds), 
                       by = 'state')
  prev_map = prev_map[!prev_map$state %in% c('ND', 'DC'), ]
  prev_map$survey_round = factor(prev_map$round_pl, 
                                 levels = paste('Round', rounds))

  pl_m = ggplot(prev_map, aes(fill = preds_diff)) +
    geom_sf(colour = NA) +
    scale_fill_continuous_diverging(palette = 'Blue-Red', rev = TRUE) +
    labs(fill = '') +
    theme_void() +
    facet_grid(assay_preds ~ survey_round, drop = TRUE) +
    theme(legend.position = 'top', 
          legend.key.width = unit(1.5, 'cm'),
          plot.title = element_text(hjust = 0.5))

  cols = c('#000000', brewer.pal(3, 'Set1'))

  pl_ts = ggplot(data_ts, aes(x = week, y = values, colour = assay_lab)) + 
    geom_line() +
    geom_point() +
    scale_colour_manual(values = cols) +
    geom_vline(data = subset(data_ts, survey_round %in% rounds), 
               aes(xintercept = week), alpha = 0.05, size = 5) +
    xlab('') + ylab('Seroprevalence') +
    scale_x_date(date_breaks = '3 months', date_labels = '%b %Y') +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.06, 0.72),
          legend.background = element_rect(NA)) 

  pl = ggarrange(pl_m, pl_ts, nrow = 2, ncol = 1, heights = c(2.4, 1),
                 font.label = list(face = 'plain'),
                 common.legend = FALSE)

  ggsave(paste0('../plots/figures7.pdf'), width = 9, height = 8, unit = 'in')
}


s_glm = GetPredsPerAssay(d = states_w, ex = best[[1]], f = 'f_glm')
s_all_glm = GetPredTS(dat = s_glm)
PlotMapAndTs(data_map = s_glm, data_ts = s_all_glm)



################################################################################
# Fig S8: Complete version of model metrics tile plots.
################################################################################

# The metric value corresponding to the bottom 5 percentile, for each metric.
w_b_aic = quantile(out_3w_glm$aic, prob = 0.05) - m_ref_glm$aic
w_b_rms = quantile(out_3w_glm$rmse, prob = 0.05) - m_ref_glm$rmse
w_b_loo = quantile(out_3w_glm$rmse_loov, prob = 0.05) - m_ref_glm$rmse_loov

out_3w_glm_long = out_3w_glm
# For plot panel text.
w = which(abs(out_3w_glm_long$detects) == 1)
out_3w_glm_long$detects[w]  = paste(out_3w_glm_long$detects[w], 'week')
out_3w_glm_long$detects[-w] = paste(out_3w_glm_long$detects[-w], 'weeks')

# Express metrics relative to that of the main reference model.
out_3w_glm_long$aic       = out_3w_glm_long$aic - m_ref_glm$aic
out_3w_glm_long$rmse      = out_3w_glm_long$rmse - m_ref_glm$rmse
out_3w_glm_long$rmse_loov = out_3w_glm_long$rmse_loov - m_ref_glm$rmse_loov

out_3w_glm_long = pivot_longer(out_3w_glm_long,
                               cols = c('aic', 'rmse', 'rmse_loov'),
                               names_to = 'metric', )
# Subset to point estimates for the Ortho assay, per each metric
out_3w_glm_long = subset(out_3w_glm_long,
                         (metric == 'aic' & waning_wks_ortho == 49) | 
                         (metric == 'rmse' & waning_wks_ortho == 49) | 
                         (metric == 'rmse_loov' & waning_wks_ortho == 10))

cols = c('detects', 'waning_wks_abbott', 'waning_wks_ortho')
b_aic       = subset(out_3w_glm_long, metric == 'aic')
b_aic       = b_aic[which.min(b_aic$value), ]
b_rmse      = subset(out_3w_glm_long, metric == 'rmse')
b_rmse      = b_rmse[which.min(b_rmse$value), ]
b_rmse_loov = subset(out_3w_glm_long, metric == 'rmse_loov')
b_rmse_loov = b_rmse_loov[which.min(b_rmse_loov$value), ]


p1 = ggplot(subset(out_3w_glm_long, metric == 'aic'), 
            aes(x = waning_wks_roche, y = waning_wks_abbott)) + 
  geom_tile(aes(fill = value)) +
  geom_line(data = data.frame(x = c(7.5, 70.5), y = c(7.5, 70.5)), 
            aes(x = x, y = y), colour = 'white', linetype = '33') +
  geom_point(data = b_aic, size = 2, colour = 'green',
             aes(x = waning_wks_roche, y = waning_wks_abbott, fill = NULL)) + 
  geom_contour(aes(z = value), breaks = w_b_aic, col = 'black', alpha = 0.6) +
  coord_fixed() +
  scale_fill_continuous_diverging(palette = 'Blue-Red 3') +
  labs(x = NULL, 
       y = NULL,
       fill = 'AIC') +
  facet_grid(~ detects) +
  theme_bw() +
  theme(strip.background = element_blank())

p2 = ggplot(subset(out_3w_glm_long, metric == 'rmse'),
            aes(x = waning_wks_roche, y = waning_wks_abbott)) + 
  geom_tile(aes(fill = value)) +
  geom_line(data = data.frame(x = c(7.5, 70.5), y = c(7.5, 70.5)), 
            aes(x = x, y = y), colour = 'white', linetype = '33') +
  geom_point(data = b_rmse, size = 2, colour = 'green',
             aes(x = waning_wks_roche, y = waning_wks_abbott, fill = NULL)) + 
  geom_contour(aes(z = value), breaks = w_b_rms, col = 'black', alpha = 0.6) +
  coord_fixed() +
  scale_fill_continuous_diverging(palette = 'Blue-Red 3') +
  labs(x = NULL, 
       y = NULL,
       fill = 'RMSE') +
  facet_grid(~ detects) +
  theme_bw() +
  theme(strip.background = element_blank())

p3 = ggplot(subset(out_3w_glm_long, metric == 'rmse_loov'),
            aes(x = waning_wks_roche, y = waning_wks_abbott)) + 
  geom_tile() +
  geom_tile(aes(fill = value)) +
  geom_line(data = data.frame(x = c(7.5, 70.5), y = c(7.5, 70.5)), 
            aes(x = x, y = y), colour = 'white', linetype = '33') +
  geom_point(data = b_rmse_loov, size = 2, colour = 'green',
             aes(x = waning_wks_roche, y = waning_wks_abbott, fill = NULL)) + 
  geom_contour(aes(z = value), breaks = w_b_loov, col = 'black', alpha = 0.6) +
  coord_fixed() +
  scale_fill_continuous_diverging(palette = 'Blue-Red 3') +
  labs(x = NULL, 
       y = NULL,
       fill = 'LOO RMSE') +
  facet_grid(~ detects) +
  theme_bw() +
  theme(strip.background = element_blank())

p = ggarrange(p1, p2, p3, ncol = 1, align = 'v')
p = annotate_figure(p, bottom = text_grob('Roche time to seroreversion (weeks)'),
                    left = text_grob('Abbott time to seroreversion (weeks)', 
                                     rot = 90, vjust = 1))

ggsave('../plots/figures8.pdf', width = 8, height = 5, unit = 'in')
              


###############################################################################
# SI Fig S9. Estimated proportion infected, vaccination coverage, per state.
###############################################################################

preds_glm_all = subset(preds_glm_all, state != 'ND' & 
                       week >= as.Date('2020-08-01') &
                       week <= as.Date('2022-01-15'))

pl = ggplot(subset(preds_glm, state != 'US-wide'), 
            aes(x = week, y = seroprevalence / 100)) +
  geom_ribbon(aes(x = week, ymin = pred_seroprev_lo, ymax = pred_seroprev_hi),
              alpha = 0.2) +
  geom_ribbon(aes(ymin = pred_incidence_lo, ymax = pred_incidence_hi), 
              alpha = 0.2, fill = 'royalblue2') +
  geom_point() +
  geom_line(aes(y = pred_seroprev, colour = 'Fitted seroprevalence')) +
  geom_line(data = preds_glm_all,
            aes(y = pred_incidence, colour = 'Estimated proportion infected')) +
  geom_line(aes(y = vaccinated_cumulative_perc / 100, 
                colour = 'Proportion vaccinated')) +
  facet_wrap(~ state, ncol = 5) +
  scale_x_date(date_breaks = '6 months', date_labels = '%b %Y') +
  theme_bw() +
  labs(x = NULL, y = 'Seroprevalence, estimated proportion infected, vaccination coverage') +
  scale_colour_manual('', 
                      breaks = c('Fitted seroprevalence',
                                 'Estimated proportion infected',
                                 'Proportion vaccinated'),
                      values = c('Fitted seroprevalence' = 'darkgray',
                                 'Estimated proportion infected' = 'royalblue2',
                                 'Proportion vaccinated' = 'red')) +
  guides(colour = guide_legend(byrow = TRUE)) + 
  theme(legend.position = 'bottom', 
        strip.background = element_blank(),
        legend.key.width = unit(1, 'cm'),
        axis.text = element_text(size = 7))

ggsave('../plots/figures9.pdf', width = 8, height = 12, unit = 'in')



################################################################################
# SI Fig S10. Variation in seroprevalences and incidences across the US, per round.
################################################################################

#' Calculate and plot variation for seroprevalences and predicted incidences
#' across the US, per round.
#'
#' @description The coefficients of variation should indicate the spatial
#' heterogeneity in seroprevalences and incidences across the US for each point
#' in time.
#'
#' @param preds :data.table: augmented data.table with predicted values (output
#' from GetPredictions above).
#' @param ex :list: list of length 4 numeric vectors, with the four elements
#' being the times to seroreversion for Abbott, Ortho, and Roche assays, and the
#' lead/lag time.
#'
#' @return A plot with two panels, one for observed seroprevalences, and another
#' for estimated proportions infected (both panels with and without
#' vaccinations.

PlotVars = function(preds, ex)
{
  # COV for seroprevalence without vaccination.
  cov_med_s  = preds[, .(par = 'No vacc.', 
                         week = mean(week),
                         cov_med = sd(seroprevalence / 100) / 
                           mean(seroprevalence / 100),
                         cov_med_neg = NA), 
                     by = 'survey_round']

  # COV for seroprevalence including vaccination (with both no assumed
  # correlation, and negative correlation).
  cov_med_sv = preds[, .(par = 'With vacc.', 
                         week = mean(week),
                         cov_med = sd(seroprev_v) / mean(seroprev_v),
                         cov_med_neg = sd(seroprev_v_neg) / mean(seroprev_v_neg)), 
                     by = 'survey_round']

  cov_med_s = rbind(cov_med_s, cov_med_sv)
  cov_med_s = cov_med_s[, par_facet := 'Observed seroprevalence']

  # COV for estimated proportion infected without vaccination.
  cov_med_s_ex  = preds[, .(par = 'No vacc.', 
                            week = mean(week),
                            cov_med = sd(pred_incidence) / mean(pred_incidence),
                            cov_med_neg = NA), 
                        by = c('survey_round')]

  # COV for EPIV (with both no assumed correlation, and negative correlation).
  cov_med_sv_ex = preds[, .(par = 'With vacc.', 
                            week = mean(week),
                            cov_med = sd(pred_epiv) / mean(pred_epiv),
                            cov_med_neg = sd(pred_epiv_neg) /
                              mean(pred_epiv_neg)), 
                        by = c('survey_round')]

  cov_med_s_ex = rbind(cov_med_s_ex, cov_med_sv_ex)
  cov_med_s_ex = cov_med_s_ex[, `:=`(par_facet = 'Estimated proportion infected')]

  cov_meds = rbind(cov_med_s, cov_med_s_ex)

  cov_meds$par_facet = factor(cov_meds$par_facet, 
                              levels = c('Observed seroprevalence', 
                                         'Estimated proportion infected'))
  b = brewer.pal(3, 'Set1')

  p = ggplot(cov_meds, aes(x = week, y = cov_med, colour = par)) + 
    geom_line(show.legend = TRUE, aes(fill = 'Ind. prob.')) +
    geom_line(data = cov_meds, aes(y = cov_med_neg, colour = par, 
                               fill = 'Neg. corr.'), linetype = '33') +
    geom_line(data = subset(cov_meds, par == 'No vacc.'), aes(y = cov_med)) + 
    scale_colour_brewer(palette = 'Set1', guide = guide_legend(order = 1)) +
    scale_fill_manual('', values = c(1:2),
                      guide = guide_legend(override.aes = list(colour = b[2], 
                                                               linetype = c(1, 2)))) +
    labs(colour = '') +
    scale_x_date(date_breaks = '6 months', date_labels = '%b %Y') +
    xlab(NULL) + ylab('Coefficient of variation') +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.spacing = unit(1, 'lines')) +
    facet_wrap(~ par_facet, ncol = 2)
}

pl_cov_meds = PlotVars(preds = preds_glm, ex = best)
ggsave('../plots/figures10.pdf', width = 6, height = 2.5, unit = 'in')



################################################################################
# Fig S11: Maps of pre-infection vaccination coverage at two time points.
################################################################################

preds_glm_all = GetPredictions(d = states_w, ex = best, f = 'f_glm_3w')

max_date = max(subset(preds_glm_all, survey_round == 29)$week)
preds_glm_all = subset(preds_glm_all, week <= max_date)

preds_glm_all = preds_glm_all[, vacc_d := c(NA, diff(vaccinated_cumulative_perc / 100)), 
                              by = 'state']

preds_glm_all = preds_glm_all[, pred_v_lag := shift(pred_epiv, n = 1, fill = 0)]
preds_glm_all = preds_glm_all[, metric_pre := vacc_d * (1 - pred_v_lag)]


#' Get cumulative sum excluding NAs.

CumSumNoNA = function(x)
{
  cumsum(ifelse(is.na(x), 0, x)) + x * 0
}

preds_glm_all = preds_glm_all[, metric := CumSumNoNA(metric_pre), by = 'state']

sub_date = subset(preds_glm_all, week %in% c(as.Date('2021-03-20'), max_date))
sub_date = rbind(sub_date, sub_date[(nrow(sub_date) + 1):(nrow(sub_date) + 2), ])

sub_date$state[(nrow(sub_date) - 1):nrow(sub_date)] = 'ND'
sub_date$week[(nrow(sub_date) - 1):nrow(sub_date)]  = c(as.Date('2021-03-20'), max_date)


map = left_join(state_laea, unique(fips_codes[, c('state', 'state_code')]), 
                by = c('GEOID' = 'state_code'))
map = left_join(map, sub_date, by = 'state')
map = map[!map$state %in% c('DC'), ]


pl = ggplot(map, aes(fill = metric)) +
  geom_sf(colour = NA) +
  scale_fill_continuous_sequential(palette = 'Oranges', rev = TRUE) +
  labs(fill = NULL) +
  facet_wrap(~ week) +
  theme_void() +
  theme(legend.position = 'top', 
        legend.title = element_text(size = 9),
        legend.key.width = unit(0.8, 'cm'),
        plot.title = element_text(hjust = 0.5))

ggsave('../plots/figures11.pdf', width = 7, height = 3.2, unit = 'in')



################################################################################
# SI Fig S12 (version of Fig 5 but looking at anti-N seroprevalence in the blood
# donors data). 
################################################################################

l = loess((seroprevalence / 100) ~ blood_prevalence_no_vacc, data = plot_preds)
plot_preds$serosurvey_blood[!is.na(plot_preds$blood_prevalence_no_vacc)] = predict(l)

l = loess(pred_incidence ~ blood_prevalence_no_vacc, data = plot_preds)
plot_preds$pred_inc_blood[!is.na(plot_preds$blood_prevalence_no_vacc)] = predict(l)

p1 = ggplot(plot_preds, aes(x = blood_prevalence_no_vacc, y = serosurvey_blood)) +
  geom_abline(intercept = 0, slope = 1, lty = '33') +
  geom_point(aes(x = blood_prevalence_no_vacc, y = seroprevalence / 100), 
             colour = 'black', alpha = 0.2, size = 0.6) +
  coord_fixed() + 
  lims(x = c(0, 0.6), y = c(0, 0.6)) +
  geom_line(colour = 'black', size = 1) + 
  labs(x = NULL, 
       y = 'Seroprevalence from CDC serosurveys') +
  theme_bw() 

p2 = ggplot(plot_preds, aes(x = blood_prevalence_no_vacc, y = pred_inc_blood)) +
  geom_abline(intercept = 0, slope = 1, lty = '33') +
  geom_point(aes(x = blood_prevalence_no_vacc, y = pred_incidence), 
             colour = 'royalblue2',
             alpha = 0.2, size = 0.6) +
  coord_fixed() + 
  lims(x = c(0, 0.6), y = c(0, 0.6)) +
  geom_line(colour = 'royalblue2', size = 1) + 
  labs(x = NULL, 
       y = 'Estimated proportion infected') +
  theme_bw() 

pl = ggarrange(p1, p2, nrow = 1, ncol = 2)
pl = annotate_figure(pl, 
                     left = text_grob(NULL),
                     bottom = text_grob('Seroprevalence (no vaccinations) from the blood donors dataset'))

ggsave('../plots/figures12.pdf', width = 6, height = 3.2)



################################################################################
# SI Fig S13: Time series of blood donor prevalence compared to predictions.
################################################################################

preds_glm = GetPredictions(d = states_w, ex = best, f = 'f_glm_3w')
preds_glm = subset(preds_glm, state != 'ND' &
                   week >= as.Date('2020-08-01') &
                   week <= as.Date('2022-01-15'))

b = brewer.pal(3, 'Set1')

don_sub = subset(don, !state %in% c('Al', 'CR', 'ND', 'PR'))

p = ggplot(na.omit(don_sub[, c('state', 'week', 'blood_prevalence')])) +
  geom_line(aes(x = week, y = blood_prevalence, 
                colour = 'Anti-S seroprevalence (blood donors)'), 
            size = 0.8) +
  geom_line(data = preds_glm,
            aes(x = week, y = pred_epiv, 
                colour = 'EPIV (independent prob., complete series)'),
            size = 0.8) +
  geom_line(data = preds_glm,
            aes(x = week, y = pred_epiv_d1, 
                colour = 'EPIV (independent prob., dose 1)'), 
            size = 0.8) + 
  geom_line(data = preds_glm,
            aes(x = week, y = pred_epiv_d1_neg, 
                colour = 'EPIV (negative correlation, dose 1)'), 
            size = 0.8) + 
  facet_wrap(~ state, ncol = 5) + 
  labs(x = NULL, y = 'Seroprevalence') +
  scale_colour_manual('', 
                      breaks = c('Anti-S seroprevalence (blood donors)',
                                 'EPIV (independent prob., complete series)',
                                 'EPIV (independent prob., dose 1)',
                                 'EPIV (negative correlation, dose 1)'),
                      values = c('Anti-S seroprevalence (blood donors)' = 'black',
                                 'EPIV (independent prob., complete series)' = b[2],
                                 'EPIV (independent prob., dose 1)' = b[1],
                                 'EPIV (negative correlation, dose 1)' = b[3])) +
  guides(colour = guide_legend(ncol = 2)) +
  scale_x_date(date_breaks = '6 months', date_labels = '%b %Y') +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-15, -15, -5, -15),
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 90))

ggsave('../plots/figures13.pdf', width = 8, height = 12)



###############################################################################
# SI Fig S14: Correlation matrices.
###############################################################################

s = states_w[!is.na(seroprevalence), ]
s = subset(s, state != 'ND')
s = s[, `:=`(sqrt_cases_cumulative_perc = sqrt(cases_cumulative_perc),
             sqrt_deaths_cumulative_perc = sqrt(deaths_cumulative_perc),
             ln_tests_total_cumulative_perc = log(tests_total_cumulative_perc),
             sqrt_hospitalisation_cumulative_perc = sqrt(hospitalisation_cumulative_perc))]

vars = c('sqrt_cases_cumulative_perc', 'sqrt_deaths_cumulative_perc', 
         'excess_deaths_difference_perc',
         'ln_tests_total_cumulative_perc',  
         'vaccinated_cumulative_perc', 
         'sqrt_hospitalisation_cumulative_perc',
         'cases_00_19_perc', 'cases_20_49_perc', 'cases_50_69_perc',
         'cases_70_perc', 
         'Abbott_Architect_perc', 'Ortho_VITROS_perc', 'Roche_Elecsys_perc')
labs = c('Sqrt % cases', 'Sqrt % deaths', '% exc. deaths diff.', 
         'Ln % tested', '% vaccinated', 'Sqrt % hosp.', 
         '% 0-19', '% 20-49', '% 50-69', '% >70', '% Abbott', '% Ortho', '% Roche')
cors = cor(unique(s[, ..vars]), method = 'spearman')
cors_test = cor.mtest(unique(s[, ..vars]), method = 'spearman', 
                      conf.level = 0.95, exact = TRUE)
colnames(cors) = labs
rownames(cors) = labs

pdf(file = '../plots/figures14.pdf', width = 5.5, height = 5.5, pointsize = 8)

par(mar = c(0.5, 0.5, 0.5, 0.5), mgp = c(1.2, 0.2, 0), las = 1, tcl = 0.25)

corrplot(cors, p.mat = cors_test$p, 
         type = 'upper', order = 'original', 
         tl.pos = 'lt', tl.col = 'black', insig = 'blank',  sig.level = 0.05)
corrplot(cors, p.mat = cors_test$p, 
         type = "lower", method = "number", order = 'original',
         tl.pos = 'n', diag = FALSE, cl.pos = 'n', insig = 'blank', add = TRUE)
dev.off()



################################################################################
# SI Fig S15: diagrammatic example of different waning rates and how adjustment
# of cases works.
################################################################################

pdf(file = '../plots/figures15.pdf', width = 5, height = 4, pointsize = 8)

par(mar = c(2.5, 3, 0.5, 0.5), mgp = c(1.6, 0.3, 0), tcl = 0.14, las = 1)

split.screen(c(3, 3))

screen(1)
plot(1, 1, type = 'n', xlim = c(0, 43), ylim = c(-0.08, 1.22), yaxt = 'n',
     xlab = 'Weeks', ylab = 'Prob. seropositive')
lines(c(0, 2, 2, 16, 16, 50), c(0, 0, 1, 1, 0, 0), lty = '12', lwd = 0.7)
axis(2, at = c(0, 0.5, 1))
axis(4, at = c(0, 0.5, 1), labels = FALSE)

rect(xleft = 0, xright = 2, ybottom = 1.06, ytop = 1.14, border = NA, 
     col = adjustcolor('forestgreen', alpha = 0.4))
rect(xleft = 2, xright = 16, ybottom = 1.06, ytop = 1.14, border = NA, 
     col = adjustcolor('darkred', alpha = 0.4))

screen(4)
plot(1, 1, type = 'n', xlim = c(0, 43), ylim = c(-0.08, 1.22), yaxt = 'n',
     xlab = 'Weeks', ylab = 'Prob. seropositive')
lines(c(0, 2, 2, 24, 24, 50), c(0, 0, 1, 1, 0, 0), lty = '12', lwd = 0.7)
axis(2, at = c(0, 0.5, 1))
axis(4, at = c(0, 0.5, 1), labels = FALSE)

rect(xleft = 0, xright = 2, ybottom = 1.06, ytop = 1.14, border = NA, 
     col = adjustcolor('forestgreen', alpha = 0.4))
rect(xleft = 2, xright = 24, ybottom = 1.06, ytop = 1.14, border = NA, 
     col = adjustcolor('darkred', alpha = 0.4))

screen(7)
plot(1, 1, type = 'n', xlim = c(0, 43), ylim = c(-0.08, 1.22), yaxt = 'n',
     xlab = 'Weeks', ylab = 'Prob. seropositive')
lines(c(0, 2, 2, 40, 40, 50), c(0, 0, 1, 1, 0, 0), lty = '12', lwd = 0.7)
axis(2, at = c(0, 0.5, 1))
axis(4, at = c(0, 0.5, 1), labels = FALSE)

rect(xleft = 0, xright = 2, ybottom = 1.06, ytop = 1.14, border = NA, 
     col = adjustcolor('forestgreen', alpha = 0.4))
rect(xleft = 2, xright = 40, ybottom = 1.06, ytop = 1.14, border = NA, 
     col = adjustcolor('darkred', alpha = 0.4))


# Get three examples of waning: slow, medium, fast, to illustrate the resulting
# differences of the cases with waning, actual cumulative cases, and what that
# might be expected to imply with regards to the distinction between
# seroprevalence with waning and without.
# 
# Use Florida as the basis for the example.
states_w = states_w[, `:=`(ex_1 = GetWanedCasesSimple(cases, detect = 0, 
                                                      ab_dur = 40),
                           ex_2 = GetWanedCasesSimple(cases, detect = 0, 
                                                      ab_dur = 24),
                           ex_3 = GetWanedCasesSimple(cases, detect = 0, 
                                                      ab_dur = 16)),
                    by = state]

fl = subset(states_w, state == 'FL')

screen(2)
plot(head(fl$ex_3, 78) / 1e6, type = 'l', lwd = 0.7, lty = '12',
     ylim = c(0, 2.4), xlab = 'Week', ylab = 'Cases (millions)')
lines(fl$cases_cumulative / 1e6, lwd = 0.7)

screen(5)
plot(head(fl$ex_2, 78) / 1e6, type = 'l', lwd = 0.7, lty = '12',
     ylim = c(0, 2.4), xlab = 'Week', ylab = 'Cases (millions)')
lines(fl$cases_cumulative / 1e6, lwd = 0.7)

screen(8)
plot(head(fl$ex_1, 78) / 1e6, type = 'l', lwd = 0.7, lty = '12',
     ylim = c(0, 2.4), xlab = 'Week', ylab = 'Cases (millions)')
lines(fl$cases_cumulative / 1e6, lwd = 0.7)

MakeCumulative = function(x)
{
  x1 = diff(x)
  x1[x1 < 0] = 0
  return(cumsum(c(x[1], x1)))
}

fl = subset(fl, !is.na(seroprevalence) & survey_round < 23)

screen(3)
plot(fl$seroprevalence / 100 * 0.85, type = 'l', lwd = 0.7, lty = '12', 
     ylim = c(0, 1),
     xlab = 'Round', ylab = 'Seroprevalence')
y = MakeCumulative(fl$cases_cumulative / fl$ex_3 * fl$seroprevalence / 100 * 0.85)
lines(y, lwd = 0.7)

screen(6)
plot(fl$seroprevalence / 100 * 0.85, type = 'l', lwd = 0.7, 
     lty = '12', ylim = c(0, 1),
     xlab = 'Round', ylab = 'Seroprevalence')
y = MakeCumulative(fl$cases_cumulative / fl$ex_2 * fl$seroprevalence / 100 * 0.85)
lines(y, lwd = 0.7)

screen(9)
plot(fl$seroprevalence / 100 * 0.85, type = 'l', lwd = 0.7, 
     lty = '12', ylim = c(0, 1),
     xlab = 'Round', ylab = 'Seroprevalence')
y = MakeCumulative(fl$cases_cumulative / fl$ex_1 * fl$seroprevalence / 100 * 0.85)
lines(y, lwd = 0.7)

close.screen(all.screens = TRUE)
dev.off()



################################################################################
# SI Fig S16: GAM functions.
################################################################################

#' Get the GAM functions for each variable from a model fit.
#'
#' @param m :GAM object: fit GAM model.
#' @param p :character: variable to extract.
#'
#' @return a matrix with four columns: x values, estimate, and lower and upper
#' bounds.

GetGamFuns = function(m, p)
{
  n = strsplit(as.character(m$formula)[3], '\\+')[[1]]
  n = n[setdiff(seq_along(n), grep('state', n))]
  n = n[setdiff(seq_along(n), grep('s\\(', n, invert = TRUE))]
  g = grep(p, n)[1]

  if (!is.na(g))
  {
    z = plot(m, select = 1, ylab = '', se = TRUE)
    x = cbind(z[[g]]$x, z[[g]]$fit, 
              z[[g]]$fit - z[[g]]$se, z[[g]]$fit + z[[g]]$se)
  } else
  {
    x = data.frame(NA, NA)
  }

  return(x)
}



#' Plot GAM functions extracted from GetGamFuns

PlotGamFuns = function(x, p, lab, add = FALSE, cols = 'black')
{
  xr = range(x[, 1], na.rm = TRUE)
  xr = range(xr[is.finite(xr)])

  yr = c(-2, 1.5)

  if (isFALSE(add))
  {
    plot(1, 1, type = 'n', xlim = xr, ylim = yr, xlab = lab, ylab = '')
    abline(h = 0, lty = '24')
  }

  if (nrow(x) > 1)
  {
    polygon(x = c(x[, 1], rev(x[, 1])), y = c(x[, 3], rev(x[, 4])),
            border = NA, col = adjustcolor(cols, alpha = 0.1))
    lines(x[, 1], x[, 2], lwd = 1.5, col = cols)
  }
}

vars = c('cases_cumulative_perc', 'deaths_cumulative_perc', 
         'excess_deaths_difference_perc',
         'hospitalisation_cumulative_perc', 'tests_total_cumulative_perc',
         'vaccinated_cumulative_perc', 
         'cases_00_19_perc', 'cases_20_49_perc', 'cases_50_69_perc',
         'Abbott_Architect_perc')

labs = c('sqrt(% cases)', 'sqrt(% deaths)', '% exc. deaths diff.', 
         'sqrt(% hosp.)', 'ln(% tested)', '% vaccinated', 
         '% 0-19', '% 20-49', '% 50-69', '% Abbott tests')

ds = vector('list', length = length(vars))

for (i in seq_along(vars))
{
  ds[[i]] = GetGamFuns(m = m_ref_gam, p = vars[i])
}


pdf(paste0('../plots/figures16.pdf'), width = 5, height = 2, 
    pointsize = 6)
par(mar = c(2.5, 2.0, 0.5, 0.8), mgp = c(1.3, 0.3, 0), tcl = 0.25)

split.screen(c(2, 5))

for (i in seq_along(vars))
{
  screen(i)

  PlotGamFuns(x = ds[[i]], p = vars[i], lab = labs[i])
}

close.screen(all.screens = TRUE)
dev.off()



################################################################################
# SI Fig S17: compare predicted seroprevalences and incidences for GLM and GAM.
################################################################################

pdf('../plots/figures17.pdf', width = 1.7, height = 1.7, 
    pointsize = 6)
par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.3, 0.3, 0), tcl = 0.25)

r = cor(predict(m_ref_glm, type = 'response'), predict(m_ref_gam, type = 'response'))

plot(predict(m_ref_glm, type = 'response'), predict(m_ref_gam, type = 'response'), 
     pch = 16, 
     col = adjustcolor('black', alpha = 0.2), cex = 0.7,
     xlab = 'Predicted seroprevalence (GLM)', 
     ylab = 'Predicted seroprevalence (GAM)')
abline(0, 1, lty = '33')
dev.off()


