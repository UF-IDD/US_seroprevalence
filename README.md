# Supplementary code and data for "Accounting for assay performance when estimating the temporal dynamics in SARS-CoV-2 seroprevalence in the U.S."

This repository includes all the R scripts and data necessary to reproduce the
analyses in the manuscript. 

## Data

- [supporting_data.csv](data/supporting_data.csv): The processed data required
to run the model, including the seroprevalence estimates, numbers of reported
cases and deaths, excess deaths, numbers of hospitalisations and tests, and
what assays were used.

- [blood_donors_surveys_processed.csv](data/blood_donors_surveys_processed.csv):
The processed version of the blood donors seroprevalence surveys, necessary to
make the comparisons against our estimates of proportions infected.


## Code

- [functions.R](code/functions.R): Functions necessary to run the analyses, and
model formulae.

- [01_glm_main_reference_model.R](code/seroprevalence_glm_main_model.R): Runs
the main reference model (single GLM model with no waning). This produces
[output/main_glm.rds](output/main_glm.rds).

- [02_glm_waning_models.R](code/seroprevalence_glm_3wanings.R): Runs the waning
models. This script is intended to run on an HPC due to the large number of
models being fit. This produces a large number of files, which in the analyses
are assumed to be combined into a single file, which is here provided as
[output/glm_3-wanings.rds](output/glm_3-wanings.rds).

- [03_get_overall_uncertainty.R](code/seroprevalence_get_overall_uncertainty.R):
Estimates uncertainty in our estimates of proportions infected to account for
uncertainty in the selection of times to seroreversion for the three assays,
and in the lead/lag times between a case being reported and seroconverting.
This produces a large number of files, which in the analyses are assumed to be
combined into a single file, which is here provided as
[output/cis_glm_3-wanings.rds](output/cis_glm_3-wanings.rds).

- [04_gam_main_reference_model.R](code/seroprevalence_gam_main_model.R): Runs
the main reference model (single GAM with no waning). This produces
[output/main_gam.rds](output/main_gam.rds)

- [05_figures.R](code/figures.R): Reproduces all figures in the manuscript
(including supporting figures).



## Output

- [main_glm.rds](output/main_glm.rds): GLM object of the main reference GLM
model, including RMSE and LOO median RMSE. The output is produced by
[01_glm_main_reference_model.R](code/01_glm_main_reference_model.R).

- [glm_3-wanings.rds](output/glm_3-wanings.rds): A table with summary metrics
for all the waning models. The different combinations of times to seroreversion
and lead/lag time are arranged by row, with the corresponding AIC, RMSE, and
LOO median RMSE for each combination. The output is produced by
[02_glm_waning_models.R](code/02_glm_waning_models.R).

- [cis_glm_3-wanings.rds](output/cis_glm_3-wanings.rds): A table with the lower
and upper uncertainty boundaries for both predicted seroprevalences from the
models, and the estimated proportions infected. The uncertainty accounts for
the 95% confidence interval of each individual waning GLM, and the uncertainty
due to the selection of the times to seroreversion of the three assays and the
lead/lag time (uncertainty across GLMs). The output is produced by
[03_get_overall_uncertainty.R](code/03_get_overall_uncertainty.R).

- [main_gam.rds](output/main_gam.rds): GAM object of the main reference GAM
model, including RMSE and LOO median RMSE. The output is produced by
[04_gam_main_reference_model.R](code/04_gam_main_reference_model.R).

- [proportions_infected_estimates.csv](output/proportions_infected_estimates.csv):
Estimates of proportions infected, with uncertainty, from the best waning
model.

