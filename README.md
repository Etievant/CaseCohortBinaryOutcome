Replication of the simulation study and data analysis in "Variance estimation for generalized linear regression with binary outcome when using case-cohort data", by Etiévant (2026). The method proposed in the article and the other methods it is compared to are implemented in the file `help.functions.R`.

### Required packages 

```
dplyr, glmtoolbox, ggh4x, ggplot2, openxlsx, parallel, pracma, survey, survival, tictoc, xtable.

```

### Scripts

* Scripts `simulations_logistic_exhaustive.R` and `simulations_loglinear_exhaustive.R` replicate the simulation proposed by Etiévant (2026) in Section 6.1.1.

* Scripts `analysis_logistic_exhaustive.R` and `analysis_loglinear_exhaustive.R` summarize the simulation results as shown in Section 6.2.1 and Supplementary Section E.1.1.

* Scripts `simulations_logistic_stratified.R` and `simulations_loglinear_stratified.R` replicate the simulation proposed by Etiévant (2026) in Section 6.1.2.

* Scripts `analysis_logistic_stratified.R` and `analysis_loglinear_stratified.R` summarize the simulation results as shown in Section 6.2.2 and Supplementary Section E.1.2.

* Script `nwtco_analysis.R` replicates the data illustration in Etiévant and Gail (2026) in Section 7.

Each simulation script relies on functions provided in `help.functions.R`.


### Instructions to run each script

* Save the chosen simulation script(s) and file `help.functions.R` in the same directory.

* Open and run the whole script(s).

* The results of the simulations are saved in csv tables and Rdata and xlsx files. For example, when running script `simulations_loglinear_exhaustive.R` and then script `analysis_loglinear_exhaustive.R`, file `details.beta0_loglinear_exhaustive.csv` will contain the simulation results given in Section 6.2.1 and in Supplementary Section E.1.1. Script `analysis_loglinear_exhaustive.R` will also print TeX tables in the same format as displayed in the article.


### Functions in `help.functions.R`

* **influences.glm** - computes the influences on the beta parameter estimates in a log-linear or logistic regression model, when the parameter has been estimated with the glm command from the stats R package. This function can be used under the design with unstratified sampling of the subcohort or stratified sampling based on some covariate and exhaustive sampling of cases, or under the design with stratified sampling based on case status and possibly some covariate.

* **influences.glmgee** - computes the influences on the beta parameter estimates in a log-linear or logistic regression model, when the parameter has been estimated with the glmgee command from the glmtoolbox R package. This function can be used under the design with unstratified sampling of the subcohort or stratified sampling based on some covariate and exhaustive sampling of cases, or under the design with stratified sampling based on case status and possibly some covariate.

* **variance** - returns the influence-based variance estimate that follows the complete variance decomposition. This function can be used with phase-two data obtained from the design with unstratified sampling of the subcohort or stratified sampling based on some covariate and exhaustive sampling of cases, or under the design with stratified sampling based on case status and possibly some covariate.

* **CI** - computes the 100*(1 - alpha) percent confidence interval using the quantiles from the Normal distribution.

Other functions are provided in `help.functions.R`. Notably, analyses of cohort and case-cohort data with the different methods and readily available software considered in the article have been implemented within functions `cohort.analysis` and `casecohort.analysis`. Variance estimations from the two bootstrap schemes considered in the article have been implemented within functions `casecohort.boot.wacholder` and `casecohort.boot.huang`.


### Reference

Etiévant L (2026) Variance estimation for generalized linear regression with binary outcome when using case-cohort data. preprint

