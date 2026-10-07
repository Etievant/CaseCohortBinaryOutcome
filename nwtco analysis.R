## -----------------------------------------------------------------------------
## Description: This script performs parameter and variance estimation under
##              a weighted logistic regression model when using the cohort 
##              subsampling design on the nwtco data set with 
##              - sampling of the subcohort stratified on covariates and 
##              exhaustive sampling of the cases
##              - sampling stratified on covariates and on case status, so that 
##              only a fraction of the cases are included
##
##              Three inference methods are used:
##              - Proposed variance estimation based on the influences computed 
##              from the glm output
##              - Variance estimation from a summary of the svyglm output
##              - Bootstrap variance estimation with a scheme as in Wacholder et
##              al. (Biometrika, 1989) with B = 2,000 bootstrap samples               
##
##              Sampling from the cohort and subsequent parameter and variance
##              estimation is also replicated 50 times to investigate running
##              time with the three methods
## -----------------------------------------------------------------------------
## Required Package: glmtoolbox, survey, survival, tictoc, xtable
## -----------------------------------------------------------------------------
## Required Functions: help.functions.R
## -----------------------------------------------------------------------------

library(glmtoolbox)
library(survey)
library(survival)
library(tictoc)
library(xtable)

source("help.functions.R")

## Covariates and outcome ------------------------------------------------------
data("nwtco")

W     <- nwtco$instit
W.cat <- c(1,2)

X1    <- nwtco$histol
X2    <- nwtco$stage

Y     <- nwtco$rel

# ------------------------------------------------------------------------------
# Design 1: stratified sampling based on institutional histology and 
# exhaustive sampling of the cases ---------------------------------------------

n                   <- nrow(nwtco)
n.cases             <- sum(Y)
strata.n.cases      <- c(sum(W == 1 & Y == 1), sum(W == 2 & Y == 1))
strata.n.noncases   <- c(sum(W == 1 & Y == 0), sum(W == 2 & Y == 0))
strata.n            <- c(sum(W == 1), sum(W == 2))
strata.m            <- floor(strata.n.cases / (1 - strata.n.cases / strata.n))

cohort    <- cbind(id = 1:n, X1, X2, Y, strata = W, n, n.cases, 
                   strata.n = strata.n[W],
                   strata.n.cases = strata.n.cases[W], 
                   strata.n.noncases = strata.n.noncases[W], 
                   strata.m = strata.m[W])
cohort    <- as.data.frame(cohort)

## Subcohort sampling ----------------------------------------------------------
set.seed(12345)
ind.sampled <- NULL
for (w in W.cat) { 
  ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$strata == w], 
                                       size = strata.m[w], 
                                       replace = FALSE)) 
  }
subcohort               <- rep(0, n)
subcohort[ind.sampled]  <- 1
cohort$subcohort        <- subcohort

## Exhaustive sampling of the cases --------------------------------------------
casecohort                      <- subcohort
casecohort[cohort$Y == 1]       <- 1
cohort$casecohort               <- casecohort

## Design weights --------------------------------------------------------------
cohort$weights                  <- 0
cohort$weights[subcohort == 1]  <- cohort$strata.n[subcohort == 1] / cohort$strata.m[subcohort == 1] # 1 / stratum-specific sampling fraction
cohort$weights[cohort$Y == 1]   <- 1

casecohort <- cohort
casecohort <- casecohort[-which((cohort$casecohort == 0)), ]

## Parameter and variance estimates obtained from the full cohort --------------
mod.cohort          <- glm(Y ~ X1 * factor(X2), 
                           family = binomial(link = "logit"), 
                           data = cohort, x = TRUE, y = TRUE)
beta.cohort.glm     <- mod.cohort$coefficients
var.cohort.glm      <- diag(summary(mod.cohort)$cov.unscaled)

## Parameter and variance estimates obtained from the phase-two data and with 
## the glm function ------------------------------------------------------------
mod.casecohort      <- glm(Y ~ X1 * factor(X2), 
                           family = binomial(link = "logit"), 
                           data = casecohort, weights = weights,  x = TRUE, 
                           y = TRUE)
beta.casecohort.glm <- mod.casecohort$coefficients

## Using the output from the glm function, derivation of the influences by hand
## and computation of proposed variance estimates ------------------------------
est.casecohort            <- influences.glm(mod = mod.casecohort)
var.casecohort.infl.glm   <- variance(casecohort, est.casecohort$infl.beta, 
                                        stratified = TRUE) 

## Parameter and variance estimates obtained from the phase-two data and with 
## the survey package ----------------------------------------------------------
des         <- twophase(id = list(~id,~id), 
                        strata = list(NULL, ~interaction(strata, Y)), 
                        subset = ~I(subcohort | Y), data = cohort)

mod.survey  <- svyglm(Y ~ X1 * factor(X2), 
                      family = quasibinomial(link = "logit"), design = des)

beta.casecohort.survey  <- mod.survey$coefficients
var.casecohort.survey   <- diag(summary(mod.survey)$cov.scaled)

## Variance estimates obtained from the phase-two data and with Wacholder's 
## bootstrap scheme and the glm function using ---------------------------------
B <- 2000
strata            <- levels(as.factor(casecohort$strata))
strata.n.noncases <- integer(length(strata))
strata.noncases   <- vector("list", length(strata))
for (i in seq_along(strata)) {
  strata.noncases[[i]] <- which(casecohort$strata == strata[i] & casecohort$Y == 0)
  strata.n.noncases[i] <- sum(casecohort$strata == strata[i] & casecohort$Y == 0)
}
n.cases <- sum(casecohort$Y == 1)
cases   <- which(casecohort$Y == 1)
b.boot  <- NULL
for (b in 1:B) {
  noncase.boot <- NULL
  for (i in seq_along(strata)) {
    noncase.boot <- c(noncase.boot, sample(strata.noncases[[i]], 
                                           size = strata.n.noncases[i], 
                                           replace = TRUE))
  }
  cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
  casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
  
  trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                           family = binomial(link = "logit"), 
                           data = casecohort.boot, weights = weights, 
                           x = TRUE, y = TRUE), silent = TRUE)
  
  while(inherits(trial.boot, "try-error")){
    noncase.boot <- NULL
    for (i in seq_along(strata)) {
      noncase.boot <- c(noncase.boot, sample(strata.noncases[[i]], 
                                             size = strata.n.noncases[i], 
                                             replace = TRUE))
    }
    cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
    casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
    
    trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                             family = binomial(link = "logit"), 
                             data = casecohort.boot, weights = weights, 
                             x = TRUE, y = TRUE), silent = TRUE)
  }
  
  mod.casecohort  <- glm(Y ~ X1 * factor(X2), family = binomial(link = "logit"), 
                         data = casecohort.boot, weights = weights, 
                         x = TRUE, y = TRUE)
  b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
}
var.casecohort.wacholder  <- apply(b.boot, MARGIN = 2, FUN = var)


## Printing the Latex tables for the manuscript --------------------------------
res <- cbind(beta.casecohort.glm, var.casecohort.infl.glm, 
             var.casecohort.wacholder, beta.casecohort.survey, 
             var.casecohort.survey)
res2 <- res
res2[c(1:2,8),c(1,4)]     <- round(res[c(1:2,8),c(1,4)], digits = 1)
res2[c(1:2,3:8),c(2:3,5)] <- format(res[c(1:2,3:8),c(2:3,5)], scientific = TRUE, 
                                    digits = 2)
res2[3:7,c(1,4)]          <- format(res[3:7,c(1,4)], scientific = TRUE, 
                                    digits = 2)
res <- cbind(beta.cohort.glm, var.cohort.glm)
res3 <- res
res3[c(1:2,8),1]  <- round(res[c(1:2,8),1], digits = 1)
res3[c(1:2,8),2]  <- format(res[c(1:2,8),2], scientific = TRUE, digits = 2)
res3[3:7,1:2]     <- format(res[3:7,1:2], scientific = TRUE, digits = 2)
res_exhaustive    <- cbind(res3,res2)
res_exhaustive    <- res_exhaustive[,c(1:3,6,4,5,7)]
print(xtable(res_exhaustive), include.rownames = TRUE)

myfile  <- paste0("DataIllustration.nwtco.logistic.exhaustiveCC.Rdata")
save(res_exhaustive, file = myfile)


# ------------------------------------------------------------------------------
# Design 2: stratified sampling based on case status and institutional 
# histology --------------------------------------------------------------------

W.Y               <- interaction(W, Y) 
n                 <- length(Y)
n.cases           <- sum(Y)
strata.n.cases    <- c(sum(W == 1 & Y == 1), sum(W == 2 & Y == 1))
strata.n.noncases <- c(sum(W == 1 & Y == 0), sum(W == 2 & Y == 0))
strata.n          <- c(strata.n.noncases, strata.n.cases)
strata.m          <- floor(rep(strata.n.cases * 0.8, 2))

cohort    <- cbind(id = 1:n, X1, X2, Y, W, strata = W.Y, n, n.cases, 
                   strata.n = strata.n[2*Y + W], 
                   strata.m = strata.m[W])
cohort    <- as.data.frame(cohort)


## Sampling of the phase-two sample --------------------------------------------
set.seed(12345)
ind.sampled <- NULL
for (w in W.cat) { m <- strata.m[w]
ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 0], 
                                     size = strata.m[w], 
                                     replace = FALSE)) 
ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 1], 
                                     size = strata.m[w], 
                                     replace = FALSE)) 
}
casecohort               <- rep(0, n)
casecohort[ind.sampled]  <- 1
cohort$casecohort        <- casecohort

## Design weights --------------------------------------------------------------
cohort$weights                  <- 0
cohort$weights[casecohort == 1] <- cohort$strata.n[casecohort == 1] / cohort$strata.m[casecohort == 1] # 1 / stratum-specific sampling fraction
casecohort <- cohort
casecohort <- casecohort[-which((cohort$casecohort == 0)), ]

## Parameter and variance estimates obtained from the phase-two data and with 
## the glm function ------------------------------------------------------------
mod.casecohort      <- glm(Y ~ X1 * factor(X2), 
                           family = binomial(link = "logit"), 
                           data = casecohort, weights = weights, 
                           x = TRUE, y = TRUE)
beta.casecohort.glm <- mod.casecohort$coefficients

## Using the output from the glm function, derivation of the influences by hand
## and computation of proposed variance estimates ------------------------------
est.casecohort              <- influences.glm(mod = mod.casecohort)
var.casecohort.infl.glm     <- variance(casecohort, est.casecohort$infl.beta, 
                                        stratified = TRUE) 

## Estimation and inference using the survey package ---------------------------
des         <- twophase(id = list(~id,~id), strata = list(NULL, ~strata), 
                        subset = ~as.logical(casecohort), data = cohort)
mod.survey  <- svyglm(Y ~ X1 * factor(X2), 
                      family = quasibinomial(link = "logit"), design = des)

beta.casecohort.survey  <- mod.survey$coefficients
var.casecohort.survey   <- diag(summary(mod.survey)$cov.scaled)

## Variance estimates obtained from the phase-two data and with Wacholder's 
## bootstrap scheme and the glm function using ---------------------------------
B   <- 2000
strata          <- levels(as.factor(casecohort$strata))
strata.indiv    <- vector("list", length(strata))
strata.n.indiv  <- integer(length(strata))
for (i in seq_along(strata)) {
  strata.indiv[[i]] <- which(casecohort$strata == strata[i])
  strata.n.indiv[i] <- sum(casecohort$strata == strata[i])
}
b.boot <- NULL
for (b in 1:B) {
  casecohort.boot <- NULL
  for (i in seq_along(strata)) {
    casecohort.boot <- rbind(casecohort.boot, 
                             casecohort[sample(strata.indiv[[i]], 
                                               size = strata.n.indiv[i], 
                                               replace = TRUE),])
  }
  trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                           family = binomial(link = "logit"), 
                           data = casecohort.boot, weights = weights, x = TRUE, 
                           y = TRUE), silent = TRUE)
  while(inherits(trial.boot, "try-error")){
    casecohort.boot <- NULL
    for (i in seq_along(strata)) {
      casecohort.boot <- rbind(casecohort.boot, 
                               casecohort[sample(strata.indiv[[i]], 
                                                 size = strata.n.indiv[i],
                                                 replace = TRUE),])
    }
    trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                             family = binomial(link = "logit"), 
                             data = casecohort.boot, weights = weights, x = TRUE, 
                             y = TRUE), silent = TRUE)
  }
  mod.casecohort  <- glm(Y ~ X1 * factor(X2), 
                         family = binomial(link = "logit"), 
                         data = casecohort.boot, weights = weights, x = TRUE, 
                         y = TRUE)
  b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
}
var.casecohort.wacholder  <- apply(b.boot, MARGIN = 2, FUN = var)

## Printing the Latex tables for the manuscript --------------------------------
res <- cbind(beta.casecohort.glm, var.casecohort.infl.glm, 
             var.casecohort.wacholder, beta.casecohort.survey, 
             var.casecohort.survey)
res2 <- res
res2[c(1:2,8),c(1,4)]     <- round(res[c(1:2,8),c(1,4)], digits = 1)
res2[c(1:2,3:8),c(2:3,5)] <- format(res[c(1:2,3:8),c(2:3,5)], scientific = TRUE, 
                                    digits = 2)
res2[3:7,c(1,4)]          <- format(res[3:7,c(1,4)], scientific = TRUE, 
                                    digits = 2)
res <- cbind(beta.cohort.glm, var.cohort.glm)
res3 <- res
res3[c(1:2,8),1]  <- round(res[c(1:2,8),1], digits = 1)
res3[c(1:2,8),2]  <- format(res[c(1:2,8),2], scientific = TRUE, digits = 2)
res3[3:7,1:2]     <- format(res[3:7,1:2], scientific = TRUE, digits = 2)

res_strat <- cbind(res3,res2)
res_strat <- res_strat[,c(3,6,4,5,7)]
myfile    <- paste0("DataIllustration.nwtco.logistic.stratifiedCC.Rdata")
save(res_strat, file = myfile)
print(xtable(cbind(res_exhaustive, res_strat)), include.rownames = TRUE)


# ------------------------------------------------------------------------------
# Comparing running time with the different variance estimation commands -------

library(tictoc)
set.seed(12345)
Nreplic <- 50
res     <- NULL
for (nrep in 1:Nreplic) {
  ## Derivating the influences by hand and computing the proposed variance 
  ## estimate ------------------------------------------------------------------
  tic.clear()
  tic()
  ind.sampled <- NULL
  for (w in W.cat) { m <- strata.m[w]
  ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 0], 
                                       size = strata.m[w], 
                                       replace = FALSE)) 
  ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 1], 
                                       size = strata.m[w], 
                                       replace = FALSE)) 
  }
  casecohort               <- rep(0, n)
  casecohort[ind.sampled]  <- 1
  cohort$casecohort        <- casecohort
  cohort$weights                  <- 0
  cohort$weights[casecohort == 1] <- cohort$strata.n[casecohort == 1] / cohort$strata.m[casecohort == 1] # 1 / stratum-specific sampling fraction
  casecohort <- cohort
  casecohort <- casecohort[-which((cohort$casecohort == 0)), ]
  
  mod.casecohort      <- glm(Y ~ X1 * factor(X2), 
                             family = binomial(link = "logit"), 
                             data = casecohort, weights = weights, 
                             x = TRUE, y = TRUE)
  beta.casecohort.glm         <- mod.casecohort$coefficients
  est.casecohort              <- influences.glm(mod = mod.casecohort)
  var.casecohort.infl.glm     <- variance(casecohort, est.casecohort$infl.beta, 
                                          stratified = TRUE) 
  run         <- toc(quiet = TRUE)
  time.casecohort <- run$toc - run$tic
  
  ## Using the survey package --------------------------------------------------
  tic.clear()
  tic()
  des                     <- twophase(id = list(~id,~id), 
                                      strata = list(NULL, ~strata), 
                                      subset = ~as.logical(casecohort), 
                                      data = cohort)
  mod.survey              <- svyglm(Y ~ X1 * factor(X2), 
                                    family = quasibinomial(link = "logit"), 
                                    design = des)
  beta.casecohort.survey  <- mod.survey$coefficients
  var.casecohort.survey   <- diag(summary(mod.survey)$cov.scaled)
  run         <- toc(quiet = TRUE)
  time.survey <- run$toc - run$tic
  
  ## Using Wacholder's bootstrap scheme ----------------------------------------
  strata <- levels(as.factor(casecohort$strata))
  strata.indiv <- vector("list", length(strata))
  strata.n.indiv <- integer(length(strata))
  for (i in seq_along(strata)) {
    strata.indiv[[i]] <- which(casecohort$strata == strata[i])
    strata.n.indiv[i] <- sum(casecohort$strata == strata[i])
  }
  B <- 2000
  tic.clear()
  tic()
  b.boot <- NULL
  for (b in 1:B) {
    casecohort.boot <- NULL
    for (i in seq_along(strata)) {
      casecohort.boot <- rbind(casecohort.boot, 
                               casecohort[sample(strata.indiv[[i]], 
                                                 size = strata.n.indiv[i], 
                                                 replace = TRUE),])
    }
    trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                             family = binomial(link = "logit"), 
                             data = casecohort.boot, weights = weights, 
                             x = TRUE, y = TRUE), silent = TRUE)
    while(inherits(trial.boot, "try-error")){
      casecohort.boot <- NULL
      for (i in seq_along(strata)) {
        casecohort.boot <- rbind(casecohort.boot, 
                                 casecohort[sample(strata.indiv[[i]], 
                                                   size = strata.n.indiv[i], 
                                                   replace = TRUE),])
      }
      trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                               family = binomial(link = "logit"), 
                               data = casecohort.boot, weights = weights, 
                               x = TRUE, y = TRUE), silent = TRUE)
    }
    mod.casecohort  <- glm(Y ~ X1 * factor(X2), 
                           family = binomial(link = "logit"), 
                           data = casecohort.boot, weights = weights, x = TRUE, 
                           y = TRUE)
    b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
  }
  var.casecohort.wacholder  <- apply(b.boot, MARGIN = 2, FUN = var)
  run         <- toc(quiet = TRUE)
  time.wacholder.2000  <- run$toc - run$tic
  
  B <- 500
  tic.clear()
  tic()
  b.boot <- NULL
  for (b in 1:B) {
    casecohort.boot <- NULL
    for (i in seq_along(strata)) {
      casecohort.boot <- rbind(casecohort.boot, 
                               casecohort[sample(strata.indiv[[i]], 
                                                 size = strata.n.indiv[i], 
                                                 replace = TRUE),])
    }
    trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                             family = binomial(link = "logit"), 
                             data = casecohort.boot, weights = weights, 
                             x = TRUE, y = TRUE), silent = TRUE)
    while(inherits(trial.boot, "try-error")){
      casecohort.boot <- NULL
      for (i in seq_along(strata)) {
        casecohort.boot <- rbind(casecohort.boot, 
                                 casecohort[sample(strata.indiv[[i]], 
                                                   size = strata.n.indiv[i], 
                                                   replace = TRUE),])
      }
      trial.boot    <- try(glm(Y ~ X1 * factor(X2), 
                               family = binomial(link = "logit"), 
                               data = casecohort.boot, weights = weights, 
                               x = TRUE, y = TRUE), silent = TRUE)
    }
    mod.casecohort  <- glm(Y ~ X1 * factor(X2), 
                           family = binomial(link = "logit"), 
                           data = casecohort.boot, weights = weights, x = TRUE, 
                           y = TRUE)
    b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
  }
  var.casecohort.wacholder  <- apply(b.boot, MARGIN = 2, FUN = var)
  run                       <- toc(quiet = TRUE)
  time.wacholder.500        <- run$toc - run$tic
  
  recap   <- c(time.casecohort = time.casecohort, time.survey = time.survey, 
               time.wacholder.2000 = time.wacholder.2000,
               time.wacholder.500 = time.wacholder.500)
  res     <- rbind(res, recap)
}
myfile  <- paste0("DataIllustration.nwtco.logistic.stratifiedCC.runningtime.Rdata")
save(res, file = myfile)
colMeans(res)
