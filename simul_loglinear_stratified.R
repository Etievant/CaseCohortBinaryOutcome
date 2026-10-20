## -----------------------------------------------------------------------------
## Description: This script performs parameter and variance estimation under
##              a weighted log-linear regression model when using the cohort 
##              subsampling design with 
##              - sampling stratified on covariates and on case status, so that 
##              only a fraction of the cases are included
##              - stratum-specific design weights
##              
##              Several inference methods are used:
##              - Proposed variance estimation based on the influences computed 
##              from the glm and glmgee outputs
##              - Variance estimation as returned by the vcov.glmgee command on 
##              the glmgee output
##              - Variance estimation from a summary of the svyglm output
##              - Bootstrap variance estimation with a scheme as in Wacholder et
##              al. (Biometrika, 1989) with B = 2,000 bootstrap samples
##              - Variance estimation as returned by the glm command
##              - Robust variance estimation based on the influences computed 
##              from the glm and glmgee outputs
##              - Bootstrap variance estimation with a scheme as in Huang et al.
##              (Biometrika, 2014) with B = 2,000 bootstrap samples
##              - Proposed variance estimation based on the dfbetas as returned 
##              by the dfbeta.glmgee command on the glmgee output
##              - Robust variance estimation based on the dfbetas as returned by
##              the the dfbeta.glmgee command on the glmgee output
##
##              95% CIs are computed assuming normality, for the bootstrap 
##              methods as well, in addition to Percentile CIs
## -----------------------------------------------------------------------------
## Required Package: glmtoolbox, parallel, survey
## -----------------------------------------------------------------------------
## Required Functions: help.functions.R
## -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Read functions and packages --------------------------------------------------

library(glmtoolbox)
library(parallel)
library(survey)

source("help.functions.R")

# ------------------------------------------------------------------------------
# Different scenarios (with different combinations of parameters) to be 
# investigated over simulations ------------------------------------------------

PY    <- rbind(0.2, 0.3)  # P(Y = 1)
N     <- 5000         # cohort size
sP    <- c(0.5, 0.8)  # sampling probability for the cases in the phase-two data
part  <- 5            # to break the number of replications

combinations        <- expand.grid(part = 1:part, n = N, PY = PY, sP = sP)
scenarios           <- cbind(combinations$PY, combinations$n, combinations$sP,
                             combinations$part)
colnames(scenarios) <- c("PY", "n", "sP", "part")

# ------------------------------------------------------------------------------
# Fixed parameters -------------------------------------------------------------

PX1     <- c(0.1, 0.05, 0.5, 0.35) # P(X1 = j)
X1.cat  <- 1:4

PW      <- c(sum(PX1[1:2]), 1 - sum(PX1[1:2])) # P(W = w)
W.cat   <- 0:1 # P(W = 0 | X1 = 1) = P(W = 0 | X1 = 2) = 1 and P(W = 1 | X1 = 3) = P(W = 1 | X1 = 4) = 1

PX2_X1  <- matrix(c(0.7, 0.3, 0.5, 0.5, 0.15, 0.85, 0.4, 0.6), nrow = 4, 
                  byrow = TRUE) # P(X2 = 0 | X1 = j) and P(X2 = 1 | X1 = j) on the j-th row

alpha1  <- 0.25 # parameters used for the X3 mean 
alpha2  <- -0.35

beta1  <- 0.45 
beta2  <- -0.35
beta3  <- -0.25 
names.beta <- paste0("beta", c(0,1,2,3))

EY <- 0 # will be used to compute the intercept, beta0, to control E(Y) = P(Y = 1)
for(x1 in X1.cat) {
  for(x2 in 0:1) {
    fX3_X1X2 <- function (x3) {
      return(exp(beta3 * x3) * 1 / sqrt(2 * pi) * exp(-1/2 * (x3 - (alpha1 * x1 + alpha2 * x2))^2))
    }
    EY <- EY + integrate(Vectorize(fX3_X1X2), -10, 10)$value * 
      exp(beta2 * x2) * PX2_X1[x1, x2 + 1] * exp(beta1 * x1) * PX1[x1]
  }
}

EY_W0 <- 0 # E(Y | W = 0), will be used to choose the stratum-specific sampling fractions
for(x1 in 1:2) {
  for(x2 in 0:1) {
    fX3_X1X2 <- function (x3) {
      return(exp(beta3 * x3) * 1 / sqrt(2 * pi) * exp(-1/2 * (x3 - (alpha1 * x1 + alpha2 * x2))^2))
    }
    EY_W0 <- EY_W0 + integrate(Vectorize(fX3_X1X2), -10, 10)$value * 
      exp(beta2 * x2) * PX2_X1[x1, x2 + 1] * exp(beta1 * x1) * PX1[x1] / PW[1]
  }
}

EY_W1 <- 0 # E(Y | W = 1)
for(x1 in 3:4) {
  for(x2 in 0:1) {
    fX3_X1X2 <- function (x3) {
      return(exp(beta3 * x3) * 1 / sqrt(2 * pi) * exp(-1/2 * (x3 - (alpha1 * x1 + alpha2 * x2))^2))
    }
    EY_W1 <- EY_W1 + integrate(Vectorize(fX3_X1X2), -10, 10)$value * 
      exp(beta2 * x2) * PX2_X1[x1, x2 + 1] * exp(beta1 * x1) * PX1[x1] / PW[2]
  }
}

EY_W = c(EY_W0, EY_W1) # sum(EY_W * PW) = EY

B       <- 2000 # number of bootstrap samples
Nreplic <- 1000 # 5000 replications in total
set.seed(12345)
seed <- round(abs(rnorm(5) * 10^5))

Onerun <- function (p) {
  pY      <- scenarios[p, 1]
  n       <- scenarios[p, 2]
  sp      <- scenarios[p, 3]
  part    <- scenarios[p, 4]
  
  set.seed(seed[part])
  
  # Computation of the intercept, so that E(Y) = pY
  beta0   <- log(pY / EY)
  
  # Computation of the stratum-specific sampling fractions
  strata.m    <- floor(sp * n * PW * exp(beta0) * EY_W) # number of cases to be sampled in each stratum (based on W)
  # the same number of non-cases will be sampled
  
  res     <- NULL
  
  for (nrep in 1:Nreplic) {

    # --------------------------------------------------------------------------
    # Generation of the cohort -------------------------------------------------
    
    # Generation of the covariates ---------------------------------------------
    
    X1      <- sample(1:4, size = n, prob = PX1, replace = TRUE)
    
    X2      <- rep(NA, n)
    for (i in 1:length(X1.cat)) {
      X2[which(X1 == X1.cat[i])] <- rbinom(sum(X1 == X1.cat[i]), size = 1, 
                                           prob = PX2_X1[i,2])
    }
    
    mu3   <- alpha1 * X1 + alpha2 * X2
    X3    <- rnorm(n, mean = mu3)
    
    q     <- 3 # number of covariates
    
    # Generation of the outcome ------------------------------------------------
    
    p   <- exp(beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3)
    Y   <- rbinom(n, size = 1, prob = p) 
    
    # Generation of the strata -------------------------------------------------
    
    W   <- 0 * (X1 < 3) + 1 * (X1 >= 3) # correlated with X1
    W.Y <- interaction(W, Y) # sampling will be also stratified based on case status
    
    # Number of cases and subcohort size ---------------------------------------
    
    n.cases           <- sum(Y)
    strata.n.cases    <- c(sum(W == 0 & Y == 1), sum(W == 1 & Y == 1))
    strata.n.noncases <- c(sum(W == 0 & Y == 0), sum(W == 1 & Y == 0))
    strata.n          <- c(strata.n.noncases, strata.n.cases)
    
    while((sum(strata.n.cases < strata.m) > 0)|(max(p) > 1)){
      
      # Generate a new cohort 
      X1      <- sample(1:4, size = n, prob = PX1, replace = TRUE)
      X2      <- rep(NA, n)
      for (i in 1:length(X1.cat)) {
        X2[which(X1 == X1.cat[i])] <- rbinom(sum(X1 == X1.cat[i]), size = 1, 
                                             prob = PX2_X1[i,2])
      }
      mu3   <- alpha1 * X1 + alpha2 * X2
      X3    <- rnorm(n, mean = mu3)
      p   <- exp(beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3)
      Y   <- rbinom(n, size = 1, prob = p) 
      W   <- 0 * (X1 < 3) + 1 * (X1 >= 3) # correlated with X1
      W.Y <- interaction(W, Y) # sampling will be also stratified based on case status
      
      n.cases           <- sum(Y)
      strata.n.cases    <- c(sum(W == 0 & Y == 1), sum(W == 1 & Y == 1))
      strata.n.noncases <- c(sum(W == 0 & Y == 0), sum(W == 1 & Y == 0))
      strata.n          <- c(strata.n.noncases, strata.n.cases)
    }
    
    cohort    <- cbind(id = 1:n, X1, X2, X3, Y, W, strata = W.Y, n, n.cases, 
                       strata.n.cases = strata.n.cases[W + 1],
                       strata.n.noncases = strata.n.noncases[W + 1],
                       strata.n = strata.n[2*Y + W + 1], 
                       strata.m = strata.m[W + 1])
    cohort    <- as.data.frame(cohort)
    
    # --------------------------------------------------------------------------
    # Sampling of the phase-two sample -----------------------------------------
    
    ind.sampled <- NULL
    for (w in W.cat) { 
      m <- strata.m[w + 1]
      ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 0], 
                                           size = strata.m[w + 1], 
                                           replace = FALSE)) 
      ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 1], 
                                           size = strata.m[w + 1], 
                                           replace = FALSE)) 
    }
    casecohort               <- rep(0, n)
    casecohort[ind.sampled]  <- 1
    cohort$casecohort        <- casecohort
    
    # Design weights -----------------------------------------------------------
    
    cohort$weights                  <- 0
    cohort$weights[casecohort == 1] <- cohort$strata.n[casecohort == 1] / cohort$strata.m[casecohort == 1] # 1 / stratum-specific sampling fraction
    
    casecohort <- cohort
    casecohort <- casecohort[-which((cohort$casecohort == 0)), ]
    
    # --------------------------------------------------------------------------
    # Checking the convergence of the glm and glmgee functions -----------------
    
    trial.casecohort  <- try(casecohort.analysis(casecohort = casecohort, 
                                                 link = "log", 
                                                 stratified = TRUE), 
                             silent = TRUE)
    trial.cohort      <- try(cohort.analysis(cohort = cohort, link = "log"), 
                             silent = TRUE)
    
    while(inherits(trial.casecohort, "try-error") | inherits(trial.cohort, "try-error")){

      # Generate a new cohort and phase-two sample -----------------------------
      X1      <- sample(1:4, size = n, prob = PX1, replace = TRUE)
      X2      <- rep(NA, n)
      for (i in 1:length(X1.cat)) {
        X2[which(X1 == X1.cat[i])] <- rbinom(sum(X1 == X1.cat[i]), size = 1, 
                                             prob = PX2_X1[i,2])
      }
      mu3   <- alpha1 * X1 + alpha2 * X2
      X3    <- rnorm(n, mean = mu3)
      p   <- exp(beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3)
      Y   <- rbinom(n, size = 1, prob = p) 
      W   <- 0 * (X1 < 3) + 1 * (X1 >= 3) # correlated with X1
      W.Y <- interaction(W, Y) # sampling will be also stratified based on case status

      n.cases           <- sum(Y)
      strata.n.cases    <- c(sum(W == 0 & Y == 1), sum(W == 1 & Y == 1))
      strata.n.noncases <- c(sum(W == 0 & Y == 0), sum(W == 1 & Y == 0))
      strata.n          <- c(strata.n.noncases, strata.n.cases)
      
      while((sum(strata.n.cases < strata.m) > 0)|(max(p) > 1)){
        X1      <- sample(1:4, size = n, prob = PX1, replace = TRUE)
        X2      <- rep(NA, n)
        for (i in 1:length(X1.cat)) {
          X2[which(X1 == X1.cat[i])] <- rbinom(sum(X1 == X1.cat[i]), size = 1, 
                                               prob = PX2_X1[i,2])
        }
        mu3   <- alpha1 * X1 + alpha2 * X2
        X3    <- rnorm(n, mean = mu3)
        p   <- exp(beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3)
        Y   <- rbinom(n, size = 1, prob = p) 
        W   <- 0 * (X1 < 3) + 1 * (X1 >= 3) # correlated with X1
        W.Y <- interaction(W, Y) # sampling will be also stratified based on case status
        
        n.cases           <- sum(Y)
        strata.n.cases    <- c(sum(W == 0 & Y == 1), sum(W == 1 & Y == 1))
        strata.n.noncases <- c(sum(W == 0 & Y == 0), sum(W == 1 & Y == 0))
        strata.n          <- c(strata.n.noncases, strata.n.cases)
      }
      
      cohort    <- cbind(id = 1:n, X1, X2, X3, Y, W, strata = W.Y, n, n.cases, 
                         strata.n.cases = strata.n.cases[W + 1],
                         strata.n.noncases = strata.n.noncases[W + 1],
                         strata.n = strata.n[2*Y + W + 1], 
                         strata.m = strata.m[W + 1])
      cohort    <- as.data.frame(cohort)
  
      ind.sampled <- NULL
      for (w in W.cat) { 
        m <- strata.m[w + 1]
        ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 0], 
                                             size = strata.m[w + 1], 
                                             replace = FALSE)) 
        ind.sampled <- c(ind.sampled, sample(cohort$id[cohort$W == w & cohort$Y == 1], 
                                             size = strata.m[w + 1], 
                                             replace = FALSE)) 
      }
      casecohort               <- rep(0, n)
      casecohort[ind.sampled]  <- 1
      cohort$casecohort        <- casecohort
      cohort$weights                  <- 0
      cohort$weights[casecohort == 1] <- cohort$strata.n[casecohort == 1] / cohort$strata.m[casecohort == 1] # 1 / stratum-specific sampling fraction
      casecohort <- cohort
      casecohort <- casecohort[-which((cohort$casecohort == 0)), ]
      
      trial.casecohort  <- try(casecohort.analysis(casecohort = casecohort, 
                                                   link = "log", 
                                                   stratified = TRUE), 
                               silent = TRUE)
      trial.cohort      <- try(cohort.analysis(cohort = cohort, link = "log"), 
                               silent = TRUE)
    }
    
    # --------------------------------------------------------------------------
    # Estimation using the entire cohort ---------------------------------------
    
    esti.cohort         <- cohort.analysis(cohort, link = "log")
    
    res.cohort.glm.var              <- esti.cohort$res.cohort.glm.var
    res.cohort.glm.robvar.infl      <- esti.cohort$res.cohort.glm.robvar.infl
    res.cohort.glmgee.var           <- esti.cohort$res.cohort.glmgee.var
    res.cohort.glmgee.robvar.infl   <- esti.cohort$res.cohort.glmgee.robvar.infl
    res.cohort.glmgee.robvar.dfbeta <- esti.cohort$res.cohort.glmgee.robvar.dfbeta
    
    # --------------------------------------------------------------------------
    # Estimation using the phase-two sample with design weights ----------------
    
    esti.casecohort <- casecohort.analysis(casecohort, stratified = TRUE,
                                           link = "log")
    
    res.casecohort.glm.var              <- esti.casecohort$res.casecohort.glm.var
    res.casecohort.glm.var.infl         <- esti.casecohort$res.casecohort.glm.var.infl
    res.casecohort.glm.robvar.infl      <- esti.casecohort$res.casecohort.glm.robvar.infl
    res.casecohort.glmgee.var           <- esti.casecohort$res.casecohort.glmgee.var
    res.casecohort.glmgee.var.infl      <- esti.casecohort$res.casecohort.glmgee.var.infl
    res.casecohort.glmgee.robvar.infl   <- esti.casecohort$res.casecohort.glmgee.robvar.infl
    res.casecohort.glmgee.var.dfbeta    <- esti.casecohort$res.casecohort.glmgee.var.dfbeta
    res.casecohort.glmgee.robvar.dfbeta <- esti.casecohort$res.casecohort.glmgee.robvar.dfbeta
    
    # estimation and inference using the survey package
    des <- twophase(id = list(~id,~id), 
                    strata = list(NULL, ~strata), 
                    subset = ~as.logical(casecohort), data = cohort)
    
    mod.survey  <- svyglm(Y ~ X1 + X2 + X3, family = quasipoisson(link = "log"), 
                          design = des)
    beta.survey <- mod.survey$coefficients
    v.survey    <- diag(summary(mod.survey)$cov.scaled)
    res.casecohort.survey  <- cbind("casecohort survey", t(beta.survey), 
                                    t(v.survey), t(CI(beta.survey, v.survey)))
    
    # inference using the bootstrap schemes
    # Huang's bootstrap scheme
    esti.casecohort.boot.huang  <- casecohort.boot.huang(casecohort, 
                                                         link = "log", B = B)
    
    res.casecohort.boot.huang.wald        <- esti.casecohort.boot.huang$res.casecohort.boot.wald
    res.casecohort.boot.huang.percentile  <- esti.casecohort.boot.huang$res.casecohort.boot.percentile
    
    # Wacholder's bootstrap scheme
    esti.casecohort.boot.wacholder  <- casecohort.boot.wacholder(casecohort, B = B, 
                                                                 link = "log", 
                                                                 stratified = TRUE, 
                                                                 exhaustive = FALSE)
    
    res.casecohort.boot.wacholder.wald        <- esti.casecohort.boot.wacholder$res.casecohort.boot.wald
    res.casecohort.boot.wacholder.percentile  <- esti.casecohort.boot.wacholder$res.casecohort.boot.percentile
    
    recap <- rbind(res.cohort.glm.var, res.cohort.glm.robvar.infl, 
                   res.cohort.glmgee.var, res.cohort.glmgee.robvar.infl, 
                   res.cohort.glmgee.robvar.dfbeta, 
                   res.casecohort.glm.var, res.casecohort.glm.var.infl, 
                   res.casecohort.glm.robvar.infl, res.casecohort.glmgee.var, 
                   res.casecohort.glmgee.var.infl, 
                   res.casecohort.glmgee.robvar.infl, 
                   res.casecohort.glmgee.var.dfbeta, 
                   res.casecohort.glmgee.robvar.dfbeta, res.casecohort.survey,
                   res.casecohort.boot.huang.wald, res.casecohort.boot.huang.percentile,
                   res.casecohort.boot.wacholder.wald, res.casecohort.boot.wacholder.percentile)
    colnames(recap) <- c("Method", paste0(names.beta, ".hat"),
                         paste0("var.", names.beta, ".hat"),
                         paste0("CI.left.", names.beta, ".hat"),
                         paste0("CI.right.", names.beta, ".hat"))
    
    res <- rbind(res, cbind(recap, n, sp, strata.n_W0 = strata.n[1], 
                            strata.n_W1 = strata.n[2], n.cases, 
                            strata.n.cases_W0 = strata.n.cases[1], 
                            strata.n.cases_W1 = strata.n.cases[2],
                            strata.n.noncases_W0 = strata.n.noncases[1], 
                            strata.n.noncases_W1 = strata.n.noncases[2],
                            strata.m_W0 = strata.m[1],
                            strata.m_W1 = strata.m[2], 
                            n.casecohort = nrow(casecohort), 
                            pY, beta0, beta1, beta2, beta3))
    row.names(res) <- NULL
  }
  myfile  <- paste0("Results_stratified_loglinear-n", n, "-pY", pY, "-Sp", sp, 
                    "-B", B, "-Nreplic", 5000, "-Part", part, ".Rdata")
  save(res, file = myfile)
}

# ------------------------------------------------------------------------------
# Running the simulations in parallel for all the scenarios --------------------

result <- mclapply(1:nrow(scenarios), Onerun, mc.cores = 20)

# ------------------------------------------------------------------------------
# Combining the simulation results ---------------------------------------------

RES <- NULL
for (p in 1:nrow(scenarios)) {
  pY      <- scenarios[p, 1]
  n       <- scenarios[p, 2]
  sp      <- scenarios[p, 3]
  part    <- scenarios[p, 4]
  load(paste0("Results_stratified_loglinear-n", n, "-pY", pY, "-Sp", sp,  "-B", 
              B, "-Nreplic", 5000, "-Part", part, ".Rdata"))
  RES <- rbind(RES, res)
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(2:ncol(RECAP))])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$Method    <- as.factor(RECAP$Method)
RECAP$n         <- as.factor(RECAP$n)
RECAP$sp        <- as.factor(RECAP$sp)
RECAP$pY        <- as.factor(RECAP$pY)
myfile  <- paste0("Results_stratified_loglinear-", B, "_boot-", 5000, 
                  "_Nreplic.Rdata")
save(RECAP, file = myfile)

