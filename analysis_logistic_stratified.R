## -----------------------------------------------------------------------------
## Description: This script analyses the simulation results where we performed 
##              parameter and variance estimation under a logistic regression
##              model when using a case-cohort with 
##              - sampling stratified on covariates and on case status, so that 
##              only a fraction of the cases are included
##              - stratum-specific design weights
##
##              Several inference methods were used:
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
##              95% CIs were computed assuming normality, for the bootstrap 
##              methods as well, in addition to Percentile CIs
## -----------------------------------------------------------------------------
## Required Package: openxlsx, xtable
## -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Read functions and packages --------------------------------------------------

library(openxlsx)
library(xtable)

# Different scenarios (with different combinations of parameters) investigated
# over simulations -------------------------------------------------------------

PY    <- rbind(0.2, 0.3, 0.4)  # P(Y = 1)
N     <- 5000           # cohort size
sP    <- c(0.5, 0.8)    # sampling probability for the cases in the phase-two data
part  <- 5              # to break the number of replications

combinations        <- expand.grid(part = 1:part, n = N, PY = PY, sP = sP)
scenarios           <- cbind(combinations$PY, combinations$n, combinations$sP,
                             combinations$part)
colnames(scenarios) <- c("PY", "n", "sP", "part")

# ------------------------------------------------------------------------------
# Loading the simulation results -----------------------------------------------

B       <- 2000
Nreplic <- 5000

load(paste0("Results_stratified_logistic-", B, "_boot-", Nreplic, "_Nreplic.Rdata"))

# ------------------------------------------------------------------------------
# Details of the results  ------------------------------------------------------

coverage.left   <- round(0.95 - 1.96 * sqrt(0.95 * 0.05 / Nreplic), digits = 4)
coverage.right  <- round(0.95 + 1.96 * sqrt(0.95 * 0.05 / Nreplic), digits = 4)

scenarios <- scenarios[,-4]
scenarios <- scenarios[which(duplicated.matrix(scenarios)==FALSE),]

list.methods <- c("casecohort glm",
                  "casecohort glm var infl",
                  "casecohort glm robvar infl",
                  "casecohort glmgee",
                  "casecohort glmgee var infl",
                  "casecohort glmgee robvar infl",
                  "casecohort glmgee var dfbeta",
                  "casecohort glmgee robvar dfbeta",
                  "casecohort survey",
                  "casecohort Huang boot wald",
                  "casecohort Huang boot percentile",
                  "casecohort Wacholder boot wald",
                  "casecohort Wacholder boot percentile",
                  "cohort glm",
                  "cohort glm robvar infl",
                  "cohort glmgee",
                  "cohort glmgee robvar infl",
                  "cohort glmgee robvar dfbeta")

colname.methods <- c("CC glm",
                     "CC glm var infl",
                     "CC glm robvar infl",
                     "CC glmgee",
                     "CC glmgee var infl",
                     "CC glmgee robvar infl",
                     "CC glmgee var dfbeta",
                     "CC glmgee robvar dfbeta",
                     "CC survey",
                     "CC Huang boot wald",
                     "CC Huang boot percentile",
                     "CC Wacholder boot wald",
                     "CC Wacholder boot percentile",
                     "cohort glm",
                     "cohort glm robvar infl",
                     "cohort glmgee",
                     "cohort glmgee robvar infl",
                     "cohort glmgee robvar dfbeta")

n.method <- 18

details.beta0     <- NULL
details.beta1     <- NULL
details.beta2     <- NULL
details.beta3     <- NULL
for (i in 1:nrow(scenarios)) {
  
  RECAP1 <- RECAP[((i-1) * (n.method * Nreplic) + 1):(i * (n.method * Nreplic)), ]
  
  beta0  <- RECAP1[1,]$beta0
  beta1  <- RECAP1[1,]$beta1
  beta2  <- RECAP1[1,]$beta2
  beta3  <- RECAP1[1,]$beta3
  
  nameColResults <- function (name) {
    a <- NULL
    for (j in 1:length(list.methods)) {
      a <- c(a, paste0(name, ".", colname.methods[j]))
    }
    return(a)
  }
  
  # Bias and mean estimate
  mean.beta0  <- NULL
  mean.beta1  <- NULL
  mean.beta2  <- NULL
  mean.beta3  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    mean.beta0 <- cbind(mean.beta0, mean((RECAP1$beta0.hat[which(RECAP1$Method == method)])))
    mean.beta1 <- cbind(mean.beta1, mean((RECAP1$beta1.hat[which(RECAP1$Method == method)])))
    mean.beta2 <- cbind(mean.beta2, mean((RECAP1$beta2.hat[which(RECAP1$Method == method)])))
    mean.beta3 <- cbind(mean.beta3, mean((RECAP1$beta3.hat[which(RECAP1$Method == method)])))
  }
  mean.beta0 <- as.data.frame(mean.beta0)
  colnames(mean.beta0) <- nameColResults("mean.beta0")
  mean.beta1 <- as.data.frame(mean.beta1)
  colnames(mean.beta1) <- nameColResults("mean.beta1")
  mean.beta2 <- as.data.frame(mean.beta2)
  colnames(mean.beta2) <- nameColResults("mean.beta2")
  mean.beta3 <- as.data.frame(mean.beta3)
  colnames(mean.beta3) <- nameColResults("mean.beta3")
  
  bias.beta0 <- mean.beta0 - beta0
  colnames(bias.beta0) <- nameColResults("bias.beta0")
  bias.beta1 <- mean.beta1 - beta1
  colnames(bias.beta1) <- nameColResults("bias.beta1")
  bias.beta2 <- mean.beta2 - beta2
  colnames(bias.beta2) <- nameColResults("bias.beta2")
  bias.beta3 <- mean.beta3 - beta3
  colnames(bias.beta3) <- nameColResults("bias.beta3")
  
  # Relative bias
  relative.bias.beta0 <- abs(bias.beta0 / beta0)
  colnames(relative.bias.beta0) <- nameColResults("relative.bias.beta0")
  relative.bias.beta1 <- abs(bias.beta1 / beta1)
  colnames(relative.bias.beta1) <- nameColResults("relative.bias.beta1")
  relative.bias.beta2 <- abs(bias.beta2 / beta2)
  colnames(relative.bias.beta2) <- nameColResults("relative.bias.beta2")
  relative.bias.beta3 <- abs(bias.beta3 / beta3)
  colnames(relative.bias.beta3) <- nameColResults("relative.bias.beta3")
  
  # Empirical variance
  empir.var.beta0  <- NULL
  empir.var.beta1  <- NULL
  empir.var.beta2  <- NULL
  empir.var.beta3  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    empir.var.beta0 <- cbind(empir.var.beta0, var(RECAP1$beta0.hat[which(RECAP1$Method == method)]))
    empir.var.beta1 <- cbind(empir.var.beta1, var(RECAP1$beta1.hat[which(RECAP1$Method == method)]))
    empir.var.beta2 <- cbind(empir.var.beta2, var(RECAP1$beta2.hat[which(RECAP1$Method == method)]))
    empir.var.beta3 <- cbind(empir.var.beta3, var(RECAP1$beta3.hat[which(RECAP1$Method == method)]))
  }
  empir.var.beta0 <- as.data.frame(empir.var.beta0)
  colnames(empir.var.beta0) <- nameColResults("empir.var.beta0")
  empir.var.beta1 <- as.data.frame(empir.var.beta1)
  colnames(empir.var.beta1) <- nameColResults("empir.var.beta1")
  empir.var.beta2 <- as.data.frame(empir.var.beta2)
  colnames(empir.var.beta2) <- nameColResults("empir.var.beta2")
  empir.var.beta3 <- as.data.frame(empir.var.beta3)
  colnames(empir.var.beta3) <- nameColResults("empir.var.beta3")
  
  # Mean variance
  mean.var.beta0  <- NULL
  mean.var.beta1  <- NULL
  mean.var.beta2  <- NULL
  mean.var.beta3  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    mean.var.beta0 <- cbind(mean.var.beta0, mean(RECAP1$var.beta0.hat[which(RECAP1$Method == method)]))
    mean.var.beta1 <- cbind(mean.var.beta1, mean(RECAP1$var.beta1.hat[which(RECAP1$Method == method)]))
    mean.var.beta2 <- cbind(mean.var.beta2, mean(RECAP1$var.beta2.hat[which(RECAP1$Method == method)]))
    mean.var.beta3 <- cbind(mean.var.beta3, mean(RECAP1$var.beta3.hat[which(RECAP1$Method == method)]))
  }
  mean.var.beta0 <- as.data.frame(mean.var.beta0)
  colnames(mean.var.beta0) <- nameColResults("mean.var.beta0")
  mean.var.beta1 <- as.data.frame(mean.var.beta1)
  colnames(mean.var.beta1) <- nameColResults("mean.var.beta1")
  mean.var.beta2 <- as.data.frame(mean.var.beta2)
  colnames(mean.var.beta2) <- nameColResults("mean.var.beta2")
  mean.var.beta3 <- as.data.frame(mean.var.beta3)
  colnames(mean.var.beta3) <- nameColResults("mean.var.beta3")
  
  # Relative difference between empirical variance and mean variance estimates
  relative.diffvar.beta0 <- abs(empir.var.beta0 - mean.var.beta0) / mean.var.beta0
  colnames(relative.diffvar.beta0) <- nameColResults("relative.diffvar.beta0")
  relative.diffvar.beta1 <- abs(empir.var.beta1 - mean.var.beta1) / mean.var.beta1
  colnames(relative.diffvar.beta1) <- nameColResults("relative.diffvar.beta1")
  relative.diffvar.beta2 <- abs(empir.var.beta2 - mean.var.beta2) / mean.var.beta2
  colnames(relative.diffvar.beta2) <- nameColResults("relative.diffvar.beta2")
  relative.diffvar.beta3 <- abs(empir.var.beta3 - mean.var.beta3) / mean.var.beta3
  colnames(relative.diffvar.beta3) <- nameColResults("relative.diffvar.beta3")
  
  # Coverage of the CIs
  coverage.beta0  <- NULL
  coverage.beta1  <- NULL
  coverage.beta2  <- NULL
  coverage.beta3  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    coverage.beta0 <- cbind(coverage.beta0, mean((RECAP1$CI.left.beta0.hat[which(RECAP1$Method == method)] < beta0)&(RECAP1$CI.right.beta0.hat[which(RECAP1$Method == method)] > beta0)))
    coverage.beta1 <- cbind(coverage.beta1, mean((RECAP1$CI.left.beta1.hat[which(RECAP1$Method == method)] < beta1)&(RECAP1$CI.right.beta1.hat[which(RECAP1$Method == method)] > beta1)))
    coverage.beta2 <- cbind(coverage.beta2, mean((RECAP1$CI.left.beta2.hat[which(RECAP1$Method == method)] < beta2)&(RECAP1$CI.right.beta2.hat[which(RECAP1$Method == method)] > beta2)))
    coverage.beta3 <- cbind(coverage.beta3, mean((RECAP1$CI.left.beta3.hat[which(RECAP1$Method == method)] < beta3)&(RECAP1$CI.right.beta3.hat[which(RECAP1$Method == method)] > beta3)))
  }
  coverage.beta0 <- as.data.frame(coverage.beta0)
  colnames(coverage.beta0) <- nameColResults("coverage.beta0")
  coverage.beta1 <- as.data.frame(coverage.beta1)
  colnames(coverage.beta1) <- nameColResults("coverage.beta1")
  coverage.beta2 <- as.data.frame(coverage.beta2)
  colnames(coverage.beta2) <- nameColResults("coverage.beta2")
  coverage.beta3 <- as.data.frame(coverage.beta3)
  colnames(coverage.beta3) <- nameColResults("coverage.beta3")
  
  # nominal coverage
  nominal.coverage.beta0  <- NULL
  nominal.coverage.beta1  <- NULL
  nominal.coverage.beta2  <- NULL
  nominal.coverage.beta3  <- NULL
  for (j in 1:length(list.methods)) {
    nominal.coverage.beta0 <- cbind(nominal.coverage.beta0, 
                                    (coverage.beta0[j] > coverage.left) & (coverage.beta0[j] < coverage.right))
    nominal.coverage.beta1 <- cbind(nominal.coverage.beta1, 
                                    (coverage.beta1[j] > coverage.left) & (coverage.beta1[j] < coverage.right))
    nominal.coverage.beta2 <- cbind(nominal.coverage.beta2, 
                                    (coverage.beta2[j] > coverage.left) & (coverage.beta2[j] < coverage.right))
    nominal.coverage.beta3 <- cbind(nominal.coverage.beta3, 
                                    (coverage.beta3[j] > coverage.left) & (coverage.beta3[j] < coverage.right))
  }
  nominal.coverage.beta0 <- as.data.frame(nominal.coverage.beta0)
  colnames(nominal.coverage.beta0) <- nameColResults("nominal.coverage.beta0")
  nominal.coverage.beta1 <- as.data.frame(nominal.coverage.beta1)
  colnames(nominal.coverage.beta1) <- nameColResults("nominal.coverage.beta1")
  nominal.coverage.beta2 <- as.data.frame(nominal.coverage.beta2)
  colnames(nominal.coverage.beta2) <- nameColResults("nominal.coverage.beta2")
  nominal.coverage.beta3 <- as.data.frame(nominal.coverage.beta3)
  colnames(nominal.coverage.beta3) <- nameColResults("nominal.coverage.beta3")
  
  details.beta0   <- rbind(details.beta0, c(mean.beta0, bias.beta0, relative.bias.beta0,
                                            empir.var.beta0,
                                            mean.var.beta0, 
                                            relative.diffvar.beta0,
                                            coverage.beta0,
                                            nominal.coverage.beta0,
                                            beta0 = beta0, 
                                            beta1 = beta1,
                                            beta2 = beta2, 
                                            beta3 = beta3,
                                            n = as.character(RECAP1[1,]$n),  
                                            mean.n.cases = mean(RECAP1$n.cases),
                                            mean.n.casecohorts = mean(RECAP1$n.casecohort), 
                                            mean.strata.n_W0 = mean(RECAP1$strata.n_W0), 
                                            mean.strata.n_W1 = mean(RECAP1$strata.n_W1), 
                                            mean.strata.n.cases_W0 = mean(RECAP1$strata.n.cases_W0),
                                            mean.strata.n.cases_W1 = mean(RECAP1$strata.n.cases_W1),
                                            mean.strata.n.noncases_W0 = mean(RECAP1$strata.n.noncases_W0),
                                            mean.strata.n.noncases_W1 = mean(RECAP1$strata.n.noncases_W1),
                                            mean.strata.m_W0 = mean(RECAP1$strata.m_W0),
                                            mean.strata.m_W1 = mean(RECAP1$strata.m_W1), 
                                            sp = as.character(RECAP1[1,]$sp), 
                                            pY = as.character(RECAP1[1,]$pY)))
  
  
  details.beta1   <- rbind(details.beta1, c(mean.beta1, bias.beta1, relative.bias.beta1,
                                            empir.var.beta1,
                                            mean.var.beta1, 
                                            relative.diffvar.beta1,
                                            coverage.beta1,
                                            nominal.coverage.beta1,
                                            beta0 = beta0, 
                                            beta1 = beta1,
                                            beta2 = beta2, 
                                            beta3 = beta3,
                                            n = as.character(RECAP1[1,]$n),  
                                            mean.n.cases = mean(RECAP1$n.cases),
                                            mean.n.casecohorts = mean(RECAP1$n.casecohort), 
                                            mean.strata.n_W0 = mean(RECAP1$strata.n_W0), 
                                            mean.strata.n_W1 = mean(RECAP1$strata.n_W1), 
                                            mean.strata.n.cases_W0 = mean(RECAP1$strata.n.cases_W0),
                                            mean.strata.n.cases_W1 = mean(RECAP1$strata.n.cases_W1),
                                            mean.strata.n.noncases_W0 = mean(RECAP1$strata.n.noncases_W0),
                                            mean.strata.n.noncases_W1 = mean(RECAP1$strata.n.noncases_W1),
                                            mean.strata.m_W0 = mean(RECAP1$strata.m_W0),
                                            mean.strata.m_W1 = mean(RECAP1$strata.m_W1), 
                                            sp = as.character(RECAP1[1,]$sp), 
                                            pY = as.character(RECAP1[1,]$pY)))
  
  details.beta2   <- rbind(details.beta2, c(mean.beta2, bias.beta2, relative.bias.beta2,
                                            empir.var.beta2,
                                            mean.var.beta2, 
                                            relative.diffvar.beta2,
                                            coverage.beta2,
                                            nominal.coverage.beta2,
                                            beta0 = beta0, 
                                            beta1 = beta1,
                                            beta2 = beta2, 
                                            beta3 = beta3,
                                            n = as.character(RECAP1[1,]$n),  
                                            mean.n.cases = mean(RECAP1$n.cases),
                                            mean.n.casecohorts = mean(RECAP1$n.casecohort), 
                                            mean.strata.n_W0 = mean(RECAP1$strata.n_W0), 
                                            mean.strata.n_W1 = mean(RECAP1$strata.n_W1), 
                                            mean.strata.n.cases_W0 = mean(RECAP1$strata.n.cases_W0),
                                            mean.strata.n.cases_W1 = mean(RECAP1$strata.n.cases_W1),
                                            mean.strata.n.noncases_W0 = mean(RECAP1$strata.n.noncases_W0),
                                            mean.strata.n.noncases_W1 = mean(RECAP1$strata.n.noncases_W1),
                                            mean.strata.m_W0 = mean(RECAP1$strata.m_W0),
                                            mean.strata.m_W1 = mean(RECAP1$strata.m_W1), 
                                            sp = as.character(RECAP1[1,]$sp), 
                                            pY = as.character(RECAP1[1,]$pY)))
  
  details.beta3   <- rbind(details.beta3, c(mean.beta3, bias.beta3, relative.bias.beta3,
                                            empir.var.beta3,
                                            mean.var.beta3, 
                                            relative.diffvar.beta3,
                                            coverage.beta3,
                                            nominal.coverage.beta3,
                                            beta0 = beta0, 
                                            beta1 = beta1,
                                            beta2 = beta2, 
                                            beta3 = beta3,
                                            n = as.character(RECAP1[1,]$n),  
                                            mean.n.cases = mean(RECAP1$n.cases),
                                            mean.n.casecohorts = mean(RECAP1$n.casecohort), 
                                            mean.strata.n_W0 = mean(RECAP1$strata.n_W0), 
                                            mean.strata.n_W1 = mean(RECAP1$strata.n_W1), 
                                            mean.strata.n.cases_W0 = mean(RECAP1$strata.n.cases_W0),
                                            mean.strata.n.cases_W1 = mean(RECAP1$strata.n.cases_W1),
                                            mean.strata.n.noncases_W0 = mean(RECAP1$strata.n.noncases_W0),
                                            mean.strata.n.noncases_W1 = mean(RECAP1$strata.n.noncases_W1),
                                            mean.strata.m_W0 = mean(RECAP1$strata.m_W0),
                                            mean.strata.m_W1 = mean(RECAP1$strata.m_W1), 
                                            sp = as.character(RECAP1[1,]$sp), 
                                            pY = as.character(RECAP1[1,]$pY)))
  
}
details.beta0 <- as.data.frame(details.beta0)
details.beta1 <- as.data.frame(details.beta1)
details.beta2 <- as.data.frame(details.beta2)
details.beta3 <- as.data.frame(details.beta3)

simul.name = paste0("_logistic_stratified-", B, "_boot-", Nreplic, "_Nreplic")

save(details.beta0, file = paste0("details.beta0", simul.name, ".Rdata"))
col.names.beta0 <- colnames(details.beta0[,c(1:ncol(details.beta0))])
details.beta0[col.names.beta0] <- sapply(details.beta0[col.names.beta0], 
                                         as.numeric)
write.csv(details.beta0, file = paste0("details.beta0", simul.name, ".csv"))
write.xlsx(details.beta0, file = paste0("details.beta0", simul.name, ".xlsx"))

save(details.beta1, file = paste0("details.beta1", simul.name, ".Rdata"))
col.names.beta1 <- colnames(details.beta1[,c(1:ncol(details.beta1))])
details.beta1[col.names.beta1] <- sapply(details.beta1[col.names.beta1], 
                                         as.numeric)
write.csv(details.beta1, file = paste0("details.beta1", simul.name, ".csv"))
write.xlsx(details.beta1, file = paste0("details.beta1", simul.name, ".xlsx"))

save(details.beta2, file = paste0("details.beta2", simul.name, ".Rdata"))
col.names.beta2 <- colnames(details.beta2[,c(1:ncol(details.beta2))])
details.beta2[col.names.beta2] <- sapply(details.beta2[col.names.beta2], 
                                         as.numeric)
write.csv(details.beta2, file = paste0("details.beta2", simul.name, ".csv"))
write.xlsx(details.beta2, file = paste0("details.beta2", simul.name, ".xlsx"))

save(details.beta3, file = paste0("details.beta3", simul.name, ".Rdata"))
col.names.beta3 <- colnames(details.beta3[,c(1:ncol(details.beta3))])
details.beta3[col.names.beta3] <- sapply(details.beta3[col.names.beta3], 
                                         as.numeric)
write.csv(details.beta3, file = paste0("details.beta3", simul.name, ".csv"))
write.xlsx(details.beta3, file = paste0("details.beta3", simul.name, ".xlsx"))

# ------------------------------------------------------------------------------
## Printing the Latex tables for the manuscript --------------------------------

### beta0
a.beta0.glm             <- cbind(details.beta0[, c(158:159, 55, 74, 110)])
scientific.glm          <- a.beta0.glm
scientific.glm[,c(3:4)] <- data.frame(lapply(a.beta0.glm[,c(3:4)], 
                                             function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5)]   <- format(a.beta0.glm[,c(5)], digits = 3)
a.beta0.glm             <- data.frame(lapply(scientific.glm, 
                                             function(x) format(x, digits = 3)))

a.beta0.boot            <- cbind(details.beta0[, c(84, 120)]) 
scientific.boot         <- a.beta0.boot
scientific.boot[,c(1)]  <- format(a.beta0.boot[,c(1)], scientific = TRUE, digits = 2)
scientific.boot[,c(2)]  <- format(a.beta0.boot[,c(2)], digits = 3)
a.beta0.boot            <- data.frame(lapply(scientific.boot, 
                                             function(x) format(x, digits = 3)))

a.beta0.survey              <- cbind(details.beta0[, c(63, 81, 117)]) 
scientific.survey           <- a.beta0.survey
scientific.survey[,c(1:2)]  <- data.frame(lapply(a.beta0.survey[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.survey[,3]       <- format(a.beta0.survey[,3], digits = 3)
a.beta0.survey              <- data.frame(lapply(scientific.survey, 
                                                 function(x) format(x, digits = 3)))

a.beta0.glmgee              <- cbind(details.beta0[, c(58, 76, 112)]) 
scientific.glmgee           <- a.beta0.glmgee
scientific.glmgee[,c(1:2)]  <- data.frame(lapply(a.beta0.glmgee[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3)]    <- format(a.beta0.glmgee[,c(3)], digits = 3)
a.beta0.glmgee              <- data.frame(lapply(scientific.glmgee , 
                                                 function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta0.glm, a.beta0.boot, a.beta0.survey, a.beta0.glmgee)), 
      include.rownames = FALSE) 

### beta1
a.beta1.glm             <- cbind(details.beta1[, c(158:159, 55, 74, 110)]) 
scientific.glm          <- a.beta1.glm
scientific.glm[,c(3:4)] <- data.frame(lapply(a.beta1.glm[,c(3:4)], 
                                             function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5)]   <- format(a.beta1.glm[,c(5)], digits = 3)
a.beta1.glm             <- data.frame(lapply(scientific.glm, 
                                             function(x) format(x, digits = 3)))

a.beta1.boot            <- cbind(details.beta1[, c(84, 120)]) 
scientific.boot         <- a.beta1.boot
scientific.boot[,c(1)]  <- format(a.beta1.boot[,c(1)], scientific = TRUE, digits = 2)
scientific.boot[,c(2)]  <- format(a.beta1.boot[,c(2)], digits = 3)
a.beta1.boot            <- data.frame(lapply(scientific.boot, 
                                             function(x) format(x, digits = 3)))

a.beta1.survey              <- cbind(details.beta1[, c(63, 81, 117)]) 
scientific.survey           <- a.beta1.survey
scientific.survey[,c(1:2)]  <- data.frame(lapply(a.beta1.survey[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.survey[,3]       <- format(a.beta1.survey[,3], digits = 3)
a.beta1.survey              <- data.frame(lapply(scientific.survey, 
                                                 function(x) format(x, digits = 3)))

a.beta1.glmgee              <- cbind(details.beta1[, c(58, 76, 112)]) 
scientific.glmgee           <- a.beta1.glmgee
scientific.glmgee[,c(1:2)]  <- data.frame(lapply(a.beta1.glmgee[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3)]    <- format(a.beta1.glmgee[,c(3)], digits = 3)
a.beta1.glmgee              <- data.frame(lapply(scientific.glmgee, 
                                                 function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta1.glm, a.beta1.boot, a.beta1.survey, a.beta1.glmgee)), 
      include.rownames = FALSE) 

### beta2
a.beta2.glm             <- cbind(details.beta2[, c(158:159, 55, 74, 110)]) 
scientific.glm          <- a.beta2.glm
scientific.glm[,c(3:4)] <- data.frame(lapply(a.beta2.glm[,c(3:4)], 
                                             function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5)]   <- format(a.beta2.glm[,c(5)], digits = 3)
a.beta2.glm             <- data.frame(lapply(scientific.glm,
                                             function(x) format(x, digits = 3)))

a.beta2.boot            <- cbind(details.beta2[, c(84, 120)]) 
scientific.boot         <- a.beta2.boot
scientific.boot[,c(1)]  <- format(a.beta2.boot[,c(1)], scientific = TRUE, digits = 2)
scientific.boot[,c(2)]  <- format(a.beta2.boot[,c(2)], digits = 3)
a.beta2.boot            <- data.frame(lapply(scientific.boot, 
                                             function(x) format(x, digits = 3)))

a.beta2.survey              <- cbind(details.beta2[, c(63, 81, 117)]) 
scientific.survey           <- a.beta2.survey
scientific.survey[,c(1:2)]  <- data.frame(lapply(a.beta2.survey[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.survey[,3]       <- format(a.beta2.survey[,3], digits = 3)
a.beta2.survey              <- data.frame(lapply(scientific.survey,
                                                 function(x) format(x, digits = 3)))

a.beta2.glmgee              <- cbind(details.beta2[, c(58, 76, 112)]) 
scientific.glmgee           <- a.beta2.glmgee
scientific.glmgee[,c(1:2)]  <- data.frame(lapply(a.beta2.glmgee[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3)]    <- format(a.beta2.glmgee[,c(3)], digits = 3)
a.beta2.glmgee              <- data.frame(lapply(scientific.glmgee, 
                                                 function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta2.glm, a.beta2.boot, a.beta2.survey, a.beta2.glmgee)), include.rownames = FALSE) 

### beta3
a.beta3.glm             <- cbind(details.beta3[, c(158:159, 55, 74, 110)]) 
scientific.glm          <- a.beta3.glm
scientific.glm[,c(3:4)] <- data.frame(lapply(a.beta3.glm[,c(3:4)], 
                                             function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5)]   <- format(a.beta3.glm[,c(5)], digits = 3)
a.beta3.glm             <- data.frame(lapply(scientific.glm, 
                                             function(x) format(x, digits = 3)))

a.beta3.boot            <- cbind(details.beta3[, c(84, 120)]) 
scientific.boot         <- a.beta3.boot
scientific.boot[,c(1)]  <- format(a.beta3.boot[,c(1)], scientific = TRUE, digits = 2)
scientific.boot[,c(2)]  <- format(a.beta3.boot[,c(2)], digits = 3)
a.beta3.boot            <- data.frame(lapply(scientific.boot,
                                             function(x) format(x, digits = 3)))

a.beta3.survey              <- cbind(details.beta3[, c(63, 81, 117)]) 
scientific.survey           <- a.beta3.survey
scientific.survey[,c(1:2)]  <- data.frame(lapply(a.beta3.survey[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.survey[,3]       <- format(a.beta3.survey[,3], digits = 3)
a.beta3.survey              <- data.frame(lapply(scientific.survey,
                                                 function(x) format(x, digits = 3)))

a.beta3.glmgee              <- cbind(details.beta3[, c(58, 76, 112)]) 
scientific.glmgee           <- a.beta3.glmgee
scientific.glmgee[,c(1:2)]  <- data.frame(lapply(a.beta3.glmgee[,c(1:2)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3)]    <- format(a.beta3.glmgee[,c(3)], digits = 3)
a.beta3.glmgee              <- data.frame(lapply(scientific.glmgee, 
                                                 function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta3.glm, a.beta3.boot, a.beta3.survey, a.beta3.glmgee)), include.rownames = FALSE) 


## -----------------------------------------------------------------------------
## Printing the Latex tables for the supplementary document --------------------

### beta0
a.beta0.glm                 <- cbind(details.beta0[, c(158:159, 55, 73, 109, 75, 111)]) 
scientific.glm              <- a.beta0.glm
scientific.glm[,c(3:4,6)]   <- data.frame(lapply(a.beta0.glm[,c(3:4,6)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5,7)]     <- format(a.beta0.glm[,c(5,7)], digits = 3)
a.beta0.glm                 <- data.frame(lapply(scientific.glm,
                                                 function(x) format(x, digits = 3)))

a.beta0.boot                  <- cbind(details.beta0[, c(84, 120, 121, 82, 118, 119)]) 
scientific.boot               <- a.beta0.boot
scientific.boot[,c(1,4)]      <- format(a.beta0.boot[,c(1,4)], scientific = TRUE, digits = 2)
scientific.boot[,c(2,3,5,6)]  <- format(a.beta0.boot[,c(2,3,5,6)], digits = 3)
a.beta0.boot                  <- data.frame(lapply(scientific.boot,
                                                   function(x) format(x, digits = 3)))

a.beta0.glmgee                <- cbind(details.beta0[, c(58, 76, 112, 77, 113, 78, 114, 79, 115, 80, 116)]) 
scientific.glmgee             <- a.beta0.glmgee
scientific.glmgee[,c(1:2,4,6,8,10)] <- data.frame(lapply(a.beta0.glmgee[,c(1:2,4,6,8,10)], 
                                                         function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3,5,7,9)]      <- format(a.beta0.glmgee[,c(3,5,7,9)], digits = 3)
a.beta0.glmgee                <- data.frame(lapply(scientific.glmgee, 
                                                   function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta0.glm, a.beta0.boot)), include.rownames = FALSE) 
print(xtable(cbind(a.beta0.glm[,1:2], a.beta0.glmgee)), include.rownames = FALSE) 

### beta1
a.beta1.glm                 <- cbind(details.beta1[, c(158:159, 55, 73, 109, 75, 111)]) 
scientific.glm              <- a.beta1.glm
scientific.glm[,c(3:4,6)]   <- data.frame(lapply(a.beta1.glm[,c(3:4,6)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5,7)]     <- format(a.beta1.glm[,c(5,7)], digits = 3)
a.beta1.glm                 <- data.frame(lapply(scientific.glm,
                                                 function(x) format(x, digits = 3)))

a.beta1.boot                  <- cbind(details.beta1[, c(84, 120, 121, 82, 118, 119)]) 
scientific.boot               <- a.beta1.boot
scientific.boot[,c(1,4)]      <- format(a.beta1.boot[,c(1,4)], scientific = TRUE, digits = 2)
scientific.boot[,c(2,3,5,6)]  <- format(a.beta1.boot[,c(2,3,5,6)], digits = 3)
a.beta1.boot                  <- data.frame(lapply(scientific.boot,
                                                   function(x) format(x, digits = 3)))

a.beta1.glmgee                <- cbind(details.beta1[, c(58, 76, 112, 77, 113, 78, 114, 79, 115, 80, 116)]) 
scientific.glmgee             <- a.beta1.glmgee
scientific.glmgee[,c(1:2,4,6,8,10)] <- data.frame(lapply(a.beta1.glmgee[,c(1:2,4,6,8,10)], 
                                                         function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3,5,7,9)]      <- format(a.beta1.glmgee[,c(3,5,7,9)], digits = 3)
a.beta1.glmgee                <- data.frame(lapply(scientific.glmgee, 
                                                   function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta1.glm, a.beta1.boot)), include.rownames = FALSE) 
print(xtable(cbind(a.beta1.glm[,1:2], a.beta1.glmgee)), include.rownames = FALSE) 

### beta2
a.beta2.glm                 <- cbind(details.beta2[, c(158:159, 55, 73, 109, 75, 111)]) 
scientific.glm              <- a.beta2.glm
scientific.glm[,c(3:4,6)]   <- data.frame(lapply(a.beta2.glm[,c(3:4,6)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5,7)]     <- format(a.beta2.glm[,c(5,7)], digits = 3)
a.beta2.glm                 <- data.frame(lapply(scientific.glm,
                                                 function(x) format(x, digits = 3)))

a.beta2.boot                  <- cbind(details.beta2[, c(84, 120, 121, 82, 118, 119)]) 
scientific.boot               <- a.beta2.boot
scientific.boot[,c(1,4)]      <- format(a.beta2.boot[,c(1,4)], scientific = TRUE, digits = 2)
scientific.boot[,c(2,3,5,6)]  <- format(a.beta2.boot[,c(2,3,5,6)], digits = 3)
a.beta2.boot                  <- data.frame(lapply(scientific.boot,
                                                   function(x) format(x, digits = 3)))

a.beta2.glmgee                <- cbind(details.beta2[, c(58, 76, 112, 77, 113, 78, 114, 79, 115, 80, 116)]) 
scientific.glmgee             <- a.beta2.glmgee
scientific.glmgee[,c(1:2,4,6,8,10)] <- data.frame(lapply(a.beta2.glmgee[,c(1:2,4,6,8,10)], 
                                                         function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3,5,7,9)]      <- format(a.beta2.glmgee[,c(3,5,7,9)], digits = 3)
a.beta2.glmgee                <- data.frame(lapply(scientific.glmgee, 
                                                   function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta2.glm, a.beta2.boot)), include.rownames = FALSE) 
print(xtable(cbind(a.beta2.glm[,1:2], a.beta2.glmgee)), include.rownames = FALSE) 

### beta3
a.beta3.glm                 <- cbind(details.beta3[, c(158:159, 55, 73, 109, 75, 111)]) 
scientific.glm              <- a.beta3.glm
scientific.glm[,c(3:4,6)]   <- data.frame(lapply(a.beta3.glm[,c(3:4,6)], 
                                                 function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glm[,c(5,7)]     <- format(a.beta3.glm[,c(5,7)], digits = 3)
a.beta3.glm                 <- data.frame(lapply(scientific.glm,
                                                 function(x) format(x, digits = 3)))

a.beta3.boot                  <- cbind(details.beta3[, c(84, 120, 121, 82, 118, 119)]) 
scientific.boot               <- a.beta3.boot
scientific.boot[,c(1,4)]      <- format(a.beta3.boot[,c(1,4)], scientific = TRUE, digits = 2)
scientific.boot[,c(2,3,5,6)]  <- format(a.beta3.boot[,c(2,3,5,6)], digits = 3)
a.beta3.boot                  <- data.frame(lapply(scientific.boot,
                                                   function(x) format(x, digits = 3)))

a.beta3.glmgee                <- cbind(details.beta3[, c(58, 76, 112, 77, 113, 78, 114, 79, 115, 80, 116)]) 
scientific.glmgee             <- a.beta3.glmgee
scientific.glmgee[,c(1:2,4,6,8,10)] <- data.frame(lapply(a.beta3.glmgee[,c(1:2,4,6,8,10)], 
                                                         function(x) format(x, scientific = TRUE, digits = 2)))
scientific.glmgee[,c(3,5,7,9)]      <- format(a.beta3.glmgee[,c(3,5,7,9)], digits = 3)
a.beta3.glmgee                <- data.frame(lapply(scientific.glmgee, 
                                                   function(x) format(x, digits = 3)))

print(xtable(cbind(a.beta3.glm, a.beta3.boot)), include.rownames = FALSE) 
print(xtable(cbind(a.beta3.glm[,1:2], a.beta3.glmgee)), include.rownames = FALSE) 
