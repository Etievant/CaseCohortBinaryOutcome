## -----------------------------------------------------------------------------
## Function: influences.glm()
## -----------------------------------------------------------------------------
## Description: This function computes the influences on the beta parameter 
##              estimates in a log-linear or logistic regression model, when the 
##              parameter has been estimated with the glm function. This 
##              function can be used under the design with unstratified sampling 
##              of the subcohort or stratified sampling based on some covariate 
##              and exhaustive sampling of cases, or under the design with 
##              stratified sampling based on case status and possibly some
##              covariate
## -----------------------------------------------------------------------------
## Argument:
##
##  mod         a glm model object, result of function glm used with the 
##              binomial family and "logit" or "log" link
## -----------------------------------------------------------------------------
## Values:
##
##  beta.hat    vector with the parameter estimates obtained from the glm model
##  
##  infl.beta   data frame with the influences on the parameter estimates
##
##  link        link used in the function glm, "logit" or "log"
## -----------------------------------------------------------------------------

influences.glm  <- function(mod) {
  X           <- mod$x
  Y           <- mod$y
  q           <- ncol(X)
  n           <- nrow(X)
  weights     <- mod$prior.weights
  
  beta.hat    <- mod$coefficients
  X.weight    <- weights * X
  
  link        <- mod$family$link
  
  if(link == "log"){
    p.hat               <- mod$fitted.values
    XYp.weight          <- (Y - 1) * p.hat / (1 - p.hat)^2  * X.weight
    XXYp.weight         <- array(NA, dim = c(q, q, n))
    for(i in 1:n) { XXYp.weight[,,i] <- tcrossprod(XYp.weight[i,], X[i,]) }
    sum.XXYp.weight     <- apply(X = XXYp.weight, MARGIN = c(1,2), FUN = sum)
    sum.XXYp.weight.inv <- solve(sum.XXYp.weight)
    infl.beta <- ((p.hat - Y) / (1 - p.hat) * X.weight) %*% sum.XXYp.weight.inv
  }
  
  if(link == "logit"){
    cl.hat        <- mod$linear.predictors
    p.hat         <- sigmoid(cl.hat)
    p1p.hat       <- p.hat * (1 - p.hat)
    XX.weight     <- array(NA, dim = c(q, q, n))
    for(i in 1:n) { XX.weight[,,i] <- p1p.hat[i] * tcrossprod(X.weight[i,], X[i,]) }
    sum.XXp.weight      <- apply(X = XX.weight, MARGIN = c(1,2), FUN = sum)
    sum.XXp.weight.inv  <- solve(sum.XXp.weight)
    infl.beta           <- ((p.hat - Y) * X.weight) %*% sum.XXp.weight.inv
  }
  return(list(beta.hat = beta.hat, infl.beta = infl.beta, link = link))
}


## -----------------------------------------------------------------------------
## Function: influences.gee()
## -----------------------------------------------------------------------------
## Description: This function computes the influences on the beta parameter 
##              estimates in a log-linear or logistic regression model, when the 
##              parameter has been estimated with the glmgee function. This 
##              function can be used under the design with unstratified sampling 
##              of the subcohort or stratified sampling based on some covariate 
##              and exhaustive sampling of cases, or under the design with 
##              stratified sampling based on case status and possibly some
##              covariate
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod         a glmgee model object, result of function glmgee used with the
##              binomial family and "logit" or "log" link
##
##  X           design matrix corresponding to one used in the formula 
##              provided to the glmgee function, with columns in the same 
##              order, and including a column for the intercept
## -----------------------------------------------------------------------------
## Values:
##
##  beta.hat    vector with the parameter estimates obtained from the glmgee 
##              model
##  
##  infl.beta   data frame with the influences on the parameter estimates, with
##              as many columns as parameters, i.e., as columns in X
##
##  link        link used in the function glmgee, "logit" or "log"
## -----------------------------------------------------------------------------

influences.gee  <- function(mod, X) {
  Y             <- c(mod$y)
  q             <- ncol(X)
  n             <- nrow(X)
  weights       <- c(mod$prior.weights)
  
  beta.hat      <- mod$coefficients
  X.weight      <- weights * X
  
  link          <- mod$family$link
  
  if (link == "log") { 
    p.hat                 <- c(mod$fitted.values)
    XYp.weight          <- (Y - 1) * p.hat / (1 - p.hat)^2  * X.weight
    XXYp.weight         <- array(NA, dim = c(q, q, n))
    for(i in 1:n) { XXYp.weight[,,i] <- tcrossprod(XYp.weight[i,], X[i,]) }
    sum.XXYp.weight     <- apply(X = XXYp.weight, MARGIN = c(1,2), FUN = sum)
    sum.XXYp.weight.inv <- solve(sum.XXYp.weight)
    
    infl.beta <- ((p.hat - Y) / (1 - p.hat) * X.weight) %*% sum.XXYp.weight.inv
  }
  
  if (link == "logit") {
    cl.hat        <- c(mod$linear.predictors)
    p.hat         <- sigmoid(cl.hat)
    p1p.hat       <- p.hat * (1 - p.hat)
    XX.weight   <- array(NA, dim = c(q, q, n))
    for(i in 1:n) { XX.weight[,,i] <- p1p.hat[i] * tcrossprod(X.weight[i,], X[i,]) }
    sum.XXp.weight      <- apply(X = XX.weight, MARGIN = c(1,2), FUN = sum)
    sum.XXp.weight.inv  <- solve(sum.XXp.weight)
    infl.beta             <- ((p.hat - Y) * X.weight) %*% sum.XXp.weight.inv
  }
  return(list(beta.hat = beta.hat, infl.beta = infl.beta, link = link))
}


## -----------------------------------------------------------------------------
## Function: variance()
## -----------------------------------------------------------------------------
## Description: This function returns the influence-based variance estimate that 
##              follows the complete variance decomposition. This function can 
##              be used with phase-two data obtained from the design with 
##              unstratified sampling of the subcohort or stratified sampling 
##              based on some covariate and exhaustive sampling of cases, or 
##              under the design with stratified sampling based on case status 
##              and possibly some covariate
## -----------------------------------------------------------------------------
## Arguments:
##
## casecohort   data frame with id (identification), weights (design weights),
##              and Y (binary outcome status) for each individual in the phase
##              two data. If argument stratified below is set to FALSE, it 
##              should also contain m (subcohort size) and n (cohort size). If 
##              argument stratified below is set to TRUE, it should also contain 
##              strata (strata indicator), strata.m (subcohort or phase-two 
##              sample strata size), strata.n (cohort strata size), and n 
##              (cohort size)
##
##  infl        data frame with the influences on the parameter estimates
##
## stratified   was sampling stratified? As default, stratified is set to FALSE.
##              If stratified is set to FALSE then it is assumed that all the 
##              cases in the cohort have been sampled (i.e., exhaustive sampling
##              of the cases)
## -----------------------------------------------------------------------------
## Value:
##
##  variance    vector with the variance estimates
## -----------------------------------------------------------------------------

variance <- function(casecohort, infl, stratified = FALSE) {
  
  if (is.null(casecohort$weights)) {
    print("the subjects' weights must be provided")
    break
  }
  if (is.null(casecohort$Y)) {
    print("the subjects' binary outcome status must be provided")
    break
  }
  if (nrow(casecohort) != nrow(infl)) {
    print("the influences for all the subjects in the case-cohort must be 
          provided")}
  if (is.null(casecohort$n)) {
    print("the cohort size must be provided")
    break
  }
  
  if (stratified == FALSE) {
    if (is.null(casecohort$m)) {
      print("the subcohort size must be provided")
      break
    }
    casecohort.noncases   <- which(casecohort$Y == 0)
    casecohort.n.noncases <- length(casecohort.noncases)
    m                     <- casecohort$m[1]
    n                     <- casecohort$n[1]
    m.n                   <- m / n
    
    Product.sigma.w  <- matrix(1 - m.n * (n - 1) / (m - 1), 
                               nrow = casecohort.n.noncases, 
                               ncol = casecohort.n.noncases)
    diag(Product.sigma.w) <- 1 - m.n
    
    phase1.var  <- n / (n - 1) * crossprod((infl / sqrt(casecohort$weights)))
    phase2.var  <- t(infl[casecohort.noncases,]) %*% Product.sigma.w %*% infl[casecohort.noncases,]
    var         <- diag(phase1.var + phase2.var)
  } else {
    if (is.null(casecohort$strata.n) | is.null(casecohort$strata.m)) {
      print("the strata size in the cohort and in the subcohort or phase-two 
            sample must be provided")
      break
    }
    if (is.null(casecohort$strata)) {
      print("the subjects' strata indicator must be provided")
      break
    }
    casecohort.n            <- nrow(casecohort)
    n                       <- casecohort$n[1]
    casecohort.Y            <- casecohort$Y
    casecohort.strata       <- casecohort$strata
    casecohort.strata.m     <- casecohort$strata.m
    casecohort.strata.n     <- casecohort$strata.n
    casecohort.strata.m.n   <- casecohort.strata.m / casecohort.strata.n
    casecohort.weights      <- casecohort$weights
    
    product.sigma.w <- function (indiv) { 
      return((casecohort.weights[indiv] != 1) * 
               (casecohort.strata == casecohort.strata[indiv]) *
               (casecohort.Y[indiv] == casecohort.Y) *
               (1 - casecohort.strata.m.n * 
                  (casecohort.strata.n - 1) / (casecohort.strata.m - 1)))
    }
    Product.sigma.w       <- do.call(rbind, lapply(1:casecohort.n, product.sigma.w))
    diag(Product.sigma.w) <- 1 - 1 / casecohort.weights
    Product.sigma.w[which(is.na(Product.sigma.w))] <- 0
    
    phase1.var  <-  n / (n - 1) * crossprod((infl / sqrt(casecohort.weights)))
    phase2.var  <- t(infl) %*% Product.sigma.w %*% infl
    var         <- diag(phase1.var + phase2.var)
  }
  return(variance = var)
}


## -----------------------------------------------------------------------------
## Function: CI()
## -----------------------------------------------------------------------------
## Description: This function computes the 100*(1 - alpha)% confidence interval 
##              using the quantiles from the Normal distribution
## -----------------------------------------------------------------------------
## Arguments:
##
## b            parameter estimates
##
## v            variance estimates 
##
## alpha        1 minus the desired confidence level. Default value is 0.05
## -----------------------------------------------------------------------------

CI <- function(b, v, alpha = 0.05) {
  return(c(b - qnorm(1 - alpha / 2) * sqrt(v), b + qnorm(1 - alpha / 2) * sqrt(v)))
}


## -----------------------------------------------------------------------------
## Function: robustvariance()
## -----------------------------------------------------------------------------
## Description: This function returns the robust variance estimate, i.e., 
##              computes the sum of the squared influences
## -----------------------------------------------------------------------------
## Argument:
##
##  infl        data frame with the influences on the parameter estimates
## -----------------------------------------------------------------------------
## Value:
##
##  robust.var  vector with the variance estimates
## -----------------------------------------------------------------------------

robustvariance <- function (infl) {
  return(robust.var = diag(crossprod(infl)))
}


## -----------------------------------------------------------------------------
## Function: cohort.analysis()
## -----------------------------------------------------------------------------
## Description: This function runs the analysis on the entire cohort, as in
##              Etievant (2026), with different software and inference methods. 
##              It estimates the beta parameter in a log-linear or logistic 
##              regression model with both functions glm and glmgee. It then
##              computes the influences when using the estimates from both 
##              functions. Two robust variance estimates are finally obtained
##              from the two sets of influences. A robust variance estimate 
##              is also computed on the dfbetas obtained from the glmgee and 
##              dfbeta.glmgee functions. The two variance estimates directly 
##              obtained from the glm and glmgee functions are also returned. 
##              Corresponding 95% confidence intervals are computed, using the
##              quantiles from the Normal distribution
## -----------------------------------------------------------------------------
## Arguments:
##
## cohort       data frame with id (identification), Y (outcome), X1, X2, and X3
##              (covariates) for each individual in the cohort data
##
## link         model link function for the binomial family, either "log" or 
##              "logit"
## -----------------------------------------------------------------------------
## Values:
##
##  res.cohort.glm.var              parameter and variance estimates returned by 
##                                  the glm function, and corresponding 95% 
##                                  confidence interval obtained from the 
##                                  quantiles from the Normal distribution
##
##  res.cohort.glm.robvar.infl      parameter estimate returned by the glm 
##                                  function, robust variance estimate obtained 
##                                  from the influences derived from the glm 
##                                  output with function influences.glm, and 
##                                  corresponding 95% confidence interval 
##                                  obtained from the quantiles from the Normal 
##                                  distribution
##
##  res.cohort.glmgee.var           parameter and variance estimates returned by 
##                                  the glmgee function, and corresponding 95% 
##                                  confidence interval obtained from the 
##                                  quantiles from the Normal distribution
##
##  res.cohort.glm.robvar.infl      parameter estimate returned by the glmgee 
##                                  function, robust variance estimate obtained 
##                                  from the influences derived from the glmgee
##                                  output with function influences.gee, and 
##                                  corresponding 95% confidence interval 
##                                  obtained from the quantiles from the Normal 
##                                  distribution
##
## res.cohort.glmgee.robvar.dfbeta  parameter estimate returned by the glmgee 
##                                  function, robust variance estimate obtained 
##                                  from the dfbetas returned by the 
##                                  dfbeta.glmgee function, and corresponding 
##                                  95% confidence interval obtained from the 
##                                  quantiles from the Normal distribution
## -----------------------------------------------------------------------------

cohort.analysis <- function(cohort, link) {
  
  if (link == "log") {
    # parameter and variance estimate returned by the glm function
    mod1          <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                         data = cohort)
    mod.cohort    <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                         start = coef(mod1), data = cohort, x = TRUE, y = TRUE)
    
    beta.glm.cohort   <- mod.cohort$coefficients
    var.glm.cohort    <- diag(summary(mod.cohort)$cov.unscaled)
    
    # using the output from the glm function, derivation of the influences by 
    # hand and computation of the robust variance estimate
    est.glm.cohort          <- influences.glm(mod = mod.cohort)
    robvar.infl.glm.cohort  <- robustvariance(est.glm.cohort$infl.beta) 
    
    # parameter and variance returned by the glmgee function 
    mod.gee.cohort      <- glmgee(Y ~ X1 + X2 + X3, 
                                  family = binomial(link = "log"),
                                  data = cohort, start = coef(mod1))
    
    beta.glmgee.cohort  <- c(mod.gee.cohort$coefficients)
    var.glmgee.cohort   <- diag(mod.gee.cohort$R)
    
    # using the output from the glmgee function, derivation of the influences by 
    # hand and computation of the robust variance estimate
    est.glmgee.cohort         <- influences.gee(mod = mod.gee.cohort, 
                                                X = cbind(1, cohort$X1, 
                                                          cohort$X2, 
                                                          cohort$X3))
    robvar.infl.glmgee.cohort <- robustvariance(est.glmgee.cohort$infl.beta) 
    
    # computation of the robust variance estimate on the dfbetas obtained from 
    # the glmgee output and the dfbeta.glmgee function
    dfbs.cohort                   <- dfbeta(mod.gee.cohort, method = "full")
    robvar.dfbeta.glmgee.cohort   <- robustvariance(dfbs.cohort) 
  }
  
  if (link == "logit") {
    # parameter and variance estimate returned by the glm function
    mod.cohort        <- glm(Y ~ X1 + X2 + X3, 
                             family = binomial(link = "logit"), 
                             data = cohort, x = TRUE, y = TRUE)
    
    beta.glm.cohort   <- mod.cohort$coefficients
    var.glm.cohort    <- diag(summary(mod.cohort)$cov.unscaled)
    
    # using the output from the glm function, derivation of the influences by 
    # hand and computation of the robust variance estimate
    est.glm.cohort          <- influences.glm(mod = mod.cohort)
    robvar.infl.glm.cohort  <- robustvariance(est.glm.cohort$infl.beta) 
    
    # parameter and variance returned by the glmgee function 
    mod.gee.cohort      <- glmgee(Y ~ X1 + X2 + X3, 
                                  family = binomial(link = "logit"), 
                                  data = cohort)
    
    beta.glmgee.cohort  <- c(mod.gee.cohort$coefficients)
    var.glmgee.cohort   <- diag(mod.gee.cohort$R)
    
    # using the output from the glmgee function, derivation of the influences by 
    # hand and computation of the robust variance estimate
    est.glmgee.cohort         <- influences.gee(mod = mod.gee.cohort, 
                                                X = cbind(1, cohort$X1, 
                                                          cohort$X2, 
                                                          cohort$X3))
    robvar.infl.glmgee.cohort <- robustvariance(est.glmgee.cohort$infl.beta) 
    
    # computation of the robust variance estimate on the dfbetas obtained from 
    # the glmgee output and the dfbeta.glmgee function
    dfbs.cohort                   <- dfbeta(mod.gee.cohort, method = "full")
    robvar.dfbeta.glmgee.cohort   <- robustvariance(dfbs.cohort) 
  }
  
  res.cohort.glm.var              <- cbind("cohort glm", t(beta.glm.cohort), 
                                           t(var.glm.cohort),
                                           t(CI(beta.glm.cohort, var.glm.cohort)))
  res.cohort.glm.robvar.infl      <- cbind("cohort glm robvar infl", t(beta.glm.cohort), 
                                           t(robvar.infl.glm.cohort),
                                           t(CI(beta.glm.cohort, robvar.infl.glm.cohort)))
  res.cohort.glmgee.var           <- cbind("cohort glmgee", t(beta.glmgee.cohort), 
                                           t(var.glmgee.cohort),
                                           t(CI(beta.glmgee.cohort, var.glmgee.cohort)))
  res.cohort.glmgee.robvar.infl   <- cbind("cohort glmgee robvar infl", 
                                           t(beta.glmgee.cohort), 
                                           t(robvar.infl.glmgee.cohort),
                                           t(CI(beta.glmgee.cohort, robvar.infl.glmgee.cohort)))
  res.cohort.glmgee.robvar.dfbeta <- cbind("cohort glmgee robvar dfbeta", 
                                           t(beta.glmgee.cohort), 
                                           t(robvar.dfbeta.glmgee.cohort),
                                           t(CI(beta.glmgee.cohort, robvar.dfbeta.glmgee.cohort)))
  
  return(list(res.cohort.glm.var = res.cohort.glm.var, 
              res.cohort.glm.robvar.infl = res.cohort.glm.robvar.infl, 
              res.cohort.glmgee.var = res.cohort.glmgee.var,
              res.cohort.glmgee.robvar.infl = res.cohort.glmgee.robvar.infl,
              res.cohort.glmgee.robvar.dfbeta = res.cohort.glmgee.robvar.dfbeta))
}


## -----------------------------------------------------------------------------
## Function: casecohort.analysis()
## -----------------------------------------------------------------------------
## Description: This function runs the analysis on the case-cohort as in
##              Etievant (2026), with different software and inference methods. 
##              It estimates the beta parameter in a (weighted) log-linear or 
##              logistic regression model with both functions glm and glmgee. It 
##              then computes the influences when using the estimates from both 
##              functions. Robust variance and our proposed variance estimates
##              are then obtained from these two sets of influences. The two 
##              variance estimates are also computed on the dfbetas obtained 
##              from the glmgee and dfbeta.glmgee functions. The two variance 
##              estimates directly obtained from the glm and glmgee functions
##              are also returned. Corresponding 95% confidence intervals are 
##              computed, using the quantiles from the Normal distribution
## -----------------------------------------------------------------------------
## Arguments:
##
## casecohort   data frame with id (identification), Y (outcome), X1, X2, and X3
##              (covariates), weights (design weights) for each individual in 
##              the case-cohort data. If argument stratified below is set to 
##              FALSE, it should also contain m (subcohort size) and n 
##              (cohort size). If argument stratified below is set to TRUE, it 
##              should also contain W (strata indicator), strata.m (subcohort 
##              strata size), strata.n (cohort strata size), and n (cohort size)
##
## link         model link function for the binomial family, either "log" or 
##              "logit"
##
## stratified   was sampling stratified? As default, stratified is set to FALSE.
##              If stratified is set to FALSE then it is assumed that all the 
##              cases in the cohort have been sampled (i.e., exhaustive sampling
##              of the cases)
## -----------------------------------------------------------------------------
## Values:
##
##  res.casecohort.glm.var              parameter and variance estimates 
##                                      returned by the glm function, and 
##                                      corresponding 95% confidence interval 
##                                      obtained from the quantiles from the 
##                                      Normal distribution
##
##  res.casecohort.glm.robvar.infl      parameter estimate returned by the glm 
##                                      function, robust variance estimate 
##                                      obtained from the influences derived 
##                                      from the glm output with function 
##                                      influences.glm, and corresponding 95% 
##                                      confidence interval obtained from the 
##                                      quantiles from the Normal distribution
##
##  res.casecohort.glm.var.infl         parameter estimate returned by the glm 
##                                      function, proposed variance estimate 
##                                      obtained from the influences derived 
##                                      from the glm output with function 
##                                      influences.glm, and corresponding 95% 
##                                      confidence interval obtained from the 
##                                      quantiles from the Normal distribution
##
##  res.casecohort.glmgee.var           parameter and variance estimates 
##                                      returned by the glmgee function, and 
##                                      corresponding 95% confidence interval 
##                                      obtained from the quantiles from the 
##                                      Normal distribution
##
##
##  res.casecohort.glmgee.robvar.infl   parameter estimate returned by the 
##                                      glmgee function, robust variance 
##                                      estimate obtained from the influences 
##                                      derived from the glmgee output with 
##                                      function influences.gee, and 
##                                      corresponding 95% confidence interval 
##                                      obtained from the quantiles from the 
##                                      Normal distribution
##
##  res.casecohort.glmgee.var.infl      parameter estimate returned by the 
##                                      glmgee function, proposed variance 
##                                      estimate obtained from the influences 
##                                      derived from the glmgee output with 
##                                      function influences.gee, and 
##                                      corresponding 95% confidence interval 
##                                      obtained from the quantiles from the 
##                                      Normal distribution
##
## res.casecohort.glmgee.robvar.dfbeta  parameter estimate returned by the 
##                                      glmgee function, robust variance 
##                                      estimate obtained from the dfbetas 
##                                      returned by the dfbeta.glmgee function, 
##                                      and corresponding 95% confidence 
##                                      interval obtained from the quantiles 
##                                      from the Normal distribution
##
## res.casecohort.glmgee.var.dfbeta     parameter estimate returned by the 
##                                      glmgee function, proposed variance 
##                                      estimate obtained from the dfbetas 
##                                      returned by the dfbeta.glmgee function, 
##                                      and corresponding 95% confidence 
##                                      interval obtained from the quantiles 
##                                      from the Normal distribution
## -----------------------------------------------------------------------------

casecohort.analysis <- function(casecohort, link, stratified = FALSE) {
  
  if (link == "log") {
    # parameter and variance estimate returned by the glm function
    mod1                <- glm(Y ~ X1 + X2 + X3, 
                               family = poisson(link = "log"), 
                               data = casecohort, weights = weights)
    mod.casecohort      <- glm(Y ~ X1 + X2 + X3, 
                               family = binomial(link = "log"), 
                               data = casecohort, weights = weights, 
                               start = coef(mod1), x = TRUE, y = TRUE)
    
    beta.glm.casecohort <- mod.casecohort$coefficients
    var.glm.casecohort  <- diag(summary(mod.casecohort)$cov.unscaled)
    
    # using the output from the glm function, derivation of the influences by 
    # hand and computation of the two variance estimates
    est.casecohort              <- influences.glm(mod = mod.casecohort)
    robvar.infl.glm.casecohort  <- robustvariance(est.casecohort$infl.beta) 
    var.infl.glm.casecohort     <- variance(casecohort, est.casecohort$infl.beta, 
                                            stratified) 
    
    # parameter and variance returned by the glmgee function 
    mod1                    <- glmgee(Y ~ X1 + X2 + X3, 
                                      family = poisson(link = "log"), 
                                      data = casecohort, weights = weights)
    mod.gee.casecohort      <- glmgee(Y ~ X1 + X2 + X3, 
                                      family = binomial(link = "log"), 
                                      data = casecohort, weights = weights, 
                                      start = coef(mod1))
    
    beta.glmgee.casecohort  <- c(mod.gee.casecohort$coefficients)
    var.glmgee.casecohort   <- diag(vcov(mod.gee.casecohort, 
                                         type = "bias-corrected"))
    
    # using the output from the glmgee function, derivation of the influences by 
    # hand and computation of the two variance estimates
    est.glmgee.casecohort         <- influences.gee(mod = mod.gee.casecohort, 
                                                    X = cbind(1, casecohort$X1, 
                                                              casecohort$X2, 
                                                              casecohort$X3))
    
    robvar.infl.glmgee.casecohort <- robustvariance(est.glmgee.casecohort$infl.beta) 
    var.infl.glmgee.casecohort    <- variance(casecohort, 
                                              est.glmgee.casecohort$infl.beta, 
                                              stratified) 
    
    # computation of the two variance estimates on the dfbetas obtained from 
    # the glmgee output and the dfbeta.glmgee function
    dfbs.casecohort                 <- dfbeta(mod.gee.casecohort, 
                                              method = "full")
    robvar.dfbeta.glmgee.casecohort <- robustvariance(dfbs.casecohort) 
    var.dfbeta.glmgee.casecohort    <- variance(casecohort, dfbs.casecohort, 
                                                stratified) 
  }
  
  if(link == "logit") {
    # parameter and variance estimate returned by the glm function
    mod.casecohort      <- glm(Y ~ X1 + X2 + X3, 
                               family = binomial(link = "logit"), 
                               data = casecohort, weights = weights, x = TRUE, 
                               y = TRUE)
    
    beta.glm.casecohort <- mod.casecohort$coefficients
    var.glm.casecohort  <- diag(summary(mod.casecohort)$cov.unscaled)
    
    # using the output from the glm function, derivation of the influences by hand 
    # and computation of the two variance estimates
    est.casecohort              <- influences.glm(mod = mod.casecohort)
    robvar.infl.glm.casecohort  <- robustvariance(est.casecohort$infl.beta) 
    var.infl.glm.casecohort     <- variance(casecohort, est.casecohort$infl.beta, 
                                            stratified) 
    
    # parameter and variance returned by the glmgee function 
    mod.gee.casecohort      <- glmgee(Y ~ X1 + X2 + X3, 
                                      family = binomial(link = "logit"), 
                                      data = casecohort, weights = weights)
    beta.glmgee.casecohort  <- c(mod.gee.casecohort$coefficients)
    var.glmgee.casecohort   <- diag(vcov(mod.gee.casecohort, 
                                         type = "bias-corrected"))
    
    # using the output from the glmgee function, derivation of the influences by 
    # hand and computation of the two variance estimates
    est.glmgee.casecohort         <- influences.gee(mod = mod.gee.casecohort, 
                                                    X = cbind(1, casecohort$X1, 
                                                              casecohort$X2, 
                                                              casecohort$X3))
    robvar.infl.glmgee.casecohort <- robustvariance(est.glmgee.casecohort$infl.beta) 
    var.infl.glmgee.casecohort    <- variance(casecohort, 
                                              est.glmgee.casecohort$infl.beta, 
                                              stratified) 
    
    # computation of the two variance estimates on the dfbetas obtained from 
    # the glmgee output and the dfbeta.glmgee function
    dfbs.casecohort                 <- dfbeta(mod.gee.casecohort, 
                                              method = "full")
    robvar.dfbeta.glmgee.casecohort <- robustvariance(dfbs.casecohort) 
    var.dfbeta.glmgee.casecohort    <- variance(casecohort, dfbs.casecohort, 
                                                stratified) 
  }
  
  res.casecohort.glm.var              <- cbind("casecohort glm", t(beta.glm.casecohort), 
                                               t(var.glm.casecohort),
                                               t(CI(beta.glm.casecohort, var.glm.casecohort)))
  res.casecohort.glm.var.infl         <- cbind("casecohort glm var infl", t(beta.glm.casecohort), 
                                               t(var.infl.glm.casecohort),
                                               t(CI(beta.glm.casecohort, var.infl.glm.casecohort)))
  res.casecohort.glm.robvar.infl      <- cbind("casecohort glm robvar infl", t(beta.glm.casecohort), 
                                               t(robvar.infl.glm.casecohort),
                                               t(CI(beta.glm.casecohort, robvar.infl.glm.casecohort)))
  res.casecohort.glmgee.var           <- cbind("casecohort glmgee", 
                                               t(beta.glmgee.casecohort), 
                                               t(var.glmgee.casecohort),
                                               t(CI(beta.glmgee.casecohort, var.glmgee.casecohort)))
  res.casecohort.glmgee.var.infl      <- cbind("casecohort glmgee var infl", 
                                               t(beta.glmgee.casecohort), 
                                               t(var.infl.glmgee.casecohort),
                                               t(CI(beta.glmgee.casecohort, var.infl.glmgee.casecohort)))
  res.casecohort.glmgee.robvar.infl   <- cbind("casecohort glmgee robvar infl", 
                                               t(beta.glmgee.casecohort), 
                                               t(robvar.infl.glmgee.casecohort),
                                               t(CI(beta.glmgee.casecohort, robvar.infl.glmgee.casecohort)))
  res.casecohort.glmgee.var.dfbeta    <- cbind("casecohort glmgee var dfbeta", 
                                               t(beta.glmgee.casecohort), 
                                               t(var.dfbeta.glmgee.casecohort),
                                               t(CI(beta.glmgee.casecohort, var.dfbeta.glmgee.casecohort)))
  res.casecohort.glmgee.robvar.dfbeta <- cbind("casecohort glmgee robvar dfbeta", 
                                               t(beta.glmgee.casecohort), 
                                               t(robvar.dfbeta.glmgee.casecohort),
                                               t(CI(beta.glmgee.casecohort, robvar.dfbeta.glmgee.casecohort)))
  
  return(list(res.casecohort.glm.var = res.casecohort.glm.var, 
              res.casecohort.glm.var.infl = res.casecohort.glm.var.infl, 
              res.casecohort.glm.robvar.infl = res.casecohort.glm.robvar.infl, 
              res.casecohort.glmgee.var = res.casecohort.glmgee.var,
              res.casecohort.glmgee.var.infl = res.casecohort.glmgee.var.infl,
              res.casecohort.glmgee.robvar.infl = res.casecohort.glmgee.robvar.infl,
              res.casecohort.glmgee.var.dfbeta = res.casecohort.glmgee.var.dfbeta,
              res.casecohort.glmgee.robvar.dfbeta = res.casecohort.glmgee.robvar.dfbeta))
}


## -----------------------------------------------------------------------------
## Function: casecohort.boot.huang()
## -----------------------------------------------------------------------------
## Description: This function runs the analysis on the case-cohort with 
##              bootstrap variance estimation based on Huang's scheme 
##              (Biometrika, 2014) as in Etievant (2026). It estimates the beta 
##              parameter in a (weighted) log-linear or logistic regression
##              model with function glm. Two 95% confidence intervals are 
##              computed, using the quantiles from the Normal distribution and 
##              the standard percentile method
## -----------------------------------------------------------------------------
## Arguments:
##
## casecohort   data frame with id (identification), Y (outcome), X1, X2, and X3
##              (covariates), weights (design weights) for each individual in 
##              the case-cohort data
##
## B            number of bootstrap samples. Default value is 200
##
## link         model link function for the binomial family, either "log" or 
##              "logit"
## -----------------------------------------------------------------------------
## Values:
##
##  res.casecohort.boot.percentile  parameter estimate returned by the glm
##                                  function, variance estimate obtained from 
##                                  the bootstrap scheme, and 95% percentile 
##                                  confidence interval 
##
##  res.casecohort.boot.wald        parameter estimate returned by the glm
##                                  function, variance estimate obtained from 
##                                  the bootstrap scheme, and 95% confidence 
##                                  interval obtained from the quantiles from 
##                                  the Normal distribution
## -----------------------------------------------------------------------------

casecohort.boot.huang <- function(casecohort, B = 200, link) {
  
  if (link == "log") {
    # parameter estimate returned by the glm function
    mod1            <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                           data = casecohort, weights = weights)
    mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                           data = casecohort, weights = weights, 
                           start = coef(mod1), x = TRUE, y = TRUE)
    beta.glm        <- mod.casecohort$coefficients
    
    # Huang's bootstrap scheme
    b.boot <- NULL
    for (b in 1:B) {
      zeta        <- rexp(nrow(casecohort))
      w           <- casecohort$weights * zeta
      mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                         data = casecohort, weights = w)
      trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                             data = casecohort, weights = w, start = coef(mod1),
                             x = TRUE, y = TRUE), 
                         silent = TRUE)
      
      while(inherits(trial.boot, "try-error")){
        zeta        <- rexp(nrow(casecohort))
        w           <- casecohort$weights * zeta
        mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                           data = casecohort, weights = w)
        trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                               data = casecohort, weights = w, 
                               start = coef(mod1), x = TRUE, y = TRUE), 
                           silent = TRUE)
      }
      mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                             data = casecohort, weights = w, 
                             start = coef(mod1), x = TRUE, y = TRUE)
      b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
    }
  }
  if (link == "logit") {
    # parameter and variance estimate returned by the glm function
    mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                           data = casecohort, weights = weights, x = TRUE, 
                           y = TRUE)
    beta.glm        <- mod.casecohort$coefficients
    
    # Huang's bootstrap scheme
    b.boot <- NULL
    for (b in 1:B) {
      zeta          <- rexp(nrow(casecohort))
      w             <- casecohort$weights * zeta
      trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                               family = binomial(link = "logit"), 
                               data = casecohort, weights = w, x = TRUE, 
                               y = TRUE), silent = TRUE)
      
      while(inherits(trial.boot, "try-error")){
        zeta          <- rexp(nrow(casecohort))
        w             <- casecohort$weights * zeta
        trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                                 family = binomial(link = "logit"), 
                                 data = casecohort, weights = w, x = TRUE, 
                                 y = TRUE), silent = TRUE)
      }
      mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                             data = casecohort, weights = w, x = TRUE, y = TRUE)
      b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
    }
  }
  v.boot          <- apply(b.boot, MARGIN = 2, FUN = var)
  ci.left.boot    <- apply(b.boot, MARGIN = 2, FUN = CI.left.boot)
  ci.right.boot   <- apply(b.boot, MARGIN = 2, FUN = CI.right.boot)
  
  res.casecohort.boot.wald        <- c("casecohort Huang boot wald", beta.glm, 
                                       v.boot, t(CI(beta.glm, v.boot)))
  res.casecohort.boot.percentile  <- c("casecohort Huang boot percentile", 
                                       beta.glm, v.boot, ci.left.boot, 
                                       ci.right.boot)
  
  return(list(res.casecohort.boot.wald = res.casecohort.boot.wald, 
              res.casecohort.boot.percentile = res.casecohort.boot.percentile))
}



## -----------------------------------------------------------------------------
## Function: casecohort.boot.wacholder()
## -----------------------------------------------------------------------------
## Description: This function runs the analysis on the case-cohort with 
##              bootstrap variance estimation based on Wacholder's scheme 
##              (Biometrika, 1989) as in Etievant (2026). It estimates the beta 
##              parameter in a log-linear or logistic regression model with 
##              function glm. Two 95% confidence intervals are computed, using 
##              the quantiles from the Normal distribution and the standard 
##              percentile method
## -----------------------------------------------------------------------------
## Arguments:
##
## casecohort   data frame with id (identification), Y (outcome), X1, X2, and X3
##              (covariates), weights (design weights) for each individual in 
##              the case-cohort data
##
## B            number of bootstrap samples. Default value is 200
##
## link         model link function for the binomial family, either "log" or 
##              "logit"
##
## stratified   was sampling stratified? As default, stratified is set to FALSE.
##              If stratified is set to FALSE then argument exhaustive below is
##              set to TRUE, i.e., it is assumed that all the cases in the 
##              cohort have been sampled
##
## exhaustive   was the sampling of the cases exhaustive, i.e., where all the 
##              cases included in the case cohort? As default, exhaustive is set
##              to TRUE. If argument stratified above is set to FALSE, argument
##              exhaustive is automatically set to TRUE
## -----------------------------------------------------------------------------
## Values:
##
##  res.casecohort.boot.percentile    parameter estimate returned by the glm
##                                    function, variance estimate obtained from 
##                                    the bootstrap scheme, and 95% percentile 
##                                    confidence interval 
##
##  res.casecohort.boot.wald          parameter estimate returned by the glm
##                                    function, variance estimate obtained from 
##                                    the bootstrap scheme, and 95% confidence 
##                                    interval obtained from the quantiles from
##                                    the Normal distribution
##
## res.casecohort.boot.percentile.nw  when using link = "logit", parameter 
##                                    estimate returned by the glm function 
##                                    without using the design weights, variance 
##                                    estimate obtained from the bootstrap 
##                                    scheme, and 95% percentile confidence 
##                                    interval 
##
## res.casecohort.boot.wald.nw        when using link = "logit", parameter 
##                                    estimate returned by the glm function 
##                                    without using the design weights, variance 
##                                    estimate obtained from the bootstrap 
##                                    scheme, and 95% confidence interval 
##                                    obtained from the quantiles from the 
##                                    Normal distribution
## -----------------------------------------------------------------------------

casecohort.boot.wacholder <- function(casecohort, B = 200, link, 
                                      stratified = FALSE, exhaustive = TRUE,
                                      noweight = FALSE) {
  
  if(link == "log"){noweight = FALSE}
  
  if(stratified == FALSE){
    cases       <- which(casecohort$Y == 1)
    noncases    <- which(casecohort$Y == 0)
    n.cases     <- sum(casecohort$Y == 1)
    n.noncases  <- sum(casecohort$Y == 0)
    
    if (link == "log") {
      # parameter estimate returned by the glm function
      mod1            <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                             data = casecohort, weights = weights)
      mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                             data = casecohort, weights = weights, 
                             start = coef(mod1), x = TRUE, y = TRUE)
      beta.glm        <- mod.casecohort$coefficients
      
      # Wacholder's bootstrap scheme
      b.boot <- NULL
      for (b in 1:B) {
        cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
        noncase.boot    <- sample(noncases, size = n.noncases, replace = TRUE)
        casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
        
        mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                           data = casecohort.boot, weights = weights)
        trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                               data = casecohort.boot, weights = weights, 
                               start = coef(mod1), x = TRUE, y = TRUE), 
                           silent = TRUE)
        
        while(inherits(trial.boot, "try-error")){
          cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
          noncase.boot    <- sample(noncases, size = n.noncases, replace = TRUE)
          casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
          
          mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                             data = casecohort.boot, weights = weights)
          trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                                 data = casecohort.boot, weights = weights, 
                                 start = coef(mod1), x = TRUE, y = TRUE), 
                             silent = TRUE)
        }
        
        mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                               data = casecohort.boot, weights = weights, 
                               start = coef(mod1), x = TRUE, y = TRUE)
        b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
      }
    }
    if (link == "logit") {
      # parameter estimate returned by the glm function
      mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                             data = casecohort, weights = weights, x = TRUE, 
                             y = TRUE)
      beta.glm        <- mod.casecohort$coefficients
      
      # also running unweighted logistic regression
      if(noweight == TRUE){
        mod.casecohort.nw  <- glm(Y ~ X1 + X2 + X3, 
                                  family = binomial(link = "logit"), 
                               data = casecohort, x = TRUE, 
                               y = TRUE)
        beta.glm.nw        <- mod.casecohort.nw$coefficients
      }
      
      # Wacholder's bootstrap scheme
      b.boot <- NULL
      b.boot.nw <- NULL
      for (b in 1:B) {
        cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
        noncase.boot    <- sample(noncases, size = n.noncases, replace = TRUE)
        casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
        
        trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                                 family = binomial(link = "logit"), 
                                 data = casecohort.boot, weights = weights, 
                                 x = TRUE, y = TRUE), silent = TRUE)
        
        while(inherits(trial.boot, "try-error")){
          cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
          noncase.boot    <- sample(noncases, size = n.noncases, replace = TRUE)
          casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
          
          trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                                   family = binomial(link = "logit"), 
                                   data = casecohort.boot, weights = weights,
                                   x = TRUE, y = TRUE), silent = TRUE)
        }
        
        mod.casecohort  <- glm(Y ~ X1 + X2 + X3, 
                               family = binomial(link = "logit"), 
                               data = casecohort.boot, weights = weights, 
                               x = TRUE, y = TRUE)
        b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
        
        if(noweight == TRUE){
          mod.casecohort.nw  <- glm(Y ~ X1 + X2 + X3, 
                                    family = binomial(link = "logit"), 
                                    data = casecohort.boot, x = TRUE, y = TRUE)
          b.boot.nw          <- rbind(b.boot.nw, c(mod.casecohort.nw$coefficients))
        }
      }
    }
  }else{
    if(exhaustive == TRUE){
      strata <- levels(as.factor(casecohort$strata))
      strata.n.noncases <- integer(length(strata))
      strata.noncases   <- vector("list", length(strata))
      for (i in seq_along(strata)) {
        strata.noncases[[i]] <- which(casecohort$strata == strata[i] & casecohort$Y == 0)
        strata.n.noncases[i] <- sum(casecohort$strata == strata[i] & casecohort$Y == 0)
      }
      n.cases         <- sum(casecohort$Y == 1)
      cases           <- which(casecohort$Y == 1)
      
      if (link == "log") {
        # parameter estimate returned by the glm function
        mod1            <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                               data = casecohort, weights = weights)
        mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                               data = casecohort, weights = weights, 
                               start = coef(mod1), x = TRUE, y = TRUE)
        beta.glm        <- mod.casecohort$coefficients
        
        # Wacholder's bootstrap scheme
        b.boot <- NULL
        for (b in 1:B) {
          noncase.boot <- NULL
          for (i in seq_along(strata)) {
            noncase.boot <- c(noncase.boot, sample(strata.noncases[[i]], size = strata.n.noncases[i], replace = TRUE))
          }
          cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
          casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
          mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                             data = casecohort.boot, weights = weights)
          trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                                 data = casecohort.boot, weights = weights, 
                                 start = coef(mod1), x = TRUE, y = TRUE), 
                             silent = TRUE)
          
          while(inherits(trial.boot, "try-error")){
            noncase.boot <- NULL
            for (i in seq_along(strata)) {
              noncase.boot <- c(noncase.boot, sample(strata.noncases[[i]], size = strata.n.noncases[i], replace = TRUE))
            }
            cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
            casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
            
            mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                               data = casecohort.boot, weights = weights)
            trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                                   data = casecohort.boot, weights = weights, 
                                   start = coef(mod1), x = TRUE, y = TRUE), 
                               silent = TRUE)
          }
          mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                                 data = casecohort.boot, weights = weights, 
                                 start = coef(mod1), x = TRUE, y = TRUE)
          b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
        }
      }
      
      if (link == "logit") {
        # parameter estimate returned by the glm function
        mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                               data = casecohort, weights = weights, x = TRUE, 
                               y = TRUE)
        beta.glm        <- mod.casecohort$coefficients
        
        # also running unweighted logistic regression
        if(noweight == TRUE){
          mod.casecohort.nw  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                                 data = casecohort, x = TRUE, 
                                 y = TRUE)
          beta.glm.nw        <- mod.casecohort.nw$coefficients
        }
        
        # Wacholder's bootstrap scheme
        b.boot <- NULL
        b.boot.nw <- NULL
        for (b in 1:B) {
          noncase.boot <- NULL
          for (i in seq_along(strata)) {
            noncase.boot <- c(noncase.boot, sample(strata.noncases[[i]], size = strata.n.noncases[i], replace = TRUE))
          }
          cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
          casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
          trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                                   family = binomial(link = "logit"), 
                                   data = casecohort.boot, weights = weights, 
                                   x = TRUE, y = TRUE), silent = TRUE)
          
          while(inherits(trial.boot, "try-error")){
            noncase.boot <- NULL
            for (i in seq_along(strata)) {
              noncase.boot <- c(noncase.boot, sample(strata.noncases[[i]], size = strata.n.noncases[i], replace = TRUE))
            }
            cases.boot      <- sample(cases, size = n.cases, replace = TRUE)
            casecohort.boot <- casecohort[c(cases.boot, noncase.boot), ]
            
            trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                                     family = binomial(link = "logit"), 
                                     data = casecohort.boot, weights = weights, 
                                     x = TRUE, y = TRUE), silent = TRUE)
          }
          mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                                 data = casecohort.boot, weights = weights, 
                                 x = TRUE, y = TRUE)
          b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
          if(noweight == TRUE){
            mod.casecohort.nw  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                                   data = casecohort.boot, 
                                   x = TRUE, y = TRUE)
            b.boot.nw          <- rbind(b.boot.nw, c(mod.casecohort.nw$coefficients))
          }
        }
      }
    }else{
      strata <- levels(as.factor(casecohort$strata))
      strata.indiv <- vector("list", length(strata))
      strata.n.indiv <- integer(length(strata))
      for (i in seq_along(strata)) {
        strata.indiv[[i]] <- which(casecohort$strata == strata[i])
        strata.n.indiv[i] <- sum(casecohort$strata == strata[i])
      }
      
      if (link == "log") {
        # parameter estimate returned by the glm function
        mod1            <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                               data = casecohort, weights = weights)
        mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                               data = casecohort, weights = weights, 
                               start = coef(mod1), x = TRUE, y = TRUE)
        beta.glm        <- mod.casecohort$coefficients
        
        # Wacholder's bootstrap scheme 
        b.boot <- NULL
        for (b in 1:B) {
          
          casecohort.boot <- NULL
          for (i in seq_along(strata)) {
            casecohort.boot <- rbind(casecohort.boot, casecohort[sample(strata.indiv[[i]], size = strata.n.indiv[i], replace = TRUE),])
          }
          
          mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                             data = casecohort.boot, weights = weights)
          
          trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                                 data = casecohort.boot, weights = weights, 
                                 start = coef(mod1), x = TRUE, y = TRUE), 
                             silent = TRUE)
          
          while(inherits(trial.boot, "try-error")){
            casecohort.boot <- NULL
            for (i in seq_along(strata)) {
              casecohort.boot <- rbind(casecohort.boot, casecohort[sample(strata.indiv[[i]], size = strata.n.indiv[i], replace = TRUE),])
            }
            mod1        <- glm(Y ~ X1 + X2 + X3, family = poisson(link = "log"), 
                               data = casecohort.boot, weights = weights)
            trial.boot  <- try(glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                                   data = casecohort.boot, weights = weights, 
                                   start = coef(mod1), x = TRUE, y = TRUE), 
                               silent = TRUE)
          }
          mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "log"), 
                                 data = casecohort.boot, weights = weights, 
                                 start = coef(mod1), x = TRUE, y = TRUE)
          b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
        }
      }
      
      if (link == "logit") {
        # parameter estimate returned by the glm function
        mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                               data = casecohort, weights = weights, x = TRUE, 
                               y = TRUE)
        beta.glm        <- mod.casecohort$coefficients
        
        # also running unweighted logistic regression
        if(noweight == TRUE){
          mod.casecohort.nw  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                                 data = casecohort, x = TRUE, 
                                 y = TRUE)
          beta.glm.nw        <- mod.casecohort.nw$coefficients
        }
        
        # Wacholder's bootstrap scheme
        b.boot <- NULL
        b.boot.nw <- NULL
        for (b in 1:B) {
          casecohort.boot <- NULL
          for (i in seq_along(strata)) {
            casecohort.boot <- rbind(casecohort.boot, casecohort[sample(strata.indiv[[i]], size = strata.n.indiv[i], replace = TRUE),])
          }
          trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                                   family = binomial(link = "logit"), 
                                   data = casecohort.boot, weights = weights, x = TRUE, 
                                   y = TRUE), silent = TRUE)
          while(inherits(trial.boot, "try-error")){
            casecohort.boot <- NULL
            for (i in seq_along(strata)) {
              casecohort.boot <- rbind(casecohort.boot, casecohort[sample(strata.indiv[[i]], size = strata.n.indiv[i], replace = TRUE),])
            }
            trial.boot    <- try(glm(Y ~ X1 + X2 + X3, 
                                     family = binomial(link = "logit"), 
                                     data = casecohort.boot, weights = weights, x = TRUE, 
                                     y = TRUE), silent = TRUE)
          }
          mod.casecohort  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                                 data = casecohort.boot, weights = weights, x = TRUE, 
                                 y = TRUE)
          b.boot          <- rbind(b.boot, c(mod.casecohort$coefficients))
          if(noweight == TRUE){
            mod.casecohort.nw  <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "logit"), 
                                   data = casecohort.boot, x = TRUE, 
                                   y = TRUE)
            b.boot.nw          <- rbind(b.boot.nw, c(mod.casecohort.nw$coefficients))
          }
        }
      }
  }}
  v.boot          <- apply(b.boot, MARGIN = 2, FUN = var)
  ci.left.boot    <- apply(b.boot, MARGIN = 2, FUN = CI.left.boot)
  ci.right.boot   <- apply(b.boot, MARGIN = 2, FUN = CI.right.boot)
  
  res.casecohort.boot.wald        <- c("casecohort Wacholder boot wald", 
                                       beta.glm, v.boot, 
                                       t(CI(beta.glm, v.boot)))
  res.casecohort.boot.percentile  <- c("casecohort Wacholder boot percentile", 
                                       beta.glm, v.boot, ci.left.boot, 
                                       ci.right.boot)

  if(noweight == TRUE){
  v.boot.nw          <- apply(b.boot.nw, MARGIN = 2, FUN = var)
  ci.left.boot.nw    <- apply(b.boot.nw, MARGIN = 2, FUN = CI.left.boot)
  ci.right.boot.nw   <- apply(b.boot.nw, MARGIN = 2, FUN = CI.right.boot)
  
  res.casecohort.boot.wald.nw        <- c("casecohort Wacholder no weight boot wald", 
                                          beta.glm.nw, v.boot.nw, 
                                          t(CI(beta.glm.nw, v.boot.nw)))
  res.casecohort.boot.percentile.nw  <- c("casecohort Wacholder no weight boot percentile", 
                                          beta.glm.nw, v.boot.nw, ci.left.boot.nw, 
                                          ci.right.boot.nw)
  return(list(res.casecohort.boot.wald = res.casecohort.boot.wald, 
              res.casecohort.boot.percentile = res.casecohort.boot.percentile,
              res.casecohort.boot.wald.nw = res.casecohort.boot.wald.nw,
              res.casecohort.boot.percentile.nw = res.casecohort.boot.percentile.nw))
  }else{
  return(list(res.casecohort.boot.wald = res.casecohort.boot.wald, 
              res.casecohort.boot.percentile = res.casecohort.boot.percentile))
  }
}


## -----------------------------------------------------------------------------
## Function: CI.left.boot()
## -----------------------------------------------------------------------------
## Description: This function computes the left limit of the 100*(1 - alpha)% 
##              bootstrap percentile confidence interval
## -----------------------------------------------------------------------------
## Arguments:
##
## b            vector of bootstrap parameter estimates
##
## alpha        1 minus the desired confidence level. Default value is 0.05
## -----------------------------------------------------------------------------

CI.left.boot <- function(b, alpha = 0.05) {
  quantile(b, probs = c(alpha / 2))
}

## -----------------------------------------------------------------------------
## Function: CI.right.boot()
## -----------------------------------------------------------------------------
## Description: This function computes the right limit of the 100*(1 - alpha)% 
##              bootstrap percentile confidence interval
## -----------------------------------------------------------------------------
## Arguments:
##
## b            vector of bootstrap parameter estimates
##
## alpha        1 minus the desired confidence level. Default value is 0.05
## -----------------------------------------------------------------------------

CI.right.boot <- function(b, alpha = 0.05) {
  quantile(b, probs = c(1 - alpha / 2))
}

