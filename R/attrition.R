

#' Estimate Average Treatment Effect under Differential Attrition 
#' 
#' Estimate attrition scores and estimate ATE under differential attrition 
#'
#' @export 
#' @importFrom CBPS CBPS 
#' @importFrom Formula as.Formula
#' @param formula A formula of the form \code{outcome | treatmet ~  x1 + x2},
#'                   where \code{outcoem} is the outcome of interest with missing value \code{NA}.
#' @param data A data frame.
#' @param estimator An estimator used for estimating ATE. Default is \code{"dr"} (double robust).
attrition <- function(formula, data, estimator = 'dr', est_ps = FALSE) {
  data <- as.data.frame(data)
  
  ## check if \code{outcome | treatment}
  if (length(as.Formula(formula))[1] != 2) {
    stop('Incorrect formula specifications. See the documentation.')
  }
  
  ## check if all variables in the formula exists in the data 
  vars <- all.vars(as.Formula(formula))
  
  if (!all(vars %in% colnames(data))) {
    stop("Variable is missing from the data.")
  }
  
  ## obtain formula 
  fm_outcome <- formula(as.Formula(formula), lhs = 1, rhs = 1)
  fm_treat   <- formula(as.Formula(formula), lhs = 2, rhs = 1)
  fm_X       <- formula(as.Formula(formula), lhs = 0, rhs = 1)
  
  ## create missing outcome 
  var_outcome <- vars[1]
  var_treat   <- vars[2]
  
  ## create attrition indicator (R = 1 if observed) 
  data$R <- R <- as.numeric(!is.na(data[,var_outcome]))
  
  ## estimate attrition score 
  ascore <- attrition_score(update(fm_X, R ~ .), varname_treat = var_treat, data = data)
  
  ## propensity score 
  n1 <- sum(data[,var_treat] == 1); n <- nrow(data)
  
  if (isTRUE(est_ps)) {
    tmp <- capture.output({
      fit_ps <- CBPS(fm_treat, data = data, method = 'over')  
    })
    ps  <- fit_ps$fitted
  } else {
    ps <- rep(n1 / n, n)
  }

  
  ## outcome model 
  Ypred <- outcome_prediction(fm_outcome, varname_treat = var_treat, data = data)
  Yobs  <- pull(data, !!sym(var_outcome))
  
  ## estimate ATE
  if (estimator == 'dr') {
    out <- attrition_dr(var_treat, R, data, Yobs, Ypred, ascore, ps)
  } else if (estimator == 'ipw') {
    out <- attrition_ipw(var_treat, R, data, Yobs, ascore, ps)
  } else {
    stop("Not supported estimator")
  }
  
  return(out)
  
}


#' Estimate Attrition Score 
#' @keywords internal 
#' @importFrom CBPS CBPS
#' @importFrom dplyr pull 
#' @importFrom rlang !! sym
#' @param formula       A formula of the form \code{response ~ x1 + x2}.
#' @param varname_treat A character string of treatment variable name in \code{data}
#' @param data          A data frame.
attrition_score <- function(formula, varname_treat, data) {
  
  ## compute the attrition score 
  tmp <- capture.output(
    fit1 <- CBPS(formula, data = filter(data, !!sym(varname_treat) == 1), method = 'over') 
  )
  tmp <- capture.output(
    fit0 <- CBPS(formula, data = filter(data, !!sym(varname_treat) == 0), method = 'over')
  )
  
  ## append to the data input to ensure the ordering 
  data$attr_score <- 0 
  data$attr_score[data[,varname_treat]==1] <- fit1$fitted
  data$attr_score[data[,varname_treat]==0] <- fit0$fitted
  
  return(data$attr_score)
}

propensity_score <- function(formula, data) {
  fit <- CBPS(formula, data, method = 'over')
  return(fit$fit)
}


#' Outcome prediction 
#' @keywords internal
#' @importFrom dplyr filter
#' @importFrom rlang !! sym
outcome_prediction <- function(formula, varname_treat, data) {
  fit1 <- lm(formula, data = filter(data, !!sym(varname_treat) == 1)) 
  fit0 <- lm(formula, data = filter(data, !!sym(varname_treat) == 0))
  
  Y1pred <- predict(fit1, data)
  Y0pred <- predict(fit0, data)
  
  return(cbind(Y1pred, Y0pred))
}


#' Doubly Robust Estimator 
#' @keywords internal
attrition_dr <- function(var_treat, R, data, Yobs, Ypred, ascore, ps) {
  n <- nrow(data)
  ## estimate E[Y(1)] --------------------------------------
  D <- data[,var_treat]
  use_tr <- R == 1 & D == 1
  # outcome regression 
  EY1 <- Ypred[,1]
  # augmentation 
  EY1[use_tr] <- EY1[use_tr] + 
    (Yobs[use_tr] - Ypred[use_tr, 1]) / (ascore[use_tr] * ps[use_tr])
  
  ## estimate E[Y(0)] -------------------------------------
  use_ct <- R == 1 & D == 0
  # outcome regression 
  EY0 <- Ypred[,2]
  # augmentation 
  EY0[use_ct] <- EY0[use_ct] + 
    (Yobs[use_ct] - Ypred[use_ct, 2]) / (ascore[use_ct] * (1 - ps[use_ct]))
  
  ## estimate ATE -----------------------------------------
  EY1_est <- sum(EY1) / n 
  EY0_est <- sum(EY0) / n 
  ATE_est <- EY1_est - EY0_est
  
  ## estimate variance ------------------------------------
  var_est <- (sum((EY0 - EY0_est)^2) + sum((EY1 - EY1_est)^2)) / n 
  std.err <- sqrt(var_est) / sqrt(n)
  
  ## return object 
  out <- data.frame(
    estimate   = ATE_est, 
    std.error  = std.err, 
    statistic  = ATE_est / std.err, 
    p.value    = dnorm(ATE_est / std.err)
  )
  
  return(out)
}

#' IPW Estimator
#' @keywords internal
attrition_ipw <- function(var_treat, R, data, Yobs, ascore, ps) {
  ## treatment vector -------------------------------------
  D <- data[,var_treat]
  n <- nrow(data)
  
  ## estimate E[Y(1)] -------------------------------------
  use_tr <- D == 1 & R == 1 
  EY1_est <- sum(Yobs[use_tr] / (ascore[use_tr] * ps[use_tr])) / 
                sum(1 / (ascore[use_tr] * ps[use_tr]))
                
  ## estimate E[Y(0)] -------------------------------------
  use_ct <- D == 0 & R == 1
  EY0_est <- sum(Yobs[use_ct] / (ascore[use_ct] * (1 - ps[use_ct]))) / 
                sum(1 / (ascore[use_ct] * (1 - ps[use_ct])))

  ## estimate ATE -----------------------------------------
  ATE_est <- EY1_est - EY0_est
  
  
  ## estimate variance ------------------------------------
  var_est <- (sum(((Yobs[use_tr] - EY1_est) / (ascore[use_tr] * ps[use_tr]))^2) + 
              sum(((Yobs[use_ct] - EY0_est) / (ascore[use_ct] * (1 - ps[use_ct])))^2)) / n 
  std.err <- sqrt(var_est) / sqrt(n)
  
  ## return object 
  out <- data.frame(
    estimate  = ATE_est, 
    std.error = std.err, 
    statistic = ATE_est / std.err, 
    p.value   = dnorm(ATE_est / std.err)
  )
}
