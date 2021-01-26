

#' Estimate Average Treatment Effect under Differential Attrition
#'
#' Estimate attrition scores and estimate ATE under differential attrition.
#'
#' @export
#' @importFrom CBPS CBPS
#' @import Formula
#' @param formula   A formula of the form \code{outcome | treatmet ~  x1 + x2},
#'                  where \code{outcoem} is the outcome of interest
#'                  with missing value \code{NA}.
#'                  \code{treatment} is the binary treatment indicator.
#'                  Treatment and covariates should not have missing values.
#' @param data      A data frame.
#'                  It should contain all variables specified in the \code{formula}.
#' @param estimator An estimator used for estimating ATE.
#'                  Default is \code{"dr"} (double robust), and the other option is
#'                  \code{"ipw"} (inverse probability weighting).
#' @param est_ps    A boolean argument. If set \code{TRUE}, the propensity score is
#'                  estimated from the data. Otherwise, the complete randomization is
#'                  assumed, and the propensity score is set to the empirical propotion of
#'                  the treated units.
#'                  This option should be \code{TRUE} when the data is from
#'                  an observational study.
#' @param cbps      A boolean argument. If \code{TRUE}, the covariate balancing propensity
#'                  score (CBPS) is used to estimate the attrition score (and the propensity
#'                  score when \code{est_ps = TRUE}). Otherwise, \code{glm} is used.
#' @return A data.frame that contains treatment effect estimates as well as standard errors.
attrition <- function(formula, data, estimator = 'dr', est_ps = FALSE, cbps = TRUE) {
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
  ascore <- attrition_score(update(fm_X, R ~ .), varname_treat = var_treat,
                            data = data, cbps)

  ## propensity score
  ps <-  propensity_score(fm_treat, data, est_ps, cbps)

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
#' @param cbps          A boolean argument if CBPS is used.
attrition_score <- function(formula, varname_treat, data, cbps) {

  ## compute the attrition score
  if (isTRUE(cbps)) {
    tmp <- capture.output(
      fit1 <- CBPS(formula, data = filter(data, !!sym(varname_treat) == 1), method = 'over')
    )
    tmp <- capture.output(
      fit0 <- CBPS(formula, data = filter(data, !!sym(varname_treat) == 0), method = 'over')
    )

    fit1_fitted <- fit1$fitted
    fit0_fitted <- fit0$fitted
  } else {
    fit1 <- glm(formula, data = filter(data, !!sym(varname_treat) == 1), family = binomial())
    fit0 <- glm(formula, data = filter(data, !!sym(varname_treat) == 0), family = binomial())
    fit1_fitted <- predict(fit1, type = 'response')
    fit0_fitted <- predict(fit0, type = 'response')
  }

  ## append to the data input to ensure the ordering
  data$attr_score <- 0
  data$attr_score[data[,varname_treat]==1] <- fit1_fitted
  data$attr_score[data[,varname_treat]==0] <- fit0_fitted

  return(data$attr_score)
}


#' Estimate propensity score
#' @keywords internal
#' @importFrom CBPS CBPS
#' @param formula A formula \code{treatment ~ x1 + x2}.
#' @param data    A data frame.
#' @param est_ps  A boolean argument if the propensity score is estimated.
#'                If set \code{FALSE}, empirical fraction of treated units is used.
#' @return A vector of propensity score with length n (the number of observations).
propensity_score <- function(formula, data, est_ps, cbps) {
  if (isTRUE(est_ps)) {
    if (isTRUE(cbps)) {
      ## estimate ps with CBPS
      tmp <- capture.output({
        fit <- CBPS(formula, data, method = 'over')
      })
      ps  <- fit$fitted
    } else {
      ## estimate ps with logit
      fit <- glm(formula, data, family = binomial())
      ps  <- predict(fit, type = 'response')
    }
  } else {
    var_treat <- all.vars(formula)[1]
    n1 <- sum(data[,var_treat] == 1)
    n  <- nrow(data)
    ps <- rep(n1 / n, n)
  }

  return(ps)
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
