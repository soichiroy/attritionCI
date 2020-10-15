

#' Estimate Average Treatment Effect under Differential Attrition
#'
#' Estimate attrition scores and estimate ATE under differential attrition.
#'
#' @export
#' @importFrom CBPS CBPS
#' @importFrom Formula as.Formula
#' @param formula   A formula of the form \code{outcome | treatmet ~  x1 + x2},
#'                     where \code{outcoem} is the outcome of interest
#'                     with missing value \code{NA}.
#'                     \code{treatment} is the binary treatment indicator.
#'                     Treatment and covariates should not have missing values.
#' @param data      A data frame.
#'                    It should contain all variables specified in the \code{formula}.
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

#' Estimate Bounds for ATE
#' @inheritParams attrition
#' @export
attrition_bound <- function(formula, data, cbps = TRUE, zeta = c(1, 1.1, 1.2)) {
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

  ## -------------------------------------- ##
  ## estimate attrition score
  ## -------------------------------------- ##
  ascore <- attrition_score(update(fm_X, R ~ .), varname_treat = var_treat,
                            data = data, cbps)


  ## -------------------------------------- ##
  ## Compute bounds
  ## -------------------------------------- ##
  Yobs  <- pull(data, !!sym(var_outcome))
  Dtr   <- pull(data, !!sym(var_treat))

  LB <- UB <- rep(NA, length(zeta))
  for (z in 1:length(zeta)) {
    ## compute weights
    w_zeta     <- R * (ascore + (1 - ascore) * zeta[z]) / ascore
    w_zeta_inv <- R * (ascore + (1 - ascore) / zeta[z]) / ascore

    ## upper bound
    Pw0_ub <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = TRUE, weights = w_zeta[Dtr == 0])
    Pw1_ub <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = TRUE, weights = w_zeta_inv[Dtr == 1])

    ## lower bound
    Pw0_lb <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = TRUE, weights = w_zeta_inv[Dtr == 0])
    Pw1_lb <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = TRUE, weights = w_zeta[Dtr == 1])

    ## compute the bound
    int1_ub <- integrate(Pw1_ub, lower = min(Yobs[Dtr == 1], na.rm = TRUE),
                                 upper =  max(Yobs[Dtr == 1], na.rm = TRUE),
                                 stop.on.error = FALSE)
    int0_ub <- integrate(Pw0_ub, lower = min(Yobs[Dtr == 0], na.rm = TRUE),
                                 upper =  max(Yobs[Dtr == 0], na.rm = TRUE),
                                 stop.on.error = FALSE)

    int1_lb <- integrate(Pw1_lb, lower = min(Yobs[Dtr == 1], na.rm = TRUE),
                                 upper =  max(Yobs[Dtr == 1], na.rm = TRUE),
                                 stop.on.error = FALSE)
    int0_lb <- integrate(Pw0_lb, lower = min(Yobs[Dtr == 0], na.rm = TRUE),
                                 upper =  max(Yobs[Dtr == 0], na.rm = TRUE),
                                 stop.on.error = FALSE)

    UB[z] <- int0_ub$value - int1_ub$value
    LB[z] <- int0_lb$value - int1_lb$value

  }


  bounds <- tibble(zeta = zeta, LB = LB, UB = UB)


  ##
  ## Compute Minimum Sensitivity Parameter
  ##
  # zeta_min <- (tau + sqrt(tau^2 + 4 * int_1$value * int_0$value)) / (2 * int_1$value)
  zeta_min <- 1.1


  return(list(bouns = bounds, MSP = zeta_min))
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




#' IPW Estimator
#' @keywords internal
calc_psi <- function(var_treat, R, data, Yobs, ascore, zeta) {
  ## treatment vector -------------------------------------
  D <- data[,var_treat]
  n <- nrow(data)

  ## estimate 𝞧(1,1) and 𝞧(1,0) -------------------------
  use_tr <- D == 1 & R == 1
  n1 <- sum(D == 1)
  p11 <- sum(Yobs[use_tr] / ascore[use_tr]) / n
  p10 <- sum(Yobs[use_tr] * (1 - ascore[use_tr]) / ascore[use_tr]) / n

  ## estimate 𝞧(0,1) and 𝞧(0,0) --------------------------
  use_ct <- D == 0 & R == 1
  n0 <- sum(D == 0)
  p01 <- sum(Yobs[use_ct] / ascore[use_ct]) / n
  p00 <- sum(Yobs[use_ct] * (1 - ascore[use_ct]) / ascore[use_ct]) / n

  ## estimate bounds for ATE ------------------------------
  ATE_lb <- ATE_ub <- rep(NA, length(zeta))
  for (z in seq_along(zeta)) {
    if (zeta[z] >= 1) {
      ATE_lb[z] <- (p11 + p10 / zeta[z]) - (p01 + zeta[z] * p00)
      ATE_ub[z] <- (p11 + zeta[z] * p10) - (p01 + p00 / zeta[z])
    } else {
      zeta_inv <- 1 / zeta[z]
      ATE_lb[z] <- (p11 + p10 / zeta_inv) - (p01 + zeta_inv * p00)
      ATE_ub[z] <- (p11 + zeta_inv * p10) - (p01 + p00 / zeta_inv)
    }
  }


  ## estimate variance ------------------------------------

  ## return object
  out <- list(
    psi = c(p11, p10, p01, p00),
    lb  = ATE_lb,
    ub  = ATE_ub
  )
  return(out)
}
