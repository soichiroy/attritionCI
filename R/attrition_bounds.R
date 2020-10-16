
#' Estimate Bounds on Causal Effect
#'
#' Compute the bounds on casual effects with sensitivity parameter zeta.
#'
#' @inheritParams attrition
#' @param zeta A vector of sensitivity parameters. The values should be greater than or equal to 1.
#' @param probs A vector of quantiles. The value should be between 0 and 1.
#'              This argument is only meaninful when the quantile treatment effeect (`qoi = "qte"`)
#'              is computed.
#' @param qoi Quantity of interest. Currently supports the average treatment efeect (`"ate"`) or
#'             the quantile treatment effect (`"qte"`)
#' @importFrom dplyr left_join mutate pull rename select %>%
#' @importFrom rlang !! sym
#' @import Formula
#' @export
attrition_bound <- function(
  formula, data, qoi = "ate", cbps = TRUE,
  zeta = c(1, 1.1, 1.2), probs = c(0.25, 0.5, 0.75),
  options = list(ci = FALSE, n_boot = 500)
) {

  ## treat dataset as data.frame to avoid errors
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
  ascore <- attrition_score(update(fm_X, R ~ .),
                            varname_treat = var_treat,
                            data = data, cbps)

  ## -------------------------------------- ##
  ## Get outcome and treatment
  ## -------------------------------------- ##
  Yobs  <- pull(data, !!sym(var_outcome))  ## should have na
  Dtr   <- pull(data, !!sym(var_treat))

  ## -------------------------------------- ##
  ## Compute bounds
  ## -------------------------------------- ##
  if (qoi == "ate") {
    ## compute the identification bounds on ATE τ
    bounds <- attrition_bound_ate(Yobs, Dtr, R, ascore, zeta)

    if (isTRUE(options$ci)) {
      ## bootstrap to compute the variance of τ
      bounds_ci <- attrition_bound_ate_boot(Yobs, Dtr, R, ascore, zeta, options$n_boot)

      ## compute the CI
      # 90%
      bounds$CI90_LB <- bounds$LB + qnorm(0.1) * bounds_ci$se_lb
      bounds$CI90_UB <- bounds$UB + qnorm(0.9) * bounds_ci$se_ub

      # 95%
      bounds$CI95_LB <- bounds$LB + qnorm(0.05) * bounds_ci$se_lb
      bounds$CI95_UB <- bounds$UB + qnorm(0.95) * bounds_ci$se_ub
    }

  } else if (qoi == "qte") {
    ## compute identification bounds on QTE η(α) for α ∈ probs
    bounds <- attrition_bound_qte(Yobs, Dtr, R, ascore, zeta, probs)

    if (isTRUE(options$ci)) {
      bounds_ci <- attrition_bound_qte_boot(Yobs, Dtr, R, ascore, zeta, probs, options$n_boot)

      bounds <- left_join(bounds, bounds_ci, by = c("zeta", "probs")) %>%
                mutate(
                  CI90_LB = LB + qnorm(0.10) * se_lb,
                  CI90_UB = UB + qnorm(0.90) * se_ub,
                  CI95_LB = LB + qnorm(0.05) * se_lb,
                  CI95_UB = UB + qnorm(0.95) * se_ub) %>%
                select(-se_lb, -se_ub)
    }


  } else {
    stop("Supported QOI is either ATE or QTE.")
  }


  return(bounds)
}


#' Compute bounds on ATE
#' @inheritParams attrition
#' @importFrom spatstat ewcdf
#' @importFrom dplyr tibble
#' @keywords internal
attrition_bound_ate <- function(Yobs, Dtr, R, ascore, zeta) {
  LB <- UB <- rep(NA, length(zeta))
  for (z in 1:length(zeta)) {
    ## compute weights
    w_zeta     <- R * (ascore + (1 - ascore) * zeta[z]) / ascore
    w_zeta_inv <- R * (ascore + (1 - ascore) / zeta[z]) / ascore

    ## upper bound
    Pw0_ub <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
                              weights = w_zeta[Dtr == 0] / sum(Dtr == 0))
    Pw1_ub <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
                              weights = w_zeta_inv[Dtr == 1] / sum(Dtr == 1))

    ## lower bound
    Pw0_lb <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
                              weights = w_zeta_inv[Dtr == 0] / sum(Dtr == 0))
    Pw1_lb <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
                              weights = w_zeta[Dtr == 1] / sum(Dtr == 1))

    ## compute the bound
    int1_ub <- integrate(Pw1_ub, lower = min(Yobs, na.rm = TRUE),
                                 upper = max(Yobs, na.rm = TRUE),
                                 stop.on.error = FALSE)
    int0_ub <- integrate(Pw0_ub, lower = min(Yobs, na.rm = TRUE),
                                 upper = max(Yobs, na.rm = TRUE),
                                 stop.on.error = FALSE)

    int1_lb <- integrate(Pw1_lb, lower = min(Yobs, na.rm = TRUE),
                                 upper = max(Yobs, na.rm = TRUE),
                                 stop.on.error = FALSE)
    int0_lb <- integrate(Pw0_lb, lower = min(Yobs, na.rm = TRUE),
                                 upper = max(Yobs, na.rm = TRUE),
                                 stop.on.error = FALSE)

    UB[z] <- int0_ub$value - int1_ub$value
    LB[z] <- int0_lb$value - int1_lb$value
  }


  bounds <- tibble(zeta = zeta, LB = LB, UB = UB)
  return(bounds)
}


attrition_bound_ate_boot <- function(Yobs, Dtr, R, ascore, zeta, n_boot = 500) {

  n <- length(Dtr)

  LB <- UB <- matrix(NA, nrow = n_boot, ncol = length(zeta))
  for (i in 1:n_boot) {
    ## re-sampled weights
    w_boot <- as.vector(rmultinom(1, n, prob = rep(1/n, n)))

    for (z in 1:length(zeta)) {
      ## compute weights
      w_zeta     <- R * (ascore + (1 - ascore) * zeta[z]) / ascore
      w_zeta_inv <- R * (ascore + (1 - ascore) / zeta[z]) / ascore

      ## upper bound
      Pw0_ub <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
        weights = w_zeta[Dtr == 0] * w_boot[Dtr == 0] / sum(w_boot[Dtr == 0]))
      Pw1_ub <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
        weights = w_zeta_inv[Dtr == 1] * w_boot[Dtr == 1] / sum(w_boot[Dtr == 1]))

      ## lower bound
      Pw0_lb <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
        weights = w_zeta_inv[Dtr == 0] * w_boot[Dtr == 0] / sum(w_boot[Dtr == 0]))
      Pw1_lb <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
        weights = w_zeta[Dtr == 1] * w_boot[Dtr == 1] / sum(w_boot[Dtr == 1]))

      ## compute the bound
      int1_ub <- integrate(Pw1_ub, lower = min(Yobs, na.rm = TRUE),
                                   upper = max(Yobs, na.rm = TRUE),
                                   stop.on.error = FALSE)
      int0_ub <- integrate(Pw0_ub, lower = min(Yobs, na.rm = TRUE),
                                   upper = max(Yobs, na.rm = TRUE),
                                   stop.on.error = FALSE)

      int1_lb <- integrate(Pw1_lb, lower = min(Yobs, na.rm = TRUE),
                                   upper = max(Yobs, na.rm = TRUE),
                                   stop.on.error = FALSE)
      int0_lb <- integrate(Pw0_lb, lower = min(Yobs, na.rm = TRUE),
                                   upper = max(Yobs, na.rm = TRUE),
                                   stop.on.error = FALSE)

      UB[i, z] <- int0_ub$value - int1_ub$value
      LB[i, z] <- int0_lb$value - int1_lb$value
    }
  }

  ### compute the variance
  sigma_ub <- sqrt(apply(UB, 2, var))
  sigma_lb <- sqrt(apply(LB, 2, var))

  return(list(se_ub = sigma_ub, se_lb = sigma_lb))

}



#' Estimate Bounds for QTE
#' @inheritParams attrition
#' @importFrom spatstat ewcdf
#' @importFrom dplyr bind_rows tibble
#' @keywords internal
attrition_bound_qte <- function(Yobs, Dtr, R, ascore, zeta, probs) {

  ## -------------------------------------- ##
  ## Compute bounds
  ## -------------------------------------- ##
  res <- vector('list', length = length(zeta))
  for (z in 1:length(zeta)) {
    UB <- LB <- rep(NA, length(probs))

    ## compute weights
    w_zeta     <- R * (ascore + (1 - ascore) * zeta[z]) / ascore
    w_zeta_inv <- R * (ascore + (1 - ascore) / zeta[z]) / ascore

    ## Cumulative function: F_w(ζ)(y)
    Pw0_zeta <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
                              weights = w_zeta[Dtr == 0] / sum(Dtr == 0))
    Pw0_zinv <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
                              weights = w_zeta_inv[Dtr == 0] / sum(Dtr == 0))
    Pw1_zinv <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
                              weights = w_zeta_inv[Dtr == 1] / sum(Dtr == 1))
    Pw1_zeta <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
                              weights = w_zeta[Dtr == 1] / sum(Dtr == 1))

    ## compute F^{-1}(α)
    y_range <- c(min(Yobs, na.rm = TRUE), max(Yobs, na.rm = TRUE))
    for (alpha in 1:length(probs)) {
      ## compute LB
      Pw1_zeta_a <- function(x) { Pw1_zeta(x) - probs[alpha] }
      Pw0_zinv_a <- function(x) { Pw0_zinv(x) - probs[alpha] }
      LB[alpha]  <- uniroot(Pw1_zeta_a, interval = y_range)$root -
                    uniroot(Pw0_zinv_a, interval = y_range)$root

      ## compute UB
      Pw1_zinv_a <- function(x) { Pw1_zinv(x) - probs[alpha] }
      Pw0_zeta_a <- function(x) { Pw0_zeta(x) - probs[alpha] }
      UB[alpha]  <- uniroot(Pw1_zinv_a, interval = y_range)$root -
                    uniroot(Pw0_zeta_a, interval = y_range)$root
    }

    res[[z]] <- tibble(zeta = zeta[z], probs = probs, LB = LB, UB = UB)
  }


  bounds <- bind_rows(res)
  return(bounds)
}




#' Estimate Variance of Bounds for QTE with Bootstrap
#' @inheritParams attrition
#' @importFrom spatstat ewcdf
#' @importFrom dplyr across bind_rows group_by summarise tibble rename %>%
#' @keywords internal
attrition_bound_qte_boot <- function(Yobs, Dtr, R, ascore, zeta, probs, n_boot) {
  n0 <- sum(Dtr == 0)
  n1 <- sum(Dtr == 1)

  boot_list <- vector('list', length = n_boot)
  for (i in 1:n_boot) {
    ## re-sampled weights
    w_boot <- as.vector(rmultinom(1, n, prob = rep(1/n, n)))

    ## -------------------------------------- ##
    ## Compute bounds
    ## -------------------------------------- ##
    res <- vector('list', length = length(zeta))
    for (z in 1:length(zeta)) {
      UB <- LB <- rep(NA, length(probs))

      ## compute weights
      w_zeta     <- R * (ascore + (1 - ascore) * zeta[z]) / ascore
      w_zeta_inv <- R * (ascore + (1 - ascore) / zeta[z]) / ascore

      ## Cumulative function: F_w(ζ)(y)
      Pw0_zeta <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
        weights = w_zeta[Dtr == 0] * w_boot[Dtr == 0] / sum(w_boot[Dtr == 0]))
      Pw0_zinv <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
        weights = w_zeta_inv[Dtr == 0] * w_boot[Dtr == 0] / sum(w_boot[Dtr == 0]))
      Pw1_zinv <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
        weights = w_zeta_inv[Dtr == 1] * w_boot[Dtr == 1] / sum(w_boot[Dtr == 1]))
      Pw1_zeta <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
        weights = w_zeta[Dtr == 1] * w_boot[Dtr == 1] / sum(w_boot[Dtr == 1]))

      ## compute F^{-1}(α)
      y_range <- c(min(Yobs, na.rm = TRUE), max(Yobs, na.rm = TRUE))
      for (alpha in 1:length(probs)) {
        ## compute UB
        Pw1_zeta_a <- function(x) { Pw1_zeta(x) - probs[alpha] }
        Pw0_zinv_a <- function(x) { Pw0_zinv(x) - probs[alpha] }
        UB[alpha]  <- uniroot(Pw1_zeta_a, interval = y_range)$root -
                      uniroot(Pw0_zinv_a, interval = y_range)$root

        ## compute LB
        Pw1_zinv_a <- function(x) { Pw1_zinv(x) - probs[alpha] }
        Pw0_zeta_a <- function(x) { Pw0_zeta(x) - probs[alpha] }
        LB[alpha]  <- uniroot(Pw1_zinv_a, interval = y_range)$root -
                      uniroot(Pw0_zeta_a, interval = y_range)$root
      }

      res[[z]] <- tibble(zeta = zeta[z], probs = probs, LB = LB, UB = UB, boot_id = i)
    }

    boot_list[[i]] <- res
  }


  ## compute the variance
  bounds <- bind_rows(boot_list) %>%
    group_by(zeta, probs) %>%
    summarise(across(c(LB, UB), sd), .groups = "drop") %>%
    rename(se_lb = LB, se_ub = UB)

  return(bounds)
}
