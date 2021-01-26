
#' Estimate Bounds on Causal Effect
#'
#' Compute the bounds on casual effects with sensitivity parameter zeta.
#'
#' @inheritParams attrition
#' @param zeta    A vector of sensitivity parameters. The values should be between 1 and 2.
#' @param probs   A vector of quantiles. The value should be between 0 and 1.
#'                This argument is meaninful when the quantile treatment effeect (`qoi = "qte"`)
#'                is specified as the quantity of interest.
#' @param qoi     Quantity of interest.
#'                The function currently supports the average treatment efeect (`"ate"`) or
#'                the quantile treatment effect (`"qte"`).
#' @param options A list of the following optional arguments:
#'                \describe{
#'                  \item{ci}{A boolean argument.
#'                    If `TRUE`, the confidence intervals are computed via bootstrap.}
#'                  \item{n_boot}{The number of bootstrap iterations.}
#'                }
#' @importFrom dplyr left_join mutate pull rename select %>%
#' @importFrom rlang !! sym
#' @import Formula
#' @export
attrition_bound <- function(
  formula, data, qoi = "ate", cbps = TRUE,
  zeta = c(1, 1.1, 1.2), probs = c(0.25, 0.5, 0.75),
  ascore  = NULL,
  options = list(ci = TRUE, n_boot = 100)
) {

  ## treat dataset as data.frame to avoid errors
  data <- as.data.frame(data)

  ## -------------------------------------- ##
  ## Input checks
  ## -------------------------------------- ##
  ## check if \code{outcome | treatment}
  if (length(as.Formula(formula))[1] != 2) {
    stop('Incorrect formula specifications. See the documentation.')
  }
  if (any(zeta > 2)) {
    stop("Value of zeta is too large. Specify zeta < 2.")
  }
  if (any(zeta < 1)) {
    stop("Invalid value of zeta. Specify zeta >= 1.")
  }
  if (any(probs > 1) | any(probs < 0)) {
    stop("Invalid value of probs. Specify 0 <= probs <= 1.")
  }

  ## check if all variables in the formula exists in the data
  vars <- all.vars(as.Formula(formula))
  if (!all(vars %in% colnames(data))) {
    stop("Variable in the formula is missing from the data.")
  }

  ## option
  if (!exists('ci', options)) options$ci <- TRUE
  if (!exists('n_boot', options)) options$n_boot <- 100

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
  if (is.null(ascore)) {
    ascore <- attrition_score(
      update(fm_X, R ~ .), varname_treat = var_treat, data = data, cbps
    )
  }

  ## -------------------------------------- ##
  ## Get outcome and treatment
  ## -------------------------------------- ##
  Yobs  <- pull(data, !!sym(var_outcome))  ## should have NA
  Dtr   <- pull(data, !!sym(var_treat))
  if (sum(is.na(Yobs)) == 0) { stop("No NA found in the outcome variable.") }
  is_discrete <- ifelse( length(unique(na.omit(Yobs))) <= 10, TRUE, FALSE )
  if (qoi == 'qte' & isTRUE(is_discrete)) {
    stop("QTE option is not supproted for the discrete outcome with small number of support points")
  }

  ## -------------------------------------- ##
  ## Compute bounds
  ## -------------------------------------- ##
  if (qoi == "ate") {
    ## compute the identification bounds on ATE τ
    bounds <- attrition_bound_ate(Yobs, Dtr, R, ascore, zeta, is_discrete)

    if (isTRUE(options$ci)) {
      ## bootstrap to compute the variance of τ
      # bounds_ci <- attrition_bound_ate_boot(Yobs, Dtr, R, ascore, zeta, options$n_boot)
      bounds_ci <- attrition_bound_ate_boot_full(
        update(fm_X, R ~ .), varname_treat = var_treat, Yobs, R, data = data, cbps,
        zeta, options$n_boot, is_discrete
      )
      cv_alpha <- optimize_alpha(alpha = 0.1, bounds$UB, bounds$LB,
                              bounds_ci$se_lb, bounds_ci$se_ub, n = length(Yobs))
      ## compute the CI
      # 90%
      bounds$CI90_LB <- bounds$LB - cv_alpha * bounds_ci$se_lb
      bounds$CI90_UB <- bounds$UB + cv_alpha * bounds_ci$se_ub

      bounds$SE_UB <- bounds_ci$se_ub
      bounds$SE_LB <- bounds_ci$se_lb
    }
  } else if (qoi == "qte") {
    ## compute identification bounds on QTE η(α) for α ∈ probs
    bounds <- attrition_bound_qte(Yobs, Dtr, R, ascore, zeta, probs)

    if (isTRUE(options$ci)) {
      ## compute CI
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


optimize_alpha <- function(alpha = 0.1, UB, LB, var1, var2, n) {
  a_range <- c(1, 2)
  alpha_adaptive <- rep(NA, length(UB))
  for (i in 1:length(UB)) {
    fn <- function(x) {
      pnorm(x + sqrt(n) * (UB[i] - LB[i]) / max(var1[i], var2[i])) - pnorm(-x) - (1 - alpha)
    }
    alpha_adaptive[i] <- uniroot(fn, interval = a_range)$root
  }

  return(alpha_adaptive)
}

#' Compute bounds on ATE
#' @inheritParams attrition
#' @param is_discrete A boolean argument if the outcome is a discrete variable or not.
#' @importFrom spatstat ewcdf
#' @importFrom dplyr tibble
#' @keywords internal
attrition_bound_ate <- function(Yobs, Dtr, R, ascore, zeta, is_discrete) {
  LB <- UB <- UB2 <- LB2 <- rep(NA, length(zeta))
  for (z in 1:length(zeta)) {
    ## compute weights
    w_zeta     <- R * (ascore + (1 - ascore) * zeta[z]) / ascore
    w_zeta_inv <- R * (ascore + (1 - ascore) / zeta[z]) / ascore

    ## data summary
    n1 <- sum(Dtr == 1); n0 <- sum(Dtr == 0)
    Y1 <- Yobs[Dtr == 1]; Y0 <- Yobs[Dtr == 0]

    ## upper bound
    Pw0_ub <- spatstat::ewcdf(Y0, normalise = FALSE,
                              weights = w_zeta[Dtr == 0] / n0)
    Pw1_ub <- spatstat::ewcdf(Y1, normalise = FALSE,
                              weights = w_zeta_inv[Dtr == 1] / n1)

    ## lower bound
    Pw0_lb <- spatstat::ewcdf(Y0, normalise = FALSE,
                              weights = w_zeta_inv[Dtr == 0] / n0)
    Pw1_lb <- spatstat::ewcdf(Y1, normalise = FALSE,
                              weights = w_zeta[Dtr == 1] / n1)

    if (isTRUE(is_discrete)) {
      ## should check if it ranges from 0 -----
      y_range <- sort(unique(na.omit(Yobs)))
      UB[z] <- LB[z] <- 0
      for (i in 1:length(y_range)) {
        UB[z] <- UB[z] + Pw0_ub(y_range[i]) - Pw1_ub(y_range[i])
        LB[z] <- LB[z] + Pw0_lb(y_range[i]) - Pw1_lb(y_range[i])
      }
    } else {
      ## compute the bound
      y_max <- max(Yobs, na.rm = TRUE)
      y_min <- min(Yobs, na.rm = TRUE)
      int1_ub <- integrate(Pw1_ub, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)
      int0_ub <- integrate(Pw0_ub, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)

      int1_lb <- integrate(Pw1_lb, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)
      int0_lb <- integrate(Pw0_lb, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)

      UB[z] <- int0_ub$value - int1_ub$value
      LB[z] <- int0_lb$value - int1_lb$value


      ## alternative method
      UB2[z] <- sum(Y1 * w_zeta_inv[Dtr == 1], na.rm = TRUE) / n1 -
                sum(Y0 * w_zeta[Dtr == 0], na.rm = TRUE) / n0
      LB2[z] <- sum(Y1 * w_zeta[Dtr == 1], na.rm = TRUE) / n1 -
                sum(Y0 * w_zeta_inv[Dtr == 0], na.rm = TRUE) / n0
    }

  }


  bounds <- tibble(zeta = zeta, LB = LB, UB = UB, UB2 = UB2, LB2 = LB2)
  return(bounds)
}


#' Compute Varinace of ATE Bounds via Full Bootstrap
#' @importFrom dplyr pull
#' @importFrom rlang !! sym
attrition_bound_ate_boot_full <- function(
  formula, varname_treat, Y, R, data, cbps, zeta, n_boot, is_discrete
) {

  n_obs <- nrow(data)
  UB <- LB <- matrix(NA, nrow = n_boot, ncol = length(zeta))

  for (i in 1:n_boot) {
    ## resample data
    id_boot <- sample(1:n_obs, size = n_obs, replace = TRUE)

    Yboot <- Y[id_boot]
    Rboot <- R[id_boot]
    dat_boot <- data[id_boot, ]
    Dtr   <- pull(dat_boot, !!sym(varname_treat))

    ## reestimate  ascore
    ascore_boot <- attrition_score(formula, varname_treat,
                                   data = dat_boot, cbps)

    ## estimate bounds
    bounds <- attrition_bound_ate(Yboot, Dtr, Rboot, ascore_boot, zeta, is_discrete)
    UB[i,] <- bounds$UB
    LB[i,] <- bounds$LB
  }

  ## compute the variance
  se_ub <- apply(UB, 2, sd)
  se_lb <- apply(LB, 2, sd)
  return(list(se_ub = se_ub, se_lb = se_lb))
}

#' Compute Variance of ATE Bound via Bootstrap
#' @keywords internal
#' @param Yobs A vector of outcomes. It should contain at least two `NA` values for each condition.
#' @param Dtr A vector of binary treatment indicator.
#' @param R A vector of missing value indicator.
#' @param ascore A vector of estimated attrition score.
#' @param zeta A vector of the sensitivity parameter. The value should be larger than 1.
attrition_bound_ate_boot <- function(Yobs, Dtr, R, ascore, zeta, n_boot = 500) {

  n <- length(Dtr)
  LB <- UB <- matrix(NA, nrow = n_boot, ncol = length(zeta))
  LB2 <- UB2 <- matrix(NA, nrow = n_boot, ncol = length(zeta))
  Y1 <- Yobs[Dtr == 1]; Y0 <- Yobs[Dtr == 0]

  y_min  <- min(Yobs, na.rm = TRUE)
  y_max  <- max(Yobs, na.rm = TRUE)

  ## bootstrap -----------------------------------------------------------
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
      int1_ub <- integrate(Pw1_ub, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)
      int0_ub <- integrate(Pw0_ub, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)

      int1_lb <- integrate(Pw1_lb, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)
      int0_lb <- integrate(Pw0_lb, lower = y_min, upper = y_max,
                                   stop.on.error = FALSE)

      UB[i, z] <- int0_ub$value - int1_ub$value
      LB[i, z] <- int0_lb$value - int1_lb$value

      ## alternative method
      LB2[i, z] <- sum(Y1 * w_boot[Dtr == 1] * w_zeta_inv[Dtr == 1], na.rm = TRUE) /
                  sum(w_boot[Dtr == 1]) -
                sum(Y0 * w_boot[Dtr == 0] * w_zeta[Dtr == 0], na.rm = TRUE) /
                  sum(w_boot[Dtr == 0])
      UB2[i, z] <- sum(Y1 * w_boot[Dtr == 1] * w_zeta[Dtr == 1], na.rm = TRUE) /
                  sum(w_boot[Dtr == 1]) -
                sum(Y0 * w_boot[Dtr == 0] * w_zeta_inv[Dtr == 0], na.rm = TRUE) /
                  sum(w_boot[Dtr == 0])

    }
  }
  ## end of bootstrap iterations -----------------------------------------

  ### compute the variance
  sigma_ub <- apply(UB, 2, sd)
  sigma_lb <- apply(LB, 2, sd)

  sigma_ub2 <- apply(UB2, 2, sd)
  sigma_lb2 <- apply(LB2, 2, sd)


  return(list(se_ub = sigma_ub, se_lb = sigma_lb,
              se_ub2 = sigma_ub2, se_lb2 = sigma_lb2))

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
attrition_bound_qte_boot <- function(Yobs, Dtr, R, ascore, zeta, probs, n_boot, verbose = FALSE) {

  ## check verbose
  if (isTRUE(verbose)) {
    iter_show <- round(n_boot * 0.1)
  }

  n <- length(Dtr)
  boot_list <- vector('list', length = n_boot)


  ## bootstrap -----------------------------------------------------------
  # reject the bootstrap replica when inversion fails
  i <- 1
  while(i <= n_boot) {
    tryCatch({
      # for (i in 1:n_boot) {
        ## re-sampled weights
        w_boot <- as.vector(rmultinom(1, n, prob = rep(1/n, n)))

        ## -------------------------------------- ##
        ## Compute bounds
        ## -------------------------------------- ##
        res <- vector('list', length = length(zeta))
        for (z in 1:length(zeta)) {
          UB <- LB <- rep(NA, length(probs))

          ## compute weights --------------------------------------------------
          w_zeta     <- R * (ascore + (1 - ascore) * zeta[z]) / ascore
          w_zeta_inv <- R * (ascore + (1 - ascore) / zeta[z]) / ascore

          ## Cumulative function: F_w(ζ)(y) -----------------------------------
          Pw0_zeta <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
            weights = w_zeta[Dtr == 0] * w_boot[Dtr == 0] / sum(w_boot[Dtr == 0]))
          Pw0_zinv <- spatstat::ewcdf(Yobs[Dtr == 0], normalise = FALSE,
            weights = w_zeta_inv[Dtr == 0] * w_boot[Dtr == 0] / sum(w_boot[Dtr == 0]))
          Pw1_zinv <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
            weights = w_zeta_inv[Dtr == 1] * w_boot[Dtr == 1] / sum(w_boot[Dtr == 1]))
          Pw1_zeta <- spatstat::ewcdf(Yobs[Dtr == 1], normalise = FALSE,
            weights = w_zeta[Dtr == 1] * w_boot[Dtr == 1] / sum(w_boot[Dtr == 1]))

          ## compute F^{-1}(α) via root finding -------------------------------
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

        # save object
        boot_list[[i]] <- res

        # update iterator
        i <- i + 1
      }, error = function(e) {
      NULL
    })

    ## verbose ----------------------------------------------------
    if (isTRUE(verbose)) {
      if ((i %% iter_show) == 0) {
          cat('\r', i, "out of", n_boot, "bootstrap iterations")
          flush.console()
      }
    }
    ## end of verbose ---------------------------------------------

  }
  ## end of bootstrap iterations -----------------------------------------



  ## compute the variance
  bounds <- bind_rows(boot_list) %>%
    group_by(zeta, probs) %>%
    summarise(across(c(LB, UB), sd), .groups = "drop") %>%
    rename(se_lb = LB, se_ub = UB)

  return(bounds)
}
