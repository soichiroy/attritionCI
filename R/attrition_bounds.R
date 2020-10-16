
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
#' @export
attrition_bound <- function(
  formula, data, qoi = "ate", cbps = TRUE, zeta = c(1, 1.1, 1.2), probs = c(0.25, 0.5, 0.75)
) {
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
  ## Get outcome and treatment
  ## -------------------------------------- ##
  Yobs  <- pull(data, !!sym(var_outcome))  ## should have na
  Dtr   <- pull(data, !!sym(var_treat))

  ## -------------------------------------- ##
  ## Compute bounds
  ## -------------------------------------- ##
  if (qoi == "ate") {
    bounds <- attrition_bound_ate(Yobs, Dtr, R, ascore, zeta)
  } else if (qoi == "qte") {
    bounds <- attrition_bound_qte(Yobs, Dtr, R, ascore, zeta, probs)
  } else {
    stop("Supported QOI is either ATE or QTE.")
  }


  return(bounds)
}


#' Compute bounds on ATE
#' @inheritParams attrition
#' @importFrom spatstat ewcdf
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


}





#' Estimate Bounds for QTE
#' @inheritParams attrition
#' @importFrom spatstat ewcdf
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

    res[[z]] <- tibble(zeta = zeta[z], probs = probs, LB = LB, UB = UB)
  }


  bounds <- bind_rows(res)
  return(bounds)
}
