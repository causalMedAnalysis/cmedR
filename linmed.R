#' Product-of-coefficients estimator for natural effects
#' 
#' @description
#' `linmed()` uses the product-of-coefficients estimator, based on linear 
#' models, to estimate the total effect (ATE), natural direct effect (NDE), 
#' natural indirect effect (NIE), and controlled direct effect (CDE). The 
#' function supports estimation of both univariate and multivariate natural 
#' effects.
#' 
#' @details
#' TEMPORARY PLACEHOLDER
#' 
#' @param data A data frame.
#' @param D A character scalar identifying the name of the exposure variable in 
#'   `data`. `D` is a character string, but the exposure variable it identifies 
#'   must be numeric.
#' @param M A character vector (of one or more elements) identifying the names 
#'   of the mediator variables in `data`. If you are estimating univariate 
#'   natural effects (with a single mediator), `M` should be a character scalar 
#'   (a vector with only one element)---e.g., `M = "ever_unemp_age3539"`. If you 
#'   are estimating multivariate natural effects (with multiple mediators), `M` 
#'   should be a character vector identifying all of the mediators---e.g., 
#'   `M = c("ever_unemp_age3539", "log_faminc_adj_age3539")`. Also note that `M` 
#'   is a character vector, but the mediator variable(s) it identifies must each 
#'   be numeric.
#' @param Y A character scalar identifying the name of the outcome variable in 
#'   `data`. `Y` is a character string, but the outcome variable it identifies 
#'   must be numeric.
#' @param C A character vector (of one or more elements) identifying the names 
#'   of the covariate variables in `data` that you wish to include in both the 
#'   mediator and outcome models. If there are no such covariates you wish to 
#'   include, leave `C` as its default null argument.
#' @param d,dstar A pair of arguments, each a numeric scalar denoting a specific 
#'   value of the exposure `D`. The exposure contrast of interest is 
#'   `d - dstar`.
#' @param m A numeric vector (of one or more elements) denoting specific values 
#'   to set each of the mediators in `M` to, for estimating the CDE. The length 
#'   of the vector MUST be exactly the same length as that of the vector `M`. 
#'   This argument is only used in the estimation of the CDE, not in any of the 
#'   other returned estimands.
#' @param interaction_DM A logical scalar indicating whether the outcome model 
#'   should include exposure-mediator interactions (interactions of the exposure 
#'   with each mediator if there is more than one mediator in `M`).
#' @param interaction_DC A logical scalar indicating whether both the outcome 
#'   model and each of the mediator models should include interactions of the 
#'   exposure with each covariate in `C`.
#' @param interaction_MC A logical scalar indicating whether the outcome model 
#'   should include interactions of each mediator in `M` with each covariate in 
#'   `C`.
#' @param weights_name A character scalar identifying the name of the weights 
#'   variable in `data`, if applicable (e.g., if you have---and want to 
#'   use---sampling weights).
#' 
#' @returns `linmed()` returns a list with the following elements:
#' \item{ATE}{A numeric scalar with the estimated total average treatment effect 
#'   for the exposure contrast `d - dstar`: ATE(`d`,`dstar`).}
#' \item{NDE}{A numeric scalar with the estimated natural direct effect for the 
#'   exposure contrast `d - dstar`: NDE(`d`,`dstar`).}
#' \item{NIE}{A numeric scalar with the estimated natural indirect effect for 
#'   the exposure contrast `d - dstar`: NIE(`d`,`dstar`).}
#' \item{CDE}{A numeric scalar with the estimated controlled direct effect for 
#'   the exposure contrast `d - dstar` and the mediator value(s) `m`: 
#'   CDE(`d`,`dstar`,`m`).}
#' \item{model_m}{A list with the model objects from each of the fitted mediator 
#'   models.}
#' \item{model_y}{The model object from the fitted outcome model.}
#' \item{miss_summary}{A data frame with counts of non-missing (`nmiss`) and 
#'   missing (`miss`) observations for each of the variables specified for `D`, 
#'   `M`, `Y`, and `C`.}
#' 
#' @export
#' 
#' @examples
#' # Example 1: Single mediator, no interactions
#' ## Prepare data
#' ## For convenience with this example, we will use complete cases
#' data(nlsy)
#' covariates <- c(
#'   "female",
#'   "black",
#'   "hispan",
#'   "paredu",
#'   "parprof",
#'   "parinc_prank",
#'   "famsize",
#'   "afqt3"
#' )
#' nlsy1 <- nlsy
#' key_variables1 <- c(
#'   "cesd_age40",
#'   "ever_unemp_age3539",
#'   "att22",
#'   covariates
#' )
#' nlsy1 <- nlsy1[complete.cases(nlsy1[,key_variables1]),]
#' nlsy1$std_cesd_age40 <- 
#'   (nlsy1$cesd_age40 - mean(nlsy1$cesd_age40)) / 
#'   sd(nlsy1$cesd_age40)
#' ## Estimate natural effects
#' linmed(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   C = c(
#'     "female",
#'     "black",
#'     "hispan",
#'     "paredu",
#'     "parprof",
#'     "parinc_prank",
#'     "famsize",
#'     "afqt3"
#'   )
#' )
#' 
#' # Example 2: With exposure-mediator, exposure-covariate, and 
#' # mediator-covariate interactions
#' linmed(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   C = c(
#'     "female",
#'     "black",
#'     "hispan",
#'     "paredu",
#'     "parprof",
#'     "parinc_prank",
#'     "famsize",
#'     "afqt3"
#'   ),
#'   interaction_DM = TRUE,
#'   interaction_DC = TRUE,
#'   interaction_MC = TRUE
#' )
#' 
#' # Example 3: Specifying a mediator value for the CDE other than (the default) 
#' # zero
#' linmed(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   C = c(
#'     "female",
#'     "black",
#'     "hispan",
#'     "paredu",
#'     "parprof",
#'     "parinc_prank",
#'     "famsize",
#'     "afqt3"
#'   ),
#'   m = 1
#' )
#' 
#' # Example 4: Incorporating sampling weights
#' linmed(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   C = c(
#'     "female",
#'     "black",
#'     "hispan",
#'     "paredu",
#'     "parprof",
#'     "parinc_prank",
#'     "famsize",
#'     "afqt3"
#'   ),
#'   weights_name = "weight"
#' )
#' 
#' # Example 5: Multiple mediators
#' ## Prepare data
#' nlsy2 <- nlsy
#' key_variables2 <- c(
#'   key_variables1,
#'   "log_faminc_adj_age3539"
#' )
#' nlsy2 <- nlsy2[complete.cases(nlsy2[,key_variables2]),]
#' nlsy2$std_cesd_age40 <- 
#'   (nlsy2$cesd_age40 - mean(nlsy2$cesd_age40)) / 
#'   sd(nlsy2$cesd_age40)
#' ## Estimate natural effects
#' linmed(
#'   data = nlsy2,
#'   D = "att22",
#'   M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
#'   Y = "std_cesd_age40",
#'   C = c(
#'     "female",
#'     "black",
#'     "hispan",
#'     "paredu",
#'     "parprof",
#'     "parinc_prank",
#'     "famsize",
#'     "afqt3"
#'   )
#' )
#' 
#' # Example 6: Specifying mediator values for the CDE, with multiple mediators
#' linmed(
#'   data = nlsy2,
#'   D = "att22",
#'   M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
#'   Y = "std_cesd_age40",
#'   C = c(
#'     "female",
#'     "black",
#'     "hispan",
#'     "paredu",
#'     "parprof",
#'     "parinc_prank",
#'     "famsize",
#'     "afqt3"
#'   ),
#'   m = c(1, 3.5)
#' )
linmed <- function(
    data,
    D,
    M,
    Y,
    C = NULL,
    d = 1, 
    dstar = 0,
    m = rep(0, length(M)),
    interaction_DM = FALSE,
    interaction_DC = FALSE,
    interaction_MC = FALSE,
    weights_name = NULL
) {
  df <- data
  
  # error checks
  if (!is.numeric(df[[D]])) {
    stop("Error: The exposure variable (identified by the string argument D in data) must be numeric.")
  }
  if (!is.numeric(df[[Y]])) {
    stop("Error: The outcome variable (identified by the string argument Y in data) must be numeric.")
  }
  for (M_k in M) {
    if (!is.numeric(df[[M_k]])) {
      stop("Error: Each of the mediator variables (identified by the string argument M in data) must be numeric.")
    }
  }
  if (length(M)!=length(m)) {
    stop("Error: The lengths of the M and m arguments must be the same.")
  }
  
  # check for missing data and create missing summary output
  key_vars <- c(D, M, Y, C)
  miss_summary <- sapply(
    key_vars,
    FUN = function(v) c(
      nmiss = sum(!is.na(df[[v]])),
      miss = sum(is.na(df[[v]]))
    )
  ) |>
    t() |>
    as.data.frame()
  if (max(miss_summary$miss)>0) {
    warning("Warning: There is missing data in at least one of the variables specified for D, M, Y, and C. See the miss_summary data frame in the output.")
  }
  
  # assign weights
  if (is.null(weights_name)) {
    weights <- rep(1, nrow(df))
  }
  else {
    weights <- df[[weights_name]]
  }
  
  # demean covariates
  # (This is useful specifically if there are any D*C or M*C interactions. 
  # However, it is harmless even if there are no such interactions, so there is 
  # no need to condition the demeaning on the interaction arguments.)
  for(covariate in C) df[[covariate]] <- demean(df[[covariate]], w=weights)
  
  # mediator model(s) predictors
  m_preds <- paste(
    c(C, D), 
    collapse = " + "
  )
  if (interaction_DC) {
    m_preds <- paste(
      m_preds,
      "+",
      paste(D, C, sep = ":", collapse = " + ")
    )
  }
  
  # outcome model predictors
  y_preds <- paste(
    c(C, D, M), 
    collapse = " + "
  )
  if (interaction_DM) {
    y_preds <- paste(
      y_preds,
      "+",
      paste(D, M, sep = ":", collapse = " + ")
    )
  }
  if (interaction_DC) {
    y_preds <- paste(
      y_preds,
      "+",
      paste(D, C, sep = ":", collapse = " + ")
    )
  }
  if (interaction_MC) {
    y_preds <- paste(
      y_preds,
      "+",
      paste(outer(M, C, FUN = "paste", sep = ":"), collapse = " + ")
    )
  }
  
  # specify formulas (as character strings) for mediator and outcome models
  m_forms <- purrr::map(M, ~ paste(.x, "~", m_preds))
  y_form <- paste(Y, "~", y_preds)
  
  # fit mediator and outcome models
  m_models <- purrr::map(m_forms, ~ lm(as.formula(.x), data = df, weights = weights))
  names(m_models) <- M
  y_model <- lm(as.formula(y_form), data = df, weights = weights)
  
  # map over K mediators (K>1 if M is multivariate)
  if (interaction_DM) {
    NDE_part <- purrr::map2_dbl(
      M, m_models, 
      function(M_k, M_model_k) y_model$coef[[paste0(D, ":", M_k)]] * (M_model_k$coef[["(Intercept)"]] + M_model_k$coef[[D]]*dstar)
    )
    NIE_part <- purrr::map2_dbl(
      M, m_models, 
      function(M_k, M_model_k) M_model_k$coef[[D]] * (y_model$coef[[M_k]] + y_model$coef[[paste0(D, ":", M_k)]]*d)
    )
    CDE_part <- purrr::map2_dbl(
      M, m, 
      function(M_k, m_k) y_model$coef[[paste0(D, ":", M_k)]] * m_k
    )
  }
  else {
    NDE_part <- 0
    NIE_part <- purrr::map2_dbl(
      M, m_models, 
      function(M_k, M_model_k) M_model_k$coef[[D]] * y_model$coef[[M_k]]
    )
    CDE_part <- 0
  }
  
  # compute NDE, NIE, ATE, and CDE estimates
  NDE <- (y_model$coef[[D]] + sum(NDE_part)) * (d - dstar)
  # ^ the above translates to the following (where k indexes the mediators and 
  # using the notation from the book):
  # 1. if there is no D*M interaction:
  # gamma_2 * (d - dstar)
  # 2. if there is a D*M interaction:
  # (gamma_2 + sum(gamma_4k * (beta_0k + beta_2k*dstar))) * (d - dstar)
  NIE <- sum(NIE_part) * (d - dstar)
  # ^ the above translates to the following:
  # 1. if there is no D*M interaction:
  # sum(beta_2k * gamma_3k) * (d - dstar)
  # 2. if there is a D*M interaction:
  # sum(beta_2k * (gamma_3k + gamma_4k*d)) * (d - dstar)
  ATE <- NDE + NIE
  CDE <- (y_model$coef[[D]] + sum(CDE_part)) * (d - dstar)
  # ^ the above translates to the following:
  # 1. if there is no D*M interaction:
  # gamma_2 * (d - dstar)
  # 2. if there is a D*M interaction:
  # (gamma_2 + sum(gamma_4k * m_k)) * (d - dstar)
  
  # compile and output
  out <- list(
    ATE = ATE, NDE = NDE, NIE = NIE, CDE = CDE,
    model_m = m_models,
    model_y = y_model,
    miss_summary = miss_summary
  )
  return(out)
}

