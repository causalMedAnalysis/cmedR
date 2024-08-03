#' Controlled direct effect inverse probability weighting (IPW) estimator
#' 
#' @description
#' `ipwcde()` uses the inverse probability weighting (IPW) estimator to estimate 
#' the controlled direct effect (CDE). Note that unlike the ipwmed() function, 
#' ipwcde() requires a single mediator. Multiple mediators are not supported.
#' 
#' @details
#' TEMPORARY PLACEHOLDER
#' 
#' @param data A data frame.
#' @param D A character scalar identifying the name of the exposure variable in 
#'   `data`. `D` is a character string, but the exposure variable it identifies 
#'   must be numeric and binary, consisting only of the values 0 and 1.
#' @param M A character scalar identifying the name of the mediator variable 
#'   (only one mediator variable is supported) in `data`. `D` is a character 
#'   string, but the exposure variable it identifies must be numeric and binary, 
#'   consisting only of the values 0 and 1.
#' @param Y A character scalar identifying the name of the outcome variable in 
#'   `data`. `Y` is a character string, but the outcome variable it identifies 
#'   must be numeric.
#' @param m A numeric scalar denoting a specific value to set the mediator `M` 
#'   to, for estimating the CDE.
#' @param formula_D_string A character scalar for the formula to be fitted for a 
#'   GLM of the exposure given baseline covariates (denoted in the book as 
#'   f(D|C)). E.g., `formula_D_string = "att22~female+black+paredu"`.
#' @param formula_M_string A character scalar for the formula to be fitted for a 
#'   GLM of the mediator given baseline covariates and the exposure (denoted in 
#'   the book as g(M|C,D)). E.g., 
#'   `formula_M_string = "ever_unemp_age3539~female+black+paredu+att22"`.
#' @param base_weights_name A character scalar identifying the name of the base 
#'   weights variable in `data`, if applicable (e.g., if you have---and want to 
#'   use---sampling weights).
#' @param stabilize A logical scalar indicating whether the IPW weights should 
#'   be stabilized.
#' @param censor A logical scalar indicating whether the IPW weights should 
#'   be censored.
#' @param censor_low,censor_high A pair of arguments, each a numeric scalar 
#'   denoting a probability with values in [0,1]. If the `censor` argument is 
#'   TRUE, then IPW weights below the `censor_low` quantile will be 
#'   bottom-coded, and IPW weights above the `censor_high` quantile will be 
#'   top-coded (before multiplying by a rescaled version of the base weights, if 
#'   applicable). E.g., if the default options of `censor_low = 0.01` and 
#'   `censor_high = 0.99` are used, then the IPW weights will be censored at 
#'   their 1st and 99th percentiles in the data.
#' 
#' @returns `ipwcde()` returns a list with the following elements:
#' \item{CDE}{A numeric scalar with the estimated controlled direct effect for 
#'   the exposure contrast 1 - 0: CDE(1,0,`m`).}
#' \item{weights}{A numeric vector with the final inverse probability weights.}
#' \item{model_d}{The model object from the fitted exposure model (of the 
#'   exposure given baseline covariates, denoted in the book as f(D|C)), 
#'   corresponding to `formula_D_string`.}
#' \item{model_m}{The model object from the fitted mediator model (of the 
#'   mediator given baseline covariates and the exposure, denoted in the book 
#'   as g(M|C,D)), corresponding to `formula_M_string`.}
#' 
#' @export
#' 
#' @examples
#' # Example 1
#' ## Prepare data
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
#' ## Estimate CDE for m=1
#' out1 <- ipwcde(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   m = 1,
#'   formula_D_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
#'   formula_M_string = "ever_unemp_age3539~att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3"
#' )
#' head(out1,1)
#' 
#' # Example 2: Incorporating sampling weights
#' out2 <- ipwcde(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   m = 1,
#'   formula_D_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
#'   formula_M_string = "ever_unemp_age3539~att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
#'   base_weights_name = "weight"
#' )
#' head(out2,1)
ipwcde <- function(
    data,
    D,
    M,
    Y,
    m = 0,
    formula_D_string,
    formula_M_string,
    base_weights_name = NULL,
    stabilize = TRUE,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99
) {
  d <- 1
  dstar <- 0
  
  # error and warning checks
  if (length(M)>1) {
    stop(paste(strwrap("Error: Unlike the ipwmed() function, ipwcde() requires a single mediator. Multiple mediators are not supported."), collapse = "\n"))
  }
  if (length(m)>1) {
    stop(paste(strwrap("Error: Please specify a single value for the argument m."), collapse = "\n"))
  }
  if (!is.numeric(data[[D]])) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be numeric."), collapse = "\n"))
  }
  if (!is.numeric(data[[M]])) {
    stop(paste(strwrap("Error: The mediator variable (identified by the string argument M in data) must be numeric."), collapse = "\n"))
  }
  if (!is.numeric(data[[Y]])) {
    stop(paste(strwrap("Error: The outcome variable (identified by the string argument Y in data) must be numeric."), collapse = "\n"))
  }
  if (any(is.na(data[[D]]))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the exposure variable (identified by the string argument D in data)."), collapse = "\n"))
  }
  if (any(! data[[D]] %in% c(0,1))) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be a numeric variable consisting only of the values 0 or 1. There is at least one observation in the data that does not meet this criteria."), collapse = "\n"))
  }
  if (any(is.na(data[[M]]))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the mediator variable (identified by the string argument M in data)."), collapse = "\n"))
  }
  if (any(! data[[M]] %in% c(0,1))) {
    stop(paste(strwrap("Error: Only binary mediator variables are supported. The mediator variable (identified by the string argument M in data) must be a numeric variable consisting only of the values 0 or 1. There is at least one observation in the data that does not meet this criteria."), collapse = "\n"))
  }
  if (! m %in% c(0,1)) {
    stop(paste(strwrap("Error: Because only binary mediator variables are supported, the selected mediator value for the CDE (in the m argument) must be either 0 or 1."), collapse = "\n"))
  }
  if (grepl(pattern = M, x = formula_D_string, fixed = TRUE)) {
    warning(paste(strwrap("Warning: Check whether the mediator variable is among the predictors in the formula_D_string. The mediator should not be among the predictors in the formula_D_string."), collapse = "\n"))
  }
  if (!grepl(pattern = D, x = formula_M_string, fixed = TRUE)) {
    warning(paste(strwrap("Warning: Check whether the exposure variable is among the predictors in the formula_M_string. The exposure should be among the predictors in the formula_M_string."), collapse = "\n"))
  }
  # ^ Note that the grepl-based warning checks are fairly simple, based solely 
  # on whether the string is detected. For now, we are not using more complex 
  # checks searching for full words in the model formula.
  
  # assign base weights
  if (is.null(base_weights_name)) {
    base_weights <- rep(1, nrow(data))
  }
  else {
    base_weights <- data[[base_weights_name]]
  }
  if (!is.null(base_weights_name) & any(is.na(base_weights))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the base weights variable (identified by the string argument base_weights_name in data). If that observation should not receive a positive weight, please replace the NA value with a zero before proceeding."), collapse = "\n"))
  }
  
  # rescale base weights
  base_weights_rsc <- base_weights / mean(base_weights)
  
  # fit specified models
  d_model <- glm(
    as.formula(formula_D_string),
    data = data,
    family = binomial(link = "logit"),
    weights = base_weights_rsc
  )
  m_model <- glm(
    as.formula(formula_M_string),
    data = data,
    family = binomial(link = "logit"),
    weights = base_weights_rsc
  )
  # enforcing a no-missing-data rule
  if (nrow(d_model$model)!=nrow(data) | nrow(m_model$model)!=nrow(data)) {
    stop(paste(strwrap("Error: Please remove observations with missing/NA values for the exposure, mediator, outcome, or covariates."), collapse = "\n"))
  }
  
  # additionally fit mediator model without covariates
  m_model_no_cov <- glm(
    as.formula(paste0(M,"~",D)),
    data = data,
    family = binomial(link = "logit"),
    weights = base_weights_rsc
  )
  
  # predict exposure and mediator probabilities
  ps_D1_C <- predict(d_model, newdata = data, type = "response")
  ps_M1_CD <- predict(m_model, newdata = data, type = "response")
  ps_M1_D <- predict(m_model_no_cov, newdata = data, type = "response")
  ps_D_C <- ifelse(as.logical(data[[D]]), ps_D1_C, 1-ps_D1_C)
  ps_M_CD <- ifelse(as.logical(data[[M]]), ps_M1_CD, 1-ps_M1_CD)
  marg_prob_D1 <- weighted.mean(data[[D]], base_weights_rsc)
  marg_prob_D <- ifelse(
    as.logical(data[[D]]),
    marg_prob_D1,
    1 - marg_prob_D1
  )
  marg_prob_M_D <- ifelse(
    as.logical(data[[M]]),
    ps_M1_D,
    1 - ps_M1_D
  )
  
  # identify observations with positive base weights
  group_pos <- base_weights_rsc>0
  
  # create IPWs
  w4 <- 1 / (ps_M_CD * ps_D_C)
  
  # stabilize IPWs
  if (stabilize) {
    w4 <- w4 * marg_prob_M_D * marg_prob_D
  }
  
  # censor IPWs (among observations with positive base weights)
  if (censor) {
    w4[group_pos] <- trimQ(w4[group_pos], low = censor_low, high = censor_high)
  }
  
  # multiply IPWs by rescaled base weights
  final_w4 <- w4 * base_weights_rsc
  
  # estimate effects
  group_d_m <- data[[D]]==d & data[[M]]==m
  group_dstar_m <- data[[D]]==dstar & data[[M]]==m
  CDE <- 
    weighted.mean(data[[Y]][group_d_m], final_w4[group_d_m]) - 
    weighted.mean(data[[Y]][group_dstar_m], final_w4[group_dstar_m])
  
  # compile and output
  out <- list(
    CDE = CDE,
    weights = final_w4,
    model_d = d_model,
    model_m = m_model
  )
  return(out)
}

