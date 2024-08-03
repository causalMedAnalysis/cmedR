#' Inverse probability weighting (IPW) estimator for natural effects
#' 
#' @description
#' `ipwmed()` uses the inverse probability weighting (IPW) estimator to estimate 
#' the total effect (ATE), natural direct effect (NDE), and natural indirect 
#' effect (NIE). The function supports estimation of both univariate and 
#' multivariate natural effects.
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
#'   `M = c("ever_unemp_age3539", "log_faminc_adj_age3539")`.
#' @param Y A character scalar identifying the name of the outcome variable in 
#'   `data`. `Y` is a character string, but the outcome variable it identifies 
#'   must be numeric.
#' @param formula1_string A character scalar for the formula to be fitted for a 
#'   GLM of the exposure given baseline covariates (denoted in the book as 
#'   f(D|C)). E.g., `formula1_string = "att22~female+black+paredu"`.
#' @param formula2_string A character scalar for the formula to be fitted for a 
#'   GLM of the exposure given baseline covariates and the mediator (denoted in 
#'   the book as s(D|C,M)). E.g., 
#'   `formula2_string = "att22~female+black+paredu+ever_unemp_age3539"`.
#' @param base_weights_name A character scalar identifying the name of the base 
#'   weights variable in `data`, if applicable (e.g., if you have---and want to 
#'   use---sampling weights).
#' @param stabilize A logical scalar indicating whether the IPW weights should 
#'   be stabilized (multiplied by the marginal probabilities of the exposure).
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
#' @returns `ipwmed()` returns a list with the following elements:
#' \item{ATE}{A numeric scalar with the estimated total average treatment effect 
#'   for the exposure contrast 1 - 0: ATE(1,0).}
#' \item{NDE}{A numeric scalar with the estimated natural direct effect for the 
#'   exposure contrast 1 - 0: NDE(1,0).}
#' \item{NIE}{A numeric scalar with the estimated natural indirect effect for 
#'   the exposure contrast 1 - 0: NIE(1,0).}
#' \item{weights1}{A numeric vector with the final inverse probability weights 
#'   for w_1, as defined in the book.}
#' \item{weights2}{A numeric vector with the final inverse probability weights 
#'   for w_2, as defined in the book.}
#' \item{weights3}{A numeric vector with the final inverse probability weights 
#'   for w_3, as defined in the book.}
#' \item{model_d1}{The model object from the first fitted exposure model 
#'   (of the exposure given baseline covariates, denoted in the book as f(D|C)), 
#'   corresponding to `formula1_string`.}
#' \item{model_d2}{The model object from the second fitted exposure model 
#'   (of the exposure given baseline covariates and the mediator(s), denoted in 
#'   the book as s(D|C,M)), corresponding to `formula2_string`.}
#' 
#' @export
#' 
#' @examples
#' # Example 1: Single mediator
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
#' ## Estimate natural effects
#' out1 <- ipwmed(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
#'   formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3+ever_unemp_age3539"
#' )
#' head(out1,3)
#' 
#' # Example 2: Incorporating sampling weights
#' out2 <- ipwmed(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
#'   formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3+ever_unemp_age3539",
#'   base_weights_name = "weight"
#' )
#' head(out2,3)
#' 
#' # Example 3: Multiple mediators
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
#' out3 <- ipwmed(
#'   data = nlsy2,
#'   D = "att22",
#'   M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
#'   Y = "std_cesd_age40",
#'   formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
#'   formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3+ever_unemp_age3539+log_faminc_adj_age3539"
#' )
#' head(out3,3)
ipwmed <- function(
    data,
    D,
    M,
    Y,
    formula1_string,
    formula2_string,
    base_weights_name = NULL,
    stabilize = TRUE,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99
) {
  d <- 1
  dstar <- 0
  
  # error and warning checks
  if (!is.numeric(data[[D]])) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be numeric."), collapse = "\n"))
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
  for (M_k in M) {
    if (grepl(pattern = M_k, x = formula1_string, fixed = TRUE)) {
      warning(paste(strwrap("Warning: Check whether (each) mediator variable is among the predictors in the formula1_string. The mediator(s) should not be among the predictors in the formula1_string."), collapse = "\n"))
    }
    if (!grepl(pattern = M_k, x = formula2_string, fixed = TRUE)) {
      warning(paste(strwrap("Warning: Check whether (each) mediator variable is among the predictors in the formula2_string. The mediator(s) should be among the predictors in the formula2_string."), collapse = "\n"))
    }
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
  
  # fit exposure models
  d_model1 <- glm(
    as.formula(formula1_string),
    data = data,
    family = binomial(link = "logit"),
    weights = base_weights_rsc
  )
  d_model2 <- glm(
    as.formula(formula2_string),
    data = data,
    family = binomial(link = "logit"),
    weights = base_weights_rsc
  )
  # enforcing a no-missing-data rule
  if (nrow(d_model1$model)!=nrow(data) | nrow(d_model2$model)!=nrow(data)) {
    stop(paste(strwrap("Error: Please remove observations with missing/NA values for the exposure, mediator, outcome, or covariates."), collapse = "\n"))
  }
  
  # predict exposure probabilities
  ps1_d <- predict(d_model1, newdata = data, type = "response") 
  ps2_d <- predict(d_model2, newdata = data, type = "response")
  ps1_dstar <- 1 - ps1_d
  ps2_dstar <- 1 - ps2_d
  marg_prob_d <- weighted.mean(ps1_d, base_weights_rsc)
  marg_prob_dstar <- 1 - marg_prob_d
  
  # for convenience, create logical vectors to identify subgroups
  group_dstar <- data[[D]]==dstar & base_weights_rsc>0
  group_d <- data[[D]]==d & base_weights_rsc>0
  
  # create IPWs
  w1 <- ifelse(
    group_dstar,
    1 / ps1_dstar,
    0
  )
  w2 <- ifelse(
    group_d,
    1 / ps1_d,
    0
  )
  w3 <- ifelse(
    group_d,
    ps2_dstar / (ps2_d * ps1_dstar),
    0
  )
  
  # stabilize IPWs
  if (stabilize) {
    w1 <- w1 * marg_prob_dstar
    w2 <- w2 * marg_prob_d
    w3 <- w3 * marg_prob_d
  }
  
  # censor IPWs (among the appropriate subgroups with non-zero IPWs)
  if (censor) {
    w1[group_dstar] <- trimQ(w1[group_dstar], low = censor_low, high = censor_high)
    w2[group_d] <- trimQ(w2[group_d], low = censor_low, high = censor_high)
    w3[group_d] <- trimQ(w3[group_d], low = censor_low, high = censor_high)
  }
  
  # multiply IPWs by rescaled base weights
  final_w1 <- w1 * base_weights_rsc
  final_w2 <- w2 * base_weights_rsc
  final_w3 <- w3 * base_weights_rsc
  
  # estimate effects
  Ehat_Ydstar_Mdstar <- weighted.mean(data[[Y]], final_w1)
  Ehat_Yd_Md <- weighted.mean(data[[Y]], final_w2)
  Ehat_Yd_Mdstar <- weighted.mean(data[[Y]], final_w3)
  ATE <- Ehat_Yd_Md - Ehat_Ydstar_Mdstar
  NDE <- Ehat_Yd_Mdstar - Ehat_Ydstar_Mdstar
  NIE <- Ehat_Yd_Md - Ehat_Yd_Mdstar
  
  # compile and output
  out <- list(
    ATE = ATE,
    NDE = NDE,
    NIE = NIE,
    weights1 = final_w1,
    weights2 = final_w2,
    weights3 = final_w3,
    model_d1 = d_model1,
    model_d2 = d_model2
  )
  return(out)
}

