#' Controlled direct effect regression imputation estimator
#' 
#' @description
#' `impcde()` estimates the controlled direct effect (CDE) with a regression 
#' imputation estimator.
#' 
#' @details
#' TEMPORARY PLACEHOLDER
#' 
#' @param data A data frame.
#' @param model_y A model object for the fitted outcome model.
#' @param D A character scalar identifying the name of the exposure variable in 
#'   `data`.
#' @param M A character scalar identifying the name of the mediator variable in 
#'   `data`.
#' @param d,dstar A pair of arguments, each a numeric scalar denoting a specific 
#'   value of the exposure `D`. The exposure contrast of interest is 
#'   `d - dstar`.
#' @param m A numeric scalar denoting a specific value to set the mediator `M` 
#'   to, for estimating the CDE.
#' @param weights_name A character scalar identifying the name of the weights 
#'   variable in `data`, if applicable (e.g., if you have---and want to 
#'   use---sampling weights).
#' 
#' @returns `impcde()` returns a numeric scalar with the estimated controlled 
#' direct effect for the exposure contrast `d - dstar` and the mediator value 
#' `m`: CDE(`d`,`dstar`,`m`).
#' 
#' @export
#' 
#' @examples
#' # Example 1
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
#' key_variables <- c(
#'   "cesd_age40",
#'   "ever_unemp_age3539",
#'   "att22",
#'   covariates
#' )
#' nlsy <- nlsy[complete.cases(nlsy[,key_variables]),]
#' nlsy$std_cesd_age40 <- 
#'   (nlsy$cesd_age40 - mean(nlsy$cesd_age40)) / 
#'   sd(nlsy$cesd_age40)
#' ## Fit model
#' mod1 <- lm(
#'   std_cesd_age40 ~ ever_unemp_age3539 + att22 + 
#'     female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3,
#'   data = nlsy
#' )
#' ## Estimate CDE for m=1
#' impcde(
#'   data = nlsy,
#'   model_y = mod1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   m = 1
#' )
#' 
#' # Example 2: Incorporating sampling weights
#' impcde(
#'   data = nlsy,
#'   model_y = mod1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   m = 1,
#'   weights_name = "weight"
#' )
impcde <- function(
    data,
    model_y,
    D,
    M,
    d = 1,
    dstar = 0,
    m = 0,
    weights_name = NULL
) {
  # assign weights
  if (is.null(weights_name)) {
    weights <- rep(1, nrow(data))
  }
  else {
    weights <- data[[weights_name]]
  }
  
  # predict Y(d, m)
  gdata <- data
  gdata[[D]] <- d
  gdata[[M]] <- m
  Yhat_d_m <- predict(model_y, gdata, type="response")
  
  # predict Y(dstar, m)
  gdata[[D]] <- dstar
  Yhat_dstar_m <- predict(model_y, gdata, type="response")
  
  # estimate CDE(d, dstar, m) and output
  cde <- weighted.mean(Yhat_d_m - Yhat_dstar_m, w=weights)
  return(cde)
}

