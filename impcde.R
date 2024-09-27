#' Controlled direct effect regression imputation estimator: inner function
#' 
#' @description
#' Internal function used within `impcde()`. See the `impcde()` function 
#' documentation for a description of shared function arguments.
#' 
#' @noRd
impcde_inner <- function(
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
#' @param boot A logical scalar indicating whether the function will perform the 
#'   nonparametric bootstrap and return a two-sided confidence interval and 
#'   p-value.
#' @param boot_reps An integer scalar for the number of bootstrap replications 
#'   to perform.
#' @param boot_conf_level A numeric scalar for the confidence level of the 
#'   bootstrap interval.
#' @param boot_seed An integer scalar specifying the random-number seed used in 
#'   bootstrap resampling.
#' @param boot_parallel A logical scalar indicating whether the bootstrap will 
#'   be performed with a parallelized loop, with the goal of reducing runtime. 
#'   Parallelized computing, as implemented in this function, requires that you 
#'   have each of the following R packages installed: `doParallel`, `doRNG`, and 
#'   `foreach`. (However, you do not need to load/attach these three packages 
#'   with the `library` function prior to running this function.)
#' @param boot_cores An integer scalar specifying the number of CPU cores on 
#'   which the parallelized bootstrap will run. This argument only has an effect 
#'   if you requested a parallelized bootstrap (i.e., only if `boot` is TRUE and 
#'   `boot_parallel` is TRUE). By default, `boot_cores` is equal to the greater 
#'   of two values: (a) one and (b) the number of available CPU cores minus two. 
#'   If `boot_cores` equals one, then the bootstrap loop will not be 
#'   parallelized (regardless of whether `boot_parallel` is TRUE).
#' 
#' @returns By default, `impcde()` returns a numeric scalar with the estimated 
#' controlled direct effect for the exposure contrast `d - dstar` and the 
#' mediator value `m`: CDE(`d`,`dstar`,`m`).
#' 
#' If, however, you request the bootstrap (by setting the `boot` argument to 
#' TRUE), then the function instead returns a list with the following elements:
#' \item{CDE}{A numeric scalar with the point estimate for the CDE.}
#' \item{ci_CDE}{A numeric vector with the bootstrap confidence interval for the 
#'   CDE.}
#' \item{pvalue_CDE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the CDE is different from zero, as computed from the bootstrap.}
#' \item{boot_CDE}{A numeric vector of length `boot_reps` comprising the CDE 
#'   estimates from all replicate samples created in the bootstrap.}
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
#' ## Re-fit model with sampling weights
#' mod2 <- lm(
#'   std_cesd_age40 ~ ever_unemp_age3539 + att22 + 
#'     female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3,
#'   data = nlsy,
#'   weights = nlsy$weight
#' )
#' ## Estimate CDE for m=1
#' impcde(
#'   data = nlsy,
#'   model_y = mod2,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   m = 1,
#'   weights_name = "weight"
#' )
#' 
#' # Example 3: Perform a nonparametric bootstrap, with 2,000 replications
#' \dontrun{
#'   impcde(
#'     data = nlsy,
#'     model_y = mod1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     m = 1,
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234
#'   )
#' }
#' 
#' # Example 4: Parallelize the bootstrap, to attempt to reduce runtime
#' # Note that this requires you to have installed the `doParallel`, `doRNG`, 
#' # and `foreach` packages.
#' \dontrun{
#'   impcde(
#'     data = nlsy,
#'     model_y = mod1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     m = 1,
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234,
#'     boot_parallel = TRUE
#'   )
#' }
impcde <- function(
  data,
  model_y,
  D,
  M,
  d = 1,
  dstar = 0,
  m = 0,
  weights_name = NULL,
  boot = FALSE,
  boot_reps = 1000,
  boot_conf_level = 0.95,
  boot_seed = NULL,
  boot_parallel = FALSE,
  boot_cores = max(c(parallel::detectCores()-2,1))
) {
  # load data
  data_outer <- data
  
  
  # create adjusted boot_parallel logical
  boot_parallel_rev <- ifelse(boot_cores>1, boot_parallel, FALSE)
  
  
  # preliminary error/warning checks
  if (boot) {
    if (boot_parallel & boot_cores==1) {
      warning(paste(strwrap("Warning: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but you do not have enough cores available for parallelization. The bootstrap will proceed without parallelization."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("doParallel", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'doParallel' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("doRNG", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'doRNG' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("foreach", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'foreach' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
  }
  
  
  # compute point estimate
  est <- impcde_inner(
    data = data_outer,
    model_y = model_y,
    D = D,
    M = M,
    d = d,
    dstar = dstar,
    m = m,
    weights_name = weights_name
  )
  
  
  # bootstrap, if requested
  if (boot) {
    # bootstrap function
    boot_fnc <- function() {
      # sample from the data with replacement
      boot_data <- data_outer[sample(nrow(data_outer), size = nrow(data_outer), replace = TRUE), ]
      
      # re-fit the outcome model on the new sample
      if (is.null(weights_name)) {
        model_y_update <- update(model_y, data = boot_data)
      }
      else {
        boot_weights <- boot_data[[weights_name]]
        model_y_update <- update(model_y, data = boot_data, weights = boot_weights)
      }
      
      # compute point estimate in the replicate sample
      impcde_inner(
        data = boot_data,
        model_y = model_y_update,
        D = D,
        M = M,
        d = d,
        dstar = dstar,
        m = m,
        weights_name = weights_name
      )
    }
    
    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster, 
        varlist = c("impcde_inner"),
        envir = environment()
      )
      `%dopar%` <- foreach::`%dopar%`
    }
    
    # set seed
    if (!is.null(boot_seed)) {
      set.seed(boot_seed)
      if (boot_parallel) {
        doRNG::registerDoRNG(boot_seed)
      }
    }
    
    # compute estimates for each replicate sample
    if (boot_parallel_rev) {
      boot_CDE <- foreach::foreach(i = 1:boot_reps, .combine = c) %dopar% {
        boot_fnc()
      }
    }
    else {
      boot_CDE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_CDE[i] <- boot_fnc()
      }
    }
    
    # clean up
    if (boot_parallel_rev) {
      parallel::stopCluster(x_cluster)
      rm(x_cluster)
    }
    
    # compute bootstrap confidence intervals 
    # from percentiles of the bootstrap distributions
    boot_alpha <- 1 - boot_conf_level
    boot_ci_probs <- c(
      boot_alpha/2,
      1 - boot_alpha/2
    )
    boot_ci <- function(x) {
      quantile(x, probs=boot_ci_probs)
    }
    ci_CDE <- boot_ci(boot_CDE)
    
    # compute two-tailed bootstrap p-values
    boot_pval <- function(x) {
      2 * min(
        mean(x < 0),
        mean(x > 0)
      )
    }
    pvalue_CDE <- boot_pval(boot_CDE)
  }
  
  
  # final output
  if (boot) {
    out <- list(
      CDE = est,
      ci_CDE = ci_CDE,
      pvalue_CDE = pvalue_CDE,
      boot_CDE = boot_CDE
    )
  }
  else {
    out <- est
  }
  return(out)
}

