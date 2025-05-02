#' Linear models estimator for natural effects: inner function
#' 
#' @description
#' Internal function used within `linmed()`. See the `linmed()` function 
#' documentation for a description of shared function arguments. Here, we will 
#' only document the one argument that is not shared by `linmed_inner()` and 
#' `linmed()`: the `minimal` argument.
#' 
#' @param minimal A logical scalar indicating whether the function should 
#'   return only a minimal set of output. The `linmed()` function uses the 
#'   default of FALSE when calling `linmed_inner()` to generate the point 
#'   point estimates and sets the argument to TRUE when calling `linmed_inner()` 
#'   to perform the bootstrap.
#' 
#' @noRd
linmed_inner <- function(
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
    weights_name = NULL,
    minimal = FALSE
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
  if (!minimal) {
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
  m_forms <- lapply(M, function(x) paste(x, "~", m_preds))
  y_form <- paste(Y, "~", y_preds)
  
  # fit mediator and outcome models
  m_models <- lapply(m_forms, function(x) lm(as.formula(x), data = df, weights = weights))
  names(m_models) <- M
  y_model <- lm(as.formula(y_form), data = df, weights = weights)
  
  # map over K mediators (K>1 if M is multivariate)
  if (interaction_DM) {
    NDE_part <- mapply(
      function(M_k, M_model_k) y_model$coef[[paste0(D, ":", M_k)]] * (M_model_k$coef[["(Intercept)"]] + M_model_k$coef[[D]]*dstar),
      M, m_models
    )
    NIE_part <- mapply(
      function(M_k, M_model_k) M_model_k$coef[[D]] * (y_model$coef[[M_k]] + y_model$coef[[paste0(D, ":", M_k)]]*d),
      M, m_models
    )
    CDE_part <- mapply(
      function(M_k, m_k) y_model$coef[[paste0(D, ":", M_k)]] * m_k,
      M, m
    )
  }
  else {
    NDE_part <- 0
    NIE_part <- mapply(
      function(M_k, M_model_k) M_model_k$coef[[D]] * y_model$coef[[M_k]],
      M, m_models
    )
    CDE_part <- 0
  }
  
  # compute NDE, NIE, ATE, and CDE estimates
  NDE <- (y_model$coef[[D]] + sum(NDE_part)) * (d - dstar)
  # ^ the above translates to the following (where k indexes the mediators and 
  # using the notation from the book):
  # 1. if there is no D*M interaction:
  #    gamma_2 * (d - dstar)
  # 2. if there is a D*M interaction:
  #    (gamma_2 + sum(gamma_4k * (beta_0k + beta_2k*dstar))) * (d - dstar)
  NIE <- sum(NIE_part) * (d - dstar)
  # ^ the above translates to the following:
  # 1. if there is no D*M interaction:
  #    sum(beta_2k * gamma_3k) * (d - dstar)
  # 2. if there is a D*M interaction:
  #    sum(beta_2k * (gamma_3k + gamma_4k*d)) * (d - dstar)
  ATE <- NDE + NIE
  CDE <- (y_model$coef[[D]] + sum(CDE_part)) * (d - dstar)
  # ^ the above translates to the following:
  # 1. if there is no D*M interaction:
  #    gamma_2 * (d - dstar)
  # 2. if there is a D*M interaction:
  #    (gamma_2 + sum(gamma_4k * m_k)) * (d - dstar)
  
  # compile and output
  if (minimal) {
    out <- list(
      ATE = ATE,
      NDE = NDE,
      NIE = NIE,
      CDE = CDE
    )
  }
  else {
    out <- list(
      ATE = ATE,
      NDE = NDE,
      NIE = NIE,
      CDE = CDE,
      model_m = m_models,
      model_y = y_model,
      miss_summary = miss_summary
    )
  }
  return(out)
}






#' Linear models estimator for natural effects
#' 
#' @description
#' `linmed()` uses a product-of-coefficients estimator, based on linear 
#' models for the mediator(s) and outcome, to estimate the total effect (ATE), 
#' natural direct effect (NDE), natural indirect effect (NIE), and controlled direct effect (CDE). 
#' The function supports estimation of both univariate and multivariate natural effects.
#' 
#' @details
#' `linmed()` performs causal mediation analysis using linear models for both the mediator(s) 
#' and outcome, and it computes inferential statistics using the nonparametric bootstrap. When a 
#' single mediator is specified, it estimates total, natural direct, and natural indirect effects 
#' using two linear models: a model for the mediator conditional on treatment and baseline covariates 
#' after centering them around their sample means, and a model for the outcome conditional on treatment, 
#' the mediator, and the baseline covariates after centering them around their sample means.
#' 
#' When multiple mediators are specified, `linmed()` provides estimates for the total effect and then 
#' for the multivariate natural direct and indirect effects operating through the entire set of 
#' mediators considered together. To this end, it fits separate models for each mediator conditional 
#' on treatment and the baseline covariates after centering them around their sample means, and then 
#' a model for the outcome conditional on treatment, all the mediators, and the baseline covariates 
#' after centering them around their sample means.
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
#' @param boot A logical scalar indicating whether the function will perform the 
#'   nonparametric bootstrap and return two-sided confidence intervals and 
#'   p-values.
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
#'   with the `library` function prior to running this function.) Note that the 
#'   results of the parallelized bootstrap may differ slightly from the 
#'   non-parallelized bootstrap, even if you specify the same seed, due to 
#'   differences in how the seed is processed by the two methods.
#' @param boot_cores An integer scalar specifying the number of CPU cores on 
#'   which the parallelized bootstrap will run. This argument only has an effect 
#'   if you requested a parallelized bootstrap (i.e., only if `boot` is TRUE and 
#'   `boot_parallel` is TRUE). By default, `boot_cores` is equal to the greater 
#'   of two values: (a) one and (b) the number of available CPU cores minus two. 
#'   If `boot_cores` equals one, then the bootstrap loop will not be 
#'   parallelized (regardless of whether `boot_parallel` is TRUE).
#' 
#' @returns By default, `linmed()` returns a list with the following elements:
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
#' If you request the bootstrap (by setting the `boot` argument to TRUE), then 
#' the function returns all of the elements listed above, as well as the 
#' following additional elements:
#' \item{ci_ATE}{A numeric vector with the bootstrap confidence interval for the 
#'   total average treatment effect (ATE).}
#' \item{ci_NDE}{A numeric vector with the bootstrap confidence interval for the 
#'   natural direct effect (NDE).}
#' \item{ci_NIE}{A numeric vector with the bootstrap confidence interval for the 
#'   natural indirect effect (NIE).}
#' \item{ci_CDE}{A numeric vector with the bootstrap confidence interval for the 
#'   controlled direct effect (CDE).}
#' \item{pvalue_ATE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the ATE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_NDE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the NDE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_NIE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the NIE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_CDE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the CDE is different from zero, as computed from the bootstrap.}
#' \item{boot_ATE}{A numeric vector of length `boot_reps` comprising the ATE 
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_NDE}{A numeric vector of length `boot_reps` comprising the NDE 
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_NIE}{A numeric vector of length `boot_reps` comprising the NIE 
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_CDE}{A numeric vector of length `boot_reps` comprising the CDE 
#'   estimates from all replicate samples created in the bootstrap.}
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
#' 
#' # Example 7: Perform a nonparametric bootstrap, with 2,000 replications
#' \dontrun{
#'   linmed(
#'     data = nlsy1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     Y = "std_cesd_age40",
#'     C = c(
#'       "female",
#'       "black",
#'       "hispan",
#'       "paredu",
#'       "parprof",
#'       "parinc_prank",
#'       "famsize",
#'       "afqt3"
#'     ),
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234
#'   )
#' }
#' 
#' # Example 8: Parallelize the bootstrap, to attempt to reduce runtime
#' # Note that this requires you to have installed the `doParallel`, `doRNG`, 
#' # and `foreach` packages.
#' \dontrun{
#'   linmed(
#'     data = nlsy1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     Y = "std_cesd_age40",
#'     C = c(
#'       "female",
#'       "black",
#'       "hispan",
#'       "paredu",
#'       "parprof",
#'       "parinc_prank",
#'       "famsize",
#'       "afqt3"
#'     ),
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234,
#'     boot_parallel = TRUE
#'   )
#' }
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
    if (!is.null(weights_name)) {
      warning(paste(strwrap("Warning: You requested a bootstrap, but your design includes sampling weights. Note that this function does not internally rescale sampling weights for use with the bootstrap, and it does not account for any stratification or clustering in your sample design. Failure to properly adjust the bootstrap sampling to account for a complex sample design that requires weighting could lead to invalid inferential statistics."), collapse = "\n"))
    }
  }
  
  
  # compute point estimates
  est <- linmed_inner(
    data = data_outer,
    D = D,
    M = M,
    Y = Y,
    C = C,
    d = d, 
    dstar = dstar,
    m = m,
    interaction_DM = interaction_DM,
    interaction_DC = interaction_DC,
    interaction_MC = interaction_MC,
    weights_name = weights_name,
    minimal = FALSE
  )
  
  
  # bootstrap, if requested
  if (boot) {
    # bootstrap function
    boot_fnc <- function() {
      # sample from the data with replacement
      boot_data <- data_outer[sample(nrow(data_outer), size = nrow(data_outer), replace = TRUE), ]
      
      # compute point estimates in the replicate sample
      linmed_inner(
        data = boot_data,
        D = D,
        M = M,
        Y = Y,
        C = C,
        d = d, 
        dstar = dstar,
        m = m,
        interaction_DM = interaction_DM,
        interaction_DC = interaction_DC,
        interaction_MC = interaction_MC,
        weights_name = weights_name,
        minimal = TRUE
      )
    }
    
    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster, 
        varlist = c("linmed_inner", "demean"),
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
      boot_res <- foreach::foreach(i = 1:boot_reps, .combine = comb_list_vec) %dopar% {
        boot_fnc()
      }
      boot_ATE <- boot_res$ATE
      boot_NDE <- boot_res$NDE
      boot_NIE <- boot_res$NIE
      boot_CDE <- boot_res$CDE
    }
    else {
      boot_ATE <- rep(NA_real_, boot_reps)
      boot_NDE <- rep(NA_real_, boot_reps)
      boot_NIE <- rep(NA_real_, boot_reps)
      boot_CDE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_iter <- boot_fnc()
        boot_ATE[i] <- boot_iter$ATE
        boot_NDE[i] <- boot_iter$NDE
        boot_NIE[i] <- boot_iter$NIE
        boot_CDE[i] <- boot_iter$CDE
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
    ci_ATE <- boot_ci(boot_ATE)
    ci_NDE <- boot_ci(boot_NDE)
    ci_NIE <- boot_ci(boot_NIE)
    ci_CDE <- boot_ci(boot_CDE)
    
    # compute two-tailed bootstrap p-values
    boot_pval <- function(x) {
      2 * min(
        mean(x < 0),
        mean(x > 0)
      )
    }
    pvalue_ATE <- boot_pval(boot_ATE)
    pvalue_NDE <- boot_pval(boot_NDE)
    pvalue_NIE <- boot_pval(boot_NIE)
    pvalue_CDE <- boot_pval(boot_CDE)
    
    # compile bootstrap results
    boot_out <- list(
      ci_ATE = ci_ATE,
      ci_NDE = ci_NDE,
      ci_NIE = ci_NIE,
      ci_CDE = ci_CDE,
      pvalue_ATE = pvalue_ATE,
      pvalue_NDE = pvalue_NDE,
      pvalue_NIE = pvalue_NIE,
      pvalue_CDE = pvalue_CDE,
      boot_ATE = boot_ATE,
      boot_NDE = boot_NDE,
      boot_NIE = boot_NIE,
      boot_CDE = boot_CDE
    )
  }
  
  
  # final output
  out <- est
  if (boot) {
    out <- append(out, boot_out)
  }
  return(out)
}

