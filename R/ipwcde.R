#' Controlled direct effect inverse probability weighting (IPW) estimator: inner
#' function
#'
#' @description
#' Internal function used within `ipwcde()`. See the `ipwcde()` function
#' documentation for a description of shared function arguments. Here, we will
#' only document the one argument that is not shared by `ipwcde_inner()` and
#' `ipwcde()`: the `minimal` argument.
#'
#' @param minimal A logical scalar indicating whether the function should
#'   return only a minimal set of output. The `ipwcde()` function uses the
#'   default of FALSE when calling `ipwcde_inner()` to generate the point
#'   estimates and sets the argument to TRUE when calling `ipwcde_inner()`
#'   to perform the bootstrap.
#'
#' @noRd
ipwcde_inner <- function(
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
  censor_high = 0.99,
  minimal = FALSE
) {
  # preliminaries
  d <- 1
  dstar <- 0

  # load data
  df <- data

  # assign base weights
  if (is.null(base_weights_name)) {
    base_weights <- rep(1, nrow(df))
  }
  else {
    base_weights <- df[[base_weights_name]]
  }
  if (!is.null(base_weights_name) & any(is.na(base_weights))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the base weights variable (identified by the string argument base_weights_name in data). If that observation should not receive a positive weight, please replace the NA value with a zero before proceeding."), collapse = "\n"))
  }

  # rescale base weights
  base_weights_rsc <- base_weights / mean(base_weights)

  # fit specified models
  d_model <- glm(
    as.formula(formula_D_string),
    data = df,
    family = quasibinomial(link = "logit"),
    weights = base_weights_rsc
  )
  m_model <- glm(
    as.formula(formula_M_string),
    data = df,
    family = quasibinomial(link = "logit"),
    weights = base_weights_rsc
  )
  # enforcing a no-missing-data rule
  if (nrow(d_model$model)!=nrow(df) | nrow(m_model$model)!=nrow(df)) {
    stop(paste(strwrap("Error: Please remove observations with missing/NA values for the exposure, mediator, outcome, or covariates."), collapse = "\n"))
  }

  # additionally fit mediator model without covariates
  m_model_no_cov <- glm(
    as.formula(paste0(M,"~",D)),
    data = df,
    family = quasibinomial(link = "logit"),
    weights = base_weights_rsc
  )

  # predict exposure and mediator probabilities
  ps_D1_C <- predict(d_model, newdata = df, type = "response")
  ps_M1_CD <- predict(m_model, newdata = df, type = "response")
  ps_M1_D <- predict(m_model_no_cov, newdata = df, type = "response")
  ps_D_C <- ifelse(as.logical(df[[D]]), ps_D1_C, 1-ps_D1_C)
  ps_M_CD <- ifelse(as.logical(df[[M]]), ps_M1_CD, 1-ps_M1_CD)
  marg_prob_D1 <- weighted.mean(df[[D]], base_weights_rsc)
  marg_prob_D <- ifelse(
    as.logical(df[[D]]),
    marg_prob_D1,
    1 - marg_prob_D1
  )
  marg_prob_M_D <- ifelse(
    as.logical(df[[M]]),
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
  group_d_m <- df[[D]]==d & df[[M]]==m
  group_dstar_m <- df[[D]]==dstar & df[[M]]==m
  CDE <-
    weighted.mean(df[[Y]][group_d_m], final_w4[group_d_m]) -
    weighted.mean(df[[Y]][group_dstar_m], final_w4[group_dstar_m])

  # compile and output
  if (minimal) {
    out <- CDE
  }
  else {
    out <- list(
      CDE = CDE,
      weights = final_w4,
      model_d = d_model,
      model_m = m_model
    )
  }
  return(out)
}


#' Controlled direct effect inverse probability weighting (IPW) estimator
#'
#' @description
#' `ipwcde()` uses the inverse probability weighting (IPW) estimator to estimate
#' the controlled direct effect (CDE). Note that unlike the `ipwmed()` function,
#' `ipwcde()` requires a single mediator. Multiple mediators are not supported.
#'
#' @details
#' `ipwcde()` estimates controlled direct effects using inverse probability weights,
#' and it computes inferential statistics using the nonparametric bootstrap. To compute
#' the weights, `ipwcde()` fits the following models: (i) a logit model for exposure
#' conditional on baseline covariates, and (ii) a logit model for the mediator
#' conditional on both the exposure and baseline covariates, plus any specified
#' post-treatment covariates. These models are used to generate inverse probability
#' weights, which are then applied to fit an outcome model and estimate the controlled
#' direct effect.
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
#'   denoting a probability with values in \eqn{[0, 1]}. If the `censor` argument is
#'   TRUE, then IPW weights below the `censor_low` quantile will be
#'   bottom-coded, and IPW weights above the `censor_high` quantile will be
#'   top-coded (before multiplying by a rescaled version of the base weights, if
#'   applicable). E.g., if the default options of `censor_low = 0.01` and
#'   `censor_high = 0.99` are used, then the IPW weights will be censored at
#'   their 1st and 99th percentiles in the data.
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
#' @returns By default, `ipwcde()` returns a list with the following elements:
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
#' If you request the bootstrap (by setting the `boot` argument to TRUE), then
#' the function returns all of the elements listed above, as well as the
#' following additional elements:
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
#'   formula_D_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3",
#'   formula_M_string = "ever_unemp_age3539~att22+female+black+hispan+paredu+
#'   parprof+parinc_prank+famsize+afqt3"
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
#'   formula_D_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3",
#'   formula_M_string = "ever_unemp_age3539~att22+female+black+hispan+paredu+
#'   parprof+parinc_prank+famsize+afqt3",
#'   base_weights_name = "weight"
#' )
#' head(out2,1)
#'
#' # Example 3: Perform a nonparametric bootstrap, with 2,000 replications
#' \dontrun{
#'   out3 <- ipwcde(
#'     data = nlsy1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     Y = "std_cesd_age40",
#'     m = 1,
#'     formula_D_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'     famsize+afqt3",
#'     formula_M_string = "ever_unemp_age3539~att22+female+black+hispan+paredu+
#'     parprof+parinc_prank+famsize+afqt3",
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234
#'   )
#'   out3[c(
#'     "CDE",
#'     "ci_CDE",
#'     "pvalue_CDE"
#'   )]
#' }
#'
#' # Example 4: Parallelize the bootstrap, to attempt to reduce runtime
#' \dontrun{
#'   out4 <- ipwcde(
#'     data = nlsy1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     Y = "std_cesd_age40",
#'     m = 1,
#'     formula_D_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'     famsize+afqt3",
#'     formula_M_string = "ever_unemp_age3539~att22+female+black+hispan+paredu+
#'     parprof+parinc_prank+famsize+afqt3",
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234,
#'     boot_parallel = TRUE
#'   )
#'   out4[c(
#'     "CDE",
#'     "ci_CDE",
#'     "pvalue_CDE"
#'   )]
#' }
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
  censor_high = 0.99,
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

  # preliminary error/warning checks for the bootstrap
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
    if (!is.null(base_weights_name)) {
      warning(paste(strwrap("Warning: You requested a bootstrap, but your design includes base sampling weights. Note that this function does not internally rescale sampling weights for use with the bootstrap, and it does not account for any stratification or clustering in your sample design. Failure to properly adjust the bootstrap sampling to account for a complex sample design that requires weighting could lead to invalid inferential statistics."), collapse = "\n"))
    }
  }

  # other error/warning checks
  if (length(M)>1) {
    stop(paste(strwrap("Error: Unlike the ipwmed() function, ipwcde() requires a single mediator. Multiple mediators are not supported."), collapse = "\n"))
  }
  if (length(m)>1) {
    stop(paste(strwrap("Error: Please specify a single value for the argument m."), collapse = "\n"))
  }
  if (!is.numeric(data_outer[[D]])) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be numeric."), collapse = "\n"))
  }
  if (!is.numeric(data_outer[[M]])) {
    stop(paste(strwrap("Error: The mediator variable (identified by the string argument M in data) must be numeric."), collapse = "\n"))
  }
  if (!is.numeric(data_outer[[Y]])) {
    stop(paste(strwrap("Error: The outcome variable (identified by the string argument Y in data) must be numeric."), collapse = "\n"))
  }
  if (any(is.na(data_outer[[D]]))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the exposure variable (identified by the string argument D in data)."), collapse = "\n"))
  }
  if (any(! data_outer[[D]] %in% c(0,1))) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be a numeric variable consisting only of the values 0 or 1. There is at least one observation in the data that does not meet this criterion."), collapse = "\n"))
  }
  if (any(is.na(data_outer[[M]]))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the mediator variable (identified by the string argument M in data)."), collapse = "\n"))
  }
  if (any(! data_outer[[M]] %in% c(0,1))) {
    stop(paste(strwrap("Error: Only binary mediator variables are supported. The mediator variable (identified by the string argument M in data) must be a numeric variable consisting only of the values 0 or 1. There is at least one observation in the data that does not meet this criterion."), collapse = "\n"))
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


  # compute point estimates
  est <- ipwcde_inner(
    data = data_outer,
    D = D,
    M = M,
    Y = Y,
    m = m,
    formula_D_string = formula_D_string,
    formula_M_string = formula_M_string,
    base_weights_name = base_weights_name,
    stabilize = stabilize,
    censor = censor,
    censor_low = censor_low,
    censor_high = censor_high,
    minimal = FALSE
  )

  # bootstrap, if requested
  if (boot) {
    # bootstrap function
    boot_fnc <- function() {
      # sample from the data with replacement
      boot_data <- data_outer[sample(nrow(data_outer), size = nrow(data_outer), replace = TRUE), ]

      # compute point estimates in the replicate sample
      ipwcde_inner(
        data = boot_data,
        D = D,
        M = M,
        Y = Y,
        m = m,
        formula_D_string = formula_D_string,
        formula_M_string = formula_M_string,
        base_weights_name = base_weights_name,
        stabilize = stabilize,
        censor = censor,
        censor_low = censor_low,
        censor_high = censor_high,
        minimal = TRUE
      )
    }

    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster,
        varlist = c("ipwcde_inner", "trimQ"),
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

    # compile bootstrap results
    boot_out <- list(
      ci_CDE = ci_CDE,
      pvalue_CDE = pvalue_CDE,
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

