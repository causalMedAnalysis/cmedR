#' Inverse probability weighting (IPW) estimator for natural effects: inner
#' function
#'
#' @description
#' Internal function used within `ipwmed()`. See the `ipwmed()` function
#' documentation for a description of shared function arguments. Here, we will
#' only document the one argument that is not shared by `ipwmed_inner()` and
#' `ipwmed()`: the `minimal` argument.
#'
#' @param minimal A logical scalar indicating whether the function should
#'   return only a minimal set of outputs. The `ipwmed()` function uses the
#'   default of FALSE when calling `ipwmed_inner()` to generate the point
#'   estimates and sets the argument to TRUE when calling `ipwmed_inner()`
#'   to perform the bootstrap.
#'
#' @noRd
ipwmed_inner <- function(
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

  # fit exposure models
  d_model1 <- glm(
    as.formula(formula1_string),
    data = df,
    family = quasibinomial(link = "logit"),
    weights = base_weights_rsc
  )
  d_model2 <- glm(
    as.formula(formula2_string),
    data = df,
    family = quasibinomial(link = "logit"),
    weights = base_weights_rsc
  )
  # enforcing a no-missing-data rule
  if (nrow(d_model1$model)!=nrow(df) | nrow(d_model2$model)!=nrow(df)) {
    stop(paste(strwrap("Error: Please remove observations with missing/NA values for the exposure, mediator, outcome, or covariates."), collapse = "\n"))
  }

  # predict exposure probabilities
  ps1_d <- predict(d_model1, newdata = df, type = "response")
  ps2_d <- predict(d_model2, newdata = df, type = "response")
  ps1_dstar <- 1 - ps1_d
  ps2_dstar <- 1 - ps2_d
  marg_prob_d <- weighted.mean(ps1_d, base_weights_rsc)
  marg_prob_dstar <- 1 - marg_prob_d

  # for convenience, create logical vectors to identify subgroups
  group_dstar <- df[[D]]==dstar & base_weights_rsc>0
  group_d <- df[[D]]==d & base_weights_rsc>0

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
  Ehat_Ydstar_Mdstar <- weighted.mean(df[[Y]], final_w1)
  Ehat_Yd_Md <- weighted.mean(df[[Y]], final_w2)
  Ehat_Yd_Mdstar <- weighted.mean(df[[Y]], final_w3)
  ATE <- Ehat_Yd_Md - Ehat_Ydstar_Mdstar
  NDE <- Ehat_Yd_Mdstar - Ehat_Ydstar_Mdstar
  NIE <- Ehat_Yd_Md - Ehat_Yd_Mdstar

  # compile and output
  if (minimal) {
    out <- list(
      ATE = ATE,
      NDE = NDE,
      NIE = NIE
    )
  }
  else {
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
  }
  return(out)
}

#' Inverse probability weighting (IPW) estimator for natural effects
#'
#' @description
#' `ipwmed()` uses the inverse probability weighting (IPW) estimator to estimate
#' the total effect (ATE), natural direct effect (NDE), and natural indirect
#' effect (NIE). The function supports estimation of both univariate and
#' multivariate natural effects.
#'
#' @details
#' `ipwmed()` performs causal mediation analysis using inverse probability weighting,
#' and it computes inferential statistics using the nonparametric bootstrap. It estimates
#' two models to construct the weights: a logit model for the exposure conditional on
#' baseline covariates (if specified), and another logit model for the exposure conditional
#' on the mediator(s) and the baseline covariates. Using these weights, `ipwmed()` estimates the
#' total, natural direct, and natural indirect effects when a single mediator is specified.
#' When multiple mediators are specified, it provides estimates for the total effect and
#' the multivariate natural direct and indirect effects operating through the entire set of
#' mediators.
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
#'   denoting a probability with values in \eqn{[0, 1]}. If the `censor` argument is
#'   TRUE, then IPW weights below the `censor_low` quantile will be
#'   bottom-coded, and IPW weights above the `censor_high` quantile will be
#'   top-coded (before multiplying by a rescaled version of the base weights, if
#'   applicable). E.g., if the default options of `censor_low = 0.01` and
#'   `censor_high = 0.99` are used, then the IPW weights will be censored at
#'   their 1st and 99th percentiles in the data.
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
#' @returns By default, `ipwmed()` returns a list with the following elements:
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
#' If you request the bootstrap (by setting the `boot` argument to TRUE), then
#' the function returns all of the elements listed above, as well as the
#' following additional elements:
#' \item{ci_ATE}{A numeric vector with the bootstrap confidence interval for the
#'   total average treatment effect (ATE).}
#' \item{ci_NDE}{A numeric vector with the bootstrap confidence interval for the
#'   natural direct effect (NDE).}
#' \item{ci_NIE}{A numeric vector with the bootstrap confidence interval for the
#'   natural indirect effect (NIE).}
#' \item{pvalue_ATE}{A numeric scalar with the p-value from a two-sided test of
#'   whether the ATE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_NDE}{A numeric scalar with the p-value from a two-sided test of
#'   whether the NDE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_NIE}{A numeric scalar with the p-value from a two-sided test of
#'   whether the NIE is different from zero, as computed from the bootstrap.}
#' \item{boot_ATE}{A numeric vector of length `boot_reps` comprising the ATE
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_NDE}{A numeric vector of length `boot_reps` comprising the NDE
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_NIE}{A numeric vector of length `boot_reps` comprising the NIE
#'   estimates from all replicate samples created in the bootstrap.}
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
#'   formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3",
#'   formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3+ever_unemp_age3539"
#' )
#' head(out1,3)
#'
#' # Example 2: Incorporating sampling weights
#' out2 <- ipwmed(
#'   data = nlsy1,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3",
#'   formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3+ever_unemp_age3539",
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
#'   formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3",
#'   formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'   famsize+afqt3+ever_unemp_age3539+log_faminc_adj_age3539"
#' )
#' head(out3,3)
#'
#' # Example 4: Perform a nonparametric bootstrap, with 2,000 replications
#' \dontrun{
#'   out4 <- ipwmed(
#'     data = nlsy1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     Y = "std_cesd_age40",
#'     formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'     famsize+afqt3",
#'     formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'     famsize+afqt3+ever_unemp_age3539",
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234
#'   )
#'   out4[c(
#'     "ATE",
#'     "NDE",
#'     "NIE",
#'     "ci_ATE",
#'     "ci_NDE",
#'     "ci_NIE",
#'     "pvalue_ATE",
#'     "pvalue_NDE",
#'     "pvalue_NIE"
#'   )]
#' }
#'
#' # Example 5: Parallelize the bootstrap, to attempt to reduce runtime
#' \dontrun{
#'   out5 <- ipwmed(
#'     data = nlsy1,
#'     D = "att22",
#'     M = "ever_unemp_age3539",
#'     Y = "std_cesd_age40",
#'     formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'     famsize+afqt3",
#'     formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+
#'     famsize+afqt3+ever_unemp_age3539",
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234,
#'     boot_parallel = TRUE
#'   )
#'   out5[c(
#'     "ATE",
#'     "NDE",
#'     "NIE",
#'     "ci_ATE",
#'     "ci_NDE",
#'     "ci_NIE",
#'     "pvalue_ATE",
#'     "pvalue_NDE",
#'     "pvalue_NIE"
#'   )]
#' }

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
  if (!is.numeric(data_outer[[D]])) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be numeric."), collapse = "\n"))
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

  # compute point estimates
  est <- ipwmed_inner(
    data = data_outer,
    D = D,
    M = M,
    Y = Y,
    formula1_string = formula1_string,
    formula2_string = formula2_string,
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
      ipwmed_inner(
        data = boot_data,
        D = D,
        M = M,
        Y = Y,
        formula1_string = formula1_string,
        formula2_string = formula2_string,
        base_weights_name = base_weights_name,
        stabilize = stabilize,
        censor = censor,
        censor_low = censor_low,
        censor_high = censor_high,
        minimal = TRUE
      )
    }

    # parallelization prep, if parallelization is requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster,
        varlist = c("ipwmed_inner", "trimQ"),
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
    }
    else {
      boot_ATE <- rep(NA_real_, boot_reps)
      boot_NDE <- rep(NA_real_, boot_reps)
      boot_NIE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_iter <- boot_fnc()
        boot_ATE[i] <- boot_iter$ATE
        boot_NDE[i] <- boot_iter$NDE
        boot_NIE[i] <- boot_iter$NIE
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

    # compile bootstrap results
    boot_out <- list(
      ci_ATE = ci_ATE,
      ci_NDE = ci_NDE,
      ci_NIE = ci_NIE,
      pvalue_ATE = pvalue_ATE,
      pvalue_NDE = pvalue_NDE,
      pvalue_NIE = pvalue_NIE,
      boot_ATE = boot_ATE,
      boot_NDE = boot_NDE,
      boot_NIE = boot_NIE
    )
  }

  # final output
  out <- est
  if (boot) {
    out <- append(out, boot_out)
  }
  return(out)
}
