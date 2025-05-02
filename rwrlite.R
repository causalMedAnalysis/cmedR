#' Regression-with-residuals estimator for interventional effects: inner 
#' function
#' 
#' @description
#' Internal function used within `rwrlite()`. See the `rwrlite()` function 
#' documentation for a description of shared function arguments. Here, we will 
#' only document the one argument that is not shared by `rwrlite_inner()` and 
#' `rwrlite()`: the `minimal` argument.
#' 
#' @param minimal A logical scalar indicating whether the function should 
#'   return only a minimal set of output. The `rwrlite()` function uses the 
#'   default of FALSE when calling `rwrlite_inner()` to generate the point 
#'   point estimates and sets the argument to TRUE when calling `rwrlite_inner()` 
#'   to perform the bootstrap.
#' 
#' @noRd
rwrlite_inner <- function(
  data,
  D,
  C = NULL,
  d = 1, 
  dstar = 0,
  m = 0,
  Y_formula,
  M_formula,
  M_family = gaussian,
  L_formula_list,
  weights = NULL,
  minimal = FALSE
) {
  # fit each of the L model formulae
  if (is.null(weights)) {
    temp_weights <- rep(1, nrow(data))
  }
  else {
    temp_weights <- weights
  }
  L_func <- function(form) {
    environment(form) <- environment()
    lm(
      formula = form,
      data = data,
      weights = temp_weights
    )
  }
  L_models <- lapply(L_formula_list, L_func)
  
  # run RWR
  rwrmed_args <- list(
    treatment = D,
    zmodels = L_models,
    y_form = Y_formula,
    m_form = M_formula,
    m_family = M_family,
    data = data
  )
  if (!is.null(C)) {
    rwrmed_args <- append(rwrmed_args, list(pre_cov = C))
  }
  if (!is.null(weights)) {
    rwrmed_args <- append(rwrmed_args, list(weights = weights))
  }
  rwr <- do.call(rwrmed::rwrmed, rwrmed_args)
  
  # causal effect decomposition
  # We will only output the two-component decomposition (plus the CDE), but 
  # rwrmed::decomp also supports a four-component decomposition.
  # Also, we will not use the bootstrap feature from rwrmed::decomp, which 
  # assumes a normal distribution. Instead, we will construct a percentile 
  # bootstrap in the rwrlite function.
  dec <- rwrmed::decomp(
    rwr,
    a0 = dstar,
    a1 = d,
    m = m,
    bootstrap = FALSE
  )
  
  # compile and output
  if (minimal) {
    out <- list(
      OE = dec$twocomp["rATE", "Estimate"],
      IDE = dec$twocomp["rNDE", "Estimate"],
      IIE = dec$twocomp["rNIE", "Estimate"],
      CDE = dec$fourcomp["CDE", "Estimate"]
    )
  }
  else {
    out <- list(
      OE = dec$twocomp["rATE", "Estimate"],
      IDE = dec$twocomp["rNDE", "Estimate"],
      IIE = dec$twocomp["rNIE", "Estimate"],
      CDE = dec$fourcomp["CDE", "Estimate"],
      models_L = L_models,
      model_M = rwr$m_model,
      model_Y = rwr$y_model,
      data_ed = rwr$data_ed
    )
  }
  return(out)
}






#' Regression-with-residuals estimator for interventional effects
#' 
#' @description
#' `rwrlite()` is a wrapper for two functions from the `rwrmed` package. It 
#' implements the regression-with-residuals (RWR) estimator for interventional 
#' effects, producing estimates of the overall effect (OE), interventional 
#' direct effect (IDE), interventional indirect effect (IIE), and controlled 
#' direct effect (CDE). You may install the `rwrmed` package from the 
#' xiangzhou09/rwrmed GitHub repository, using the following code:
#' `devtools::install_github("xiangzhou09/rwrmed")`
#' 
#' @details
#' `rwrlite()` estimates interventional and controlled direct effects as follows: (i) it 
#' fits a model for the mediator conditional on treatment and the baseline covariates, after
#' centering these covariates around their sample means, and (ii) it fits a model for the 
#' outcome conditional on treatment, the mediator, the baseline covariates after centering them 
#' around their sample means, and any exposure-induced covariates after residualizing them with 
#' respect to the treatment and baseline covariates. These models then allow for estimation of 
#' controlled direct effects, interventional direct effects, interventional indirect effects, 
#' and the overall effect, using simple functions of their parameters. Inferential statistics 
#' are computed with the nonparametric bootstrap.
#' 
#' @param data A data frame.
#' @param D A character scalar identifying the name of the exposure variable in 
#'   `data`. `D` is a character string, but the exposure variable it identifies 
#'   must be numeric.
#' @param C A character vector (of one or more elements) identifying the names 
#'   of the covariate variables in `data` that you wish to include in both the 
#'   mediator and outcome models. If there are no such covariates you wish to 
#'   include, leave `C` as its default null argument.
#' @param d,dstar A pair of arguments, each a numeric scalar denoting a specific 
#'   value of the exposure `D`. The exposure contrast of interest is 
#'   `d - dstar`.
#' @param m A numeric scalar denoting a specific value to set the mediator to, 
#'   for estimating the CDE.
#' @param Y_formula A formula object for the outcome model.
#' @param M_formula A formula object for the mediator model.
#' @param M_family The family for the mediator model, to be supplied to the 
#'   `glm` function. The family describes the error distribution and link 
#'   function to be used in the mediator model. As the `glm` documentation 
#'   describes, the family can be either a character string naming a family 
#'   function, a family function, or the result of a call to a family function.
#' @param L_formula_list A list of formula objects, one for each exposure-
#'   induced confounder model.
#' @param weights An optional numeric vector of weights to be used in fitting 
#'   each of the models (the exposure-induced confounder model(s), the mediator 
#'   model, and the outcome model).
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
#' \item{OE}{A numeric scalar with the estimated overall effect for the exposure 
#'   contrast `d - dstar`: OE(`d`,`dstar`).}
#' \item{IDE}{A numeric scalar with the estimated interventional direct effect 
#'   for the exposure contrast `d - dstar`: IDE(`d`,`dstar`).}
#' \item{IIE}{A numeric scalar with the estimated interventional indirect effect 
#'   for the exposure contrast `d - dstar`: IIE(`d`,`dstar`).}
#' \item{CDE}{A numeric scalar with the estimated controlled direct effect for 
#'   the exposure contrast `d - dstar` and the mediator value `m`: 
#'   CDE(`d`,`dstar`,`m`).}
#' \item{models_L}{A list with the model objects from each of the fitted 
#'   exposure-induced confounder models.}
#' \item{model_M}{The model object from the fitted mediator model.}
#' \item{model_Y}{The model object from the fitted outcome model.}
#' \item{data_ed}{A version of the data frame supplied by the user in the `data` 
#'   argument, but with mean-centered covariates and residualized exposure-
#'   induced confounders.}
#' 
#' If you request the bootstrap (by setting the `boot` argument to TRUE), then 
#' the function returns all of the elements listed above, as well as the 
#' following additional elements:
#' \item{ci_OE}{A numeric vector with the bootstrap confidence interval for the 
#'   overall effect (OE).}
#' \item{ci_IDE}{A numeric vector with the bootstrap confidence interval for the 
#'   interventional direct effect (IDE).}
#' \item{ci_IIE}{A numeric vector with the bootstrap confidence interval for the 
#'   interventional indirect effect (IIE).}
#' \item{ci_CDE}{A numeric vector with the bootstrap confidence interval for the 
#'   controlled direct effect (CDE).}
#' \item{pvalue_OE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the OE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_IDE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the IDE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_IIE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the IIE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_CDE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the CDE is different from zero, as computed from the bootstrap.}
#' \item{boot_OE}{A numeric vector of length `boot_reps` comprising the OE 
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_IDE}{A numeric vector of length `boot_reps` comprising the IDE 
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_IIE}{A numeric vector of length `boot_reps` comprising the IIE 
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_CDE}{A numeric vector of length `boot_reps` comprising the CDE 
#'   estimates from all replicate samples created in the bootstrap.}
#' 
#' @export
#' 
#' @examples
#' # Example 1: With one exposure-induced confounder
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
#' key_variables <- c(
#'   "cesd_age40",
#'   "ever_unemp_age3539",
#'   "log_faminc_adj_age3539",
#'   "att22",
#'   covariates
#' )
#' nlsy <- nlsy[complete.cases(nlsy[,key_variables]),]
#' nlsy$std_cesd_age40 <- 
#'   (nlsy$cesd_age40 - mean(nlsy$cesd_age40)) / 
#'   sd(nlsy$cesd_age40)
#' ## Define model formulae
#' formula_L <- ever_unemp_age3539 ~ att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3
#' formula_M <- log_faminc_adj_age3539 ~ att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3
#' formula_Y <- std_cesd_age40 ~ att22*log_faminc_adj_age3539 + ever_unemp_age3539 + 
#'   female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3
#' ## Estimate interventional effects
#' out1 <- rwrlite(
#'   data = nlsy,
#'   D = "att22",
#'   C = covariates,
#'   m = log(5e4), # evaluates the CDE at log_faminc_adj_age3539 = log(5e4)
#'   Y_formula = formula_Y,
#'   M_formula = formula_M,
#'   L_formula_list = list(formula_L)
#' )
#' head(out1,4)
#' 
#' # Example 2: Incorporating sampling weights
#' out2 <- rwrlite(
#'   data = nlsy,
#'   D = "att22",
#'   C = covariates,
#'   m = log(5e4),
#'   Y_formula = formula_Y,
#'   M_formula = formula_M,
#'   L_formula_list = list(formula_L),
#'   weights = nlsy$weight
#' )
#' head(out2,4)
#' 
#' # Example 3: Perform a nonparametric bootstrap, with 2,000 replications
#' \dontrun{
#'   out3 <- rwrlite(
#'     data = nlsy,
#'     D = "att22",
#'     C = covariates,
#'     m = log(5e4),
#'     Y_formula = formula_Y,
#'     M_formula = formula_M,
#'     L_formula_list = list(formula_L),
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234
#'   )
#'   out3[c(
#'     "OE",
#'     "IDE",
#'     "IIE",
#'     "ci_OE",
#'     "ci_IDE",
#'     "ci_IIE",
#'     "pvalue_OE",
#'     "pvalue_IDE",
#'     "pvalue_IIE"
#'   )]
#' }
#' 
#' # Example 4: Parallelize the bootstrap, to attempt to reduce runtime
#' # Note that this requires you to have installed the `doParallel`, `doRNG`, 
#' # and `foreach` packages.
#' \dontrun{
#'   out4 <- rwrlite(
#'     data = nlsy,
#'     D = "att22",
#'     C = covariates,
#'     m = log(5e4),
#'     Y_formula = formula_Y,
#'     M_formula = formula_M,
#'     L_formula_list = list(formula_L),
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_seed = 1234,
#'     boot_parallel = TRUE
#'   )
#'   out4[c(
#'     "OE",
#'     "IDE",
#'     "IIE",
#'     "ci_OE",
#'     "ci_IDE",
#'     "ci_IIE",
#'     "pvalue_OE",
#'     "pvalue_IDE",
#'     "pvalue_IIE"
#'   )]
#' }
rwrlite <- function(
  data,
  D,
  C = NULL,
  d = 1, 
  dstar = 0,
  m = 0,
  Y_formula,
  M_formula,
  M_family = gaussian,
  L_formula_list,
  weights = NULL,
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
  if (!requireNamespace("rwrmed", quietly = TRUE)) {
    stop(paste(strwrap("Error: The required package 'rwrmed' has not been installed. Please install this package from the xiangzhou09/rwrmed GitHub repository."), collapse = "\n"))
  }
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
    if (!is.null(weights)) {
      warning(paste(strwrap("Warning: You requested a bootstrap, but your design includes sampling weights. Note that this function does not internally rescale sampling weights for use with the bootstrap, and it does not account for any stratification or clustering in your sample design. Failure to properly adjust the bootstrap sampling to account for a complex sample design that requires weighting could lead to invalid inferential statistics."), collapse = "\n"))
    }
  }
  
  
  # compute point estimates
  est <- rwrlite_inner(
    data = data_outer,
    D = D,
    C = C,
    d = d, 
    dstar = dstar,
    m = m,
    Y_formula = Y_formula,
    M_formula = M_formula,
    M_family = M_family,
    L_formula_list = L_formula_list,
    weights = weights,
    minimal = FALSE
  )
  
  
  # bootstrap, if requested
  if (boot) {
    # bootstrap function
    boot_fnc <- function() {
      # sample from the data with replacement
      boot_data <- data_outer[sample(nrow(data_outer), size = nrow(data_outer), replace = TRUE), ]
      
      # compute point estimates in the replicate sample
      rwrlite_inner(
        data = boot_data,
        D = D,
        C = C,
        d = d, 
        dstar = dstar,
        m = m,
        Y_formula = Y_formula,
        M_formula = M_formula,
        M_family = M_family,
        L_formula_list = L_formula_list,
        weights = weights,
        minimal = TRUE
      )
    }
    
    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster, 
        varlist = c("rwrlite_inner", "demean"),
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
      boot_OE <- boot_res$OE
      boot_IDE <- boot_res$IDE
      boot_IIE <- boot_res$IIE
      boot_CDE <- boot_res$CDE
    }
    else {
      boot_OE <- rep(NA_real_, boot_reps)
      boot_IDE <- rep(NA_real_, boot_reps)
      boot_IIE <- rep(NA_real_, boot_reps)
      boot_CDE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_iter <- boot_fnc()
        boot_OE[i] <- boot_iter$OE
        boot_IDE[i] <- boot_iter$IDE
        boot_IIE[i] <- boot_iter$IIE
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
    ci_OE <- boot_ci(boot_OE)
    ci_IDE <- boot_ci(boot_IDE)
    ci_IIE <- boot_ci(boot_IIE)
    ci_CDE <- boot_ci(boot_CDE)
    
    # compute two-tailed bootstrap p-values
    boot_pval <- function(x) {
      2 * min(
        mean(x < 0),
        mean(x > 0)
      )
    }
    pvalue_OE <- boot_pval(boot_OE)
    pvalue_IDE <- boot_pval(boot_IDE)
    pvalue_IIE <- boot_pval(boot_IIE)
    pvalue_CDE <- boot_pval(boot_CDE)
    
    # compile bootstrap results
    boot_out <- list(
      ci_OE = ci_OE,
      ci_IDE = ci_IDE,
      ci_IIE = ci_IIE,
      ci_CDE = ci_CDE,
      pvalue_OE = pvalue_OE,
      pvalue_IDE = pvalue_IDE,
      pvalue_IIE = pvalue_IIE,
      pvalue_CDE = pvalue_CDE,
      boot_OE = boot_OE,
      boot_IDE = boot_IDE,
      boot_IIE = boot_IIE,
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

