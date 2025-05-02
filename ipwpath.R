#' Inverse probability weighting (IPW) estimator for path-specific effects: 
#' inner function
#' 
#' @description
#' Internal function used within `ipwpath()`. See the `ipwpath()` function 
#' documentation for a description of shared function arguments. Here, we will 
#' only document the one argument that is not shared by `ipwpath_inner()` and 
#' `ipwpath()`: the `minimal` argument.
#' 
#' @param minimal A logical scalar indicating whether the function should 
#'   return only a minimal set of output. The `ipwpath()` function uses the 
#'   default of FALSE when calling `ipwpath_inner()` to generate the point 
#'   point estimates and sets the argument to TRUE when calling `ipwpath_inner()` 
#'   to perform the bootstrap.
#' 
#' @noRd
ipwpath_inner <- function(
    data,
    D,
    M,
    Y,
    C = NULL,
    base_weights_name = NULL,
    stabilize = TRUE,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99,
    minimal = FALSE
) {
  # prep
  K <- length(M) # number of mediators (treating multivariate mediators as a single mediator)
  n_PSE <- K + 1 # number of PSEs
  PSE <- rep(NA_real_, n_PSE) # vector to store path-specific effects
  if (!minimal) {
    # matrix to store weights3 for each k iteration
    # (weights1 and weights2 will be identical for each k iteration)
    weights3 <- matrix(NA_real_, nrow = nrow(data), ncol = K)
    # list to store the fitted second D models
    # (the first D models will be identical for each k iteration)
    models_d2 <- vector(mode = "list", length = K)
  }
  
  
  # build additive D model 1 formula: f(D|C)
  # (will be the same for all iterations of the loop over mediators below)
  if (is.null(C)) {
    predictors1_D <- "1"
  }
  else {
    predictors1_D <- paste(C, collapse = " + ")
  }
  formula1_D_string <- paste(D, "~", predictors1_D)
  
  
  # loop over mediators in reverse order to estimate PSEs
  for (k in rev(seq_len(K))) {
    ## PSE index
    PSE_index <- K - k + 1
    
    ## build additive D model 2 formula: s(D|C,M)
    predictors2_D <- paste(c(unlist(M[1:k]),C), collapse = " + ")
    formula2_D_string <- paste(D, "~", predictors2_D)
    
    ## estimate multivariate natural effects
    est <- ipwmed_inner(
      data = data,
      D = D,
      M = unlist(M[1:k]),
      Y = Y,
      formula1_string = formula1_D_string,
      formula2_string = formula2_D_string,
      base_weights_name = base_weights_name,
      stabilize = stabilize,
      censor = censor,
      censor_low = censor_low,
      censor_high = censor_high,
      minimal = minimal
    )
    
    ## special case: only one total mediator
    if (K==1) {
      PSE <- c(est$NDE, est$NIE)
      names(PSE) <- c("NDE", "NIE")
      ATE <- est$ATE
      if (!minimal) {
        weights1 <- est$weights1
        weights2 <- est$weights2
        weights3 <- est$weights3
        model_d1 <- est$model_d1
        models_d2 <- est$model_d2
      }
    }
    
    ## 2+ total mediators: last mediator
    else if (k==K) {
      PSE[[PSE_index]] <- est$NDE
      names(PSE)[PSE_index] <- "D->Y"
      if (!minimal) {
        weights1 <- est$weights1
        weights2 <- est$weights2
        weights3[,k] <- est$weights3
        model_d1 <- est$model_d1
        models_d2[[k]] <- est$model_d2
        colnames(weights3)[k] <- names(models_d2)[k] <- paste0("M1:M",k)
      }
      prev_MNDE <- est$NDE
    }
    
    ## 2+ total mediators: first mediator
    else if (k==1) {
      PSE[[PSE_index]] <- est$NDE - prev_MNDE
      PSE[[n_PSE]] <- est$NIE
      if (k+1==K) {
        names(PSE)[PSE_index] <- paste0("D->M",k+1,"->Y")
      }
      else {
        names(PSE)[PSE_index] <- paste0("D->M",k+1,"~>Y")
      }
      names(PSE)[n_PSE] <- "D->M1~>Y"
      ATE <- est$ATE
      if (!minimal) {
        weights3[,k] <- est$weights3
        models_d2[[k]] <- est$model_d2
        colnames(weights3)[k] <- names(models_d2)[k] <- "M1"
      }
    }
    
    ## 2+ total mediators: all other mediators
    else {
      PSE[[PSE_index]] <- est$NDE - prev_MNDE
      if (k+1==K) {
        names(PSE)[PSE_index] <- paste0("D->M",k+1,"->Y")
      }
      else {
        names(PSE)[PSE_index] <- paste0("D->M",k+1,"~>Y")
      }
      if (!minimal) {
        weights3[,k] <- est$weights3
        models_d2[[k]] <- est$model_d2
        colnames(weights3)[k] <- names(models_d2)[k] <- paste0("M1:M",k)
      }
      prev_MNDE <- est$NDE
    }
  }
  
  
  # compile and output
  if (minimal) {
    out <- list(
      ATE = ATE,
      PSE = PSE
    )
  }
  else {
    out <- list(
      ATE = ATE,
      PSE = PSE,
      weights1 = weights1,
      weights2 = weights2,
      weights3 = weights3,
      model_d1 = model_d1,
      models_d2 = models_d2
    )
  }
  return(out)
}






#' Inverse probability weighting (IPW) estimator for path-specific effects
#' 
#' @description
#' `ipwpath()` uses the inverse probability weighting estimator to estimate the 
#' total effect (ATE) and path-specific effects (PSEs).
#' 
#' @details
#' `ipwpath()` uses inverse probability weighting to estimate path-specific effects through 
#' multiple mediators, and it computes inferential statistics using the nonparametric bootstrap. 
#' It constructs weights using logit models for the exposure, each conditioned on different combinations 
#' of mediators and the baseline confounders. This approach allows the isolation of effects 
#' that pass through specific causal pathways defined by the ordered mediators. For example, 
#' with two causally ordered mediators (M1, M2), it estimates: a model for the exposure with 
#' no mediators, just the baseline confounders as predictors; a model for the exposure with M1 
#' and the baseline confounders as predictors; and a model for the exposure with M1, M2, and the 
#' baseline confounders as predictors. These models are used to construct a set of inverse probability 
#' weights, which then enable the calculation of direct and path-specific effects that operate through 
#' each mediator, net of the mediators that precede it causal order.
#' 
#' Specifying the `M` Argument:
#' 
#' The `M` argument is a list of character vectors identifying the names of the 
#' mediator variables. This argument is purposely a list of vectors rather than 
#' simply a vector because it accommodates both univariate and multivariate 
#' mediators. To explain, let's start with a simple example.
#' 
#' Suppose you have two single mediators, named `ever_unemp_age3539` and 
#' `log_faminc_adj_age3539`, where `ever_unemp_age3539` causally precedes 
#' `log_faminc_adj_age3539`. In this case, you would use the following syntax: 
#' `M = list("ever_unemp_age3539", "log_faminc_adj_age3539")`.
#' 
#' Now, let's say you have a third mediator, named `m3`. You believe that 
#' `ever_unemp_age3539` causally precedes both `log_faminc_adj_age3539` and 
#' `m3`. But you are unwilling to make an assumption about the relative causal 
#' order of `log_faminc_adj_age3539` and `m3` (whether `log_faminc_adj_age3539` 
#' causally precedes `m3` or vice versa). In that case, you could treat 
#' `log_faminc_adj_age3539` and `m3` as a whole, using the following syntax: 
#' `M = list("ever_unemp_age3539", c("log_faminc_adj_age3539", "m3"))`.
#' 
#' Note that the order of the elements in the `c("log_faminc_adj_age3539", "m3")` 
#' vector does not matter (it could alternatively be written as 
#' `c("m3", "log_faminc_adj_age3539")`). But the order of the vectors in the 
#' list does matter. And in this example, the mediator identified by the first 
#' element in the list, the `"ever_unemp_age3539"` scalar, is assumed to 
#' causally precede the two mediators collectively identified by the second 
#' element in the list, the `c("log_faminc_adj_age3539", "m3")` vector.
#' 
#' Finally, note that if one of your mediators is a nominal factor variable, we 
#' recommend that you dummy-encode the levels of the factor and treat the dummy 
#' variables as a multivariate whole. For instance, let's say that you have a 
#' fourth mediator, which causally follows `ever_unemp_age3539`, 
#' `log_faminc_adj_age3539`, and `m3`. This fourth mediator is a nominal 
#' variable with four levels. If you create numeric dummy variables for three 
#' levels (omitting a reference level), named `level2`, `level3`, `level4`, 
#' then you can use the following syntax for the `M` argument:
#' `M = list("ever_unemp_age3539", c("log_faminc_adj_age3539", "m3"), 
#'  c("level2","level3","level4"))`.
#' 
#' @param data A data frame.
#' @param D A character scalar identifying the name of the exposure variable in 
#'   `data`. `D` is a character string, but the exposure variable it identifies 
#'   must be numeric.
#' @param M A character vector (of one or more elements) identifying the names 
#'   of the mediator variables in `data`. The character vector MUST specify the 
#'   mediators in causal order, starting from the first in the hypothesized 
#'   causal sequence to the last. If you only specify a single mediator 
#'   variable, then the function will simply return the natural effects. Also 
#'   note that `M` is a character vector, but the mediator variable(s) it 
#'   identifies must each be numeric.
#' @param Y A character scalar identifying the name of the outcome variable in 
#'   `data`. `Y` is a character string, but the outcome variable it identifies 
#'   must be numeric.
#' @param C A character vector (of one or more elements) identifying the names 
#'   of the covariate variables in `data` that you wish to include in the 
#'   exposure models. If there are no such covariates you wish to include, leave 
#'   `C` as its default null argument.
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
#' @returns By default, `ipwpath()` returns a list with the following elements:
#' \item{ATE}{A numeric scalar with the estimated total average treatment effect 
#'   for the exposure contrast `d - dstar`: ATE(`d`,`dstar`).}
#' \item{PSE}{A numeric vector, of length `length(M)+1`, with the estimated 
#'   path-specific effects for the exposure contrast `d - dstar`. The vector is 
#'   named with the path each effect describes.}
#' \item{weights1}{A numeric vector with the final inverse probability weights 
#'   for w_1 (as defined in the book), used to estimate the multivariate natural 
#'   effects for each set of mediators, M_1 to M_k.}
#' \item{weights2}{A numeric vector with the final inverse probability weights 
#'   for w_2 (as defined in the book), used to estimate the multivariate natural 
#'   effects for each set of mediators, M_1 to M_k.}
#' \item{weights3}{A numeric matrix with the final inverse probability weights 
#'   for w_3 (as defined in the book), used to estimate the multivariate natural 
#'   effects for each set of mediators, M_1 to M_k. Unlike w_1 and w_2, the w_3 
#'   weights are unique for each set of mediators. Therefore, the matrix 
#'   consists of `length(M)` columns, with each column containing the w_3 
#'   weights for mediators M_1 to M_k.}
#' \item{model_d1}{The model object from the first fitted exposure model 
#'   (of the exposure given baseline covariates, denoted in the book as f(D|C)), 
#'   used in the estimation of the multivariate natural effects for each set of 
#'   mediators, M_1 to M_k. Since the predictors in the model do not include the 
#'   mediators, the same model is used for each set of mediators, and there is 
#'   only one model object returned.}
#' \item{models_d2}{A list of the model objects from the second fitted exposure 
#'   model (of the exposure given baseline covariates and the mediator(s), 
#'   denoted in the book as s(D|C,M)), used in the estimation of the 
#'   multivariate natural effects for each set of mediators, M_1 to M_k. Since 
#'   the predictors in the model do include the mediators, a different model is 
#'   fitted for each set of mediators, and the returned object is a list of 
#'   model objects.}
#' 
#' If you request the bootstrap (by setting the `boot` argument to TRUE), then 
#' the function returns all of the elements listed above, as well as the 
#' following additional elements:
#' \item{ci_ATE}{A numeric vector with the bootstrap confidence interval for the 
#'   total average treatment effect (ATE).}
#' \item{ci_PSE}{A numeric matrix with the bootstrap confidence interval for 
#'   each path-specific effect (PSE).}
#' \item{pvalue_ATE}{A numeric scalar with the p-value from a two-sided test of 
#'   whether the ATE is different from zero, as computed from the bootstrap.}
#' \item{pvalue_PSE}{A numeric matrix with each p-value from a two-sided test of 
#'   whether the PSE is different from zero, as computed from the bootstrap.}
#' \item{boot_ATE}{A numeric vector of length `boot_reps` comprising the ATE 
#'   estimates from all replicate samples created in the bootstrap.}
#' \item{boot_PSE}{A numeric matrix, of `length(M)+1` columns and `boot_reps` 
#'   rows, comprising all PSE estimates from all replicate samples created in 
#'   the bootstrap.}
#' 
#' @export
#' 
#' @examples
#' # Example 1: Two mediators
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
#'   "log_faminc_adj_age3539",
#'   "att22",
#'   covariates
#' )
#' nlsy <- nlsy[complete.cases(nlsy[,key_variables]),]
#' nlsy$std_cesd_age40 <- 
#'   (nlsy$cesd_age40 - mean(nlsy$cesd_age40)) / 
#'   sd(nlsy$cesd_age40)
#' ## Estimate path-specific effects
#' ipwpath(
#'   data = nlsy,
#'   D = "att22",
#'   M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
#'   # ^ note that this order encodes our assumption that ever_unemp_age3539 
#'   # causally precedes log_faminc_adj_age3539
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
#' # Example 2: If you specify only a single mediator, the function will return 
#' # the natural effects (NDE and NIE), in addition to the ATE
#' ipwpath(
#'   data = nlsy,
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
#' # Example 3: Incorporating sampling weights
#' ipwpath(
#'   data = nlsy,
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
#'   base_weights_name = "weight"
#' )
#' 
#' # Example 4: Perform a nonparametric bootstrap, with 2,000 replications
#' \dontrun{
#'   ipwpath(
#'     data = nlsy,
#'     D = "att22",
#'     M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
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
#' # Example 5: Parallelize the bootstrap, to attempt to reduce runtime
#' # Note that this requires you to have installed the `doParallel`, `doRNG`, 
#' # and `foreach` packages.
#' \dontrun{
#'   ipwpath(
#'     data = nlsy,
#'     D = "att22",
#'     M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
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
ipwpath <- function(
    data,
    D,
    M,
    Y,
    C = NULL,
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
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be a numeric variable consisting only of the values 0 or 1. There is at least one observation in the data that does not meet this criteria."), collapse = "\n"))
  }
  
  
  # compute point estimates
  est <- ipwpath_inner(
    data = data_outer,
    D = D,
    M = M,
    Y = Y,
    C = C,
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
      boot_out <- ipwpath_inner(
        data = boot_data,
        D = D,
        M = M,
        Y = Y,
        C = C,
        base_weights_name = base_weights_name,
        stabilize = stabilize,
        censor = censor,
        censor_low = censor_low,
        censor_high = censor_high,
        minimal = TRUE
      ) |>
        unlist()
      
      # adjust names
      names(boot_out) <- gsub("PSE\\.", "", names(boot_out))
      
      # output
      return(boot_out)
    }
    
    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster, 
        varlist = c("ipwpath_inner", "ipwmed_inner", "trimQ"),
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
      boot_res <- foreach::foreach(i = 1:boot_reps, .combine = rbind) %dopar% {
        boot_fnc()
      }
      boot_ATE <- boot_res[,1]
      boot_PSE <- boot_res[,-1]
    }
    else {
      boot_ATE <- rep(NA_real_, boot_reps)
      boot_PSE <- matrix(NA_real_, nrow = boot_reps, ncol = length(est$PSE))
      for (i in seq_len(boot_reps)) {
        boot_iter <- boot_fnc()
        boot_ATE[i] <- boot_iter[1]
        boot_PSE[i,] <- boot_iter[-1]
        if (i==1) {
          colnames(boot_PSE) <- names(boot_iter[-1])
        }
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
    ci_PSE <- apply(boot_PSE, MARGIN = 2, FUN = boot_ci) |>
      t()
    
    # compute two-tailed bootstrap p-values
    boot_pval <- function(x) {
      2 * min(
        mean(x < 0),
        mean(x > 0)
      )
    }
    pvalue_ATE <- boot_pval(boot_ATE)
    pvalue_PSE <- apply(boot_PSE, MARGIN = 2, FUN = boot_pval)
    
    # compile bootstrap results
    boot_out <- list(
      ci_ATE = ci_ATE,
      ci_PSE = ci_PSE,
      pvalue_ATE = pvalue_ATE,
      pvalue_PSE = pvalue_PSE,
      boot_ATE = boot_ATE,
      boot_PSE = boot_PSE
    )
  }
  
  
  # final output
  out <- est
  if (boot) {
    out <- append(out, boot_out)
  }
  return(out)
}

