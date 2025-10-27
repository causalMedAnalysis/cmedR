#' Multiply Robust (MR) Estimation for Path-specific Effects - inner function
#'
#' @description
#' Internal function used within `mrpath()`. See the `mrpath()` function
#' documentation for a description of shared function arguments. Here, we will
#' only document the one argument that is not shared by `mrpath_inner()` and
#' `mrpath()`: the `minimal` argument.
#'
#' @param minimal A logical scalar indicating whether the function should
#'   return only a minimal set of output. The `mrpath()` function uses the
#'   default of FALSE when calling `mrpath_inner()` to generate the point
#'   estimates and sets the argument to TRUE when calling `mrpath_inner()`
#'   to perform the bootstrap.
#'
#' @noRd
#---------------------- The Inner Function ----------------------#
mrpath_inner <- function(
    D,
    Y,
    M,
    C,
    data,
    d,
    dstar,
    minimal,
    censor = TRUE,
    interaction_DM = FALSE,
    interaction_DC = FALSE,
    interaction_MC = FALSE,
    censor_low = 0.01,
    censor_high = 0.99
    ){

  # Step 1: Get dimensions
  K <- length(M) # Total Number of Mediators
  PSE <- vector("list", length = K + 1) # K indirect effects through Mk; plus 1 direct effect
  models_lst_Y <- vector("list", length = K)
  models_lst_D <- vector("list", length = K)

  # Step 2: Calculate the NDE and NIE for each k
  for(k in rev(seq_len(K))) {
    # Step 2.1: Construct the estimation models based on corresponding mediators
    # D Models
    D_C_model <- as.formula(paste(D, " ~ ", paste(C, collapse= "+")))
    if(interaction_MC == FALSE) {
      D_MC_model <- as.formula(paste(D, " ~ ", paste(c(C, unlist(M[1:k])), collapse= "+")))
    }
    else {
      D_MC_model <- as.formula(
        paste(
          D,
          " ~ ",
          paste(c(C, unlist(M[1:k])), collapse= "+"),
          "+",
          paste(outer(unlist(M[1:k]), C, FUN = "paste", sep = ":"), collapse = " + ")
          )
        )
    }
    # Y Models
    if(interaction_DC == FALSE) {
      Y_DC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D), collapse= "+")))
    }
    else {
      Y_DC_model <- as.formula(
        paste(
          Y,
          " ~ ",
          paste(c(C, D), collapse= "+"),
          "+",
          paste(D, C, sep = ":", collapse = " + "))
        )
    }
    if(!any(interaction_DC, interaction_DM, interaction_MC)) {
      Y_DMC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D, unlist(M[1:k])), collapse= "+")))
    }
    else {
      Y_DMC_model <- paste(Y, " ~ ", paste(c(C, D, unlist(M[1:k])), collapse= "+"))
      if (interaction_DM) {
        Y_DMC_model <- paste(
          Y_DMC_model,
          "+",
          paste(D, unlist(M[1:k]), sep = ":", collapse = " + ")
        )
      }
      if (interaction_DC) {
        Y_DMC_model <- paste(
          Y_DMC_model,
          "+",
          paste(D, C, sep = ":", collapse = " + ")
        )
      }
      if (interaction_MC) {
        Y_DMC_model <- paste(
          Y_DMC_model,
          "+",
          paste(outer(unlist(M[1:k]), C, FUN = "paste", sep = ":"), collapse = " + ")
        )
      }
    }

    # Step 2.2: Fit the models
    est <-
      mrmed_inner(
        D = D,
        Y = Y,
        M = unlist(M[1:k]),
        C = C,
        D_C_model = D_C_model, # D ~ C
        D_MC_model = D_MC_model, # D ~ M,C
        Y_DC_model = Y_DC_model, # Y ~ D,C
        Y_DMC_model = Y_DMC_model, # Y ~ D,M,C
        M_DC_model = NULL, # M ~ D,C
        data = data,
        d = d,
        dstar = dstar,
        minimal = minimal,
        censor = censor,
        censor_low = censor_low,
        censor_high = censor_high
       )

    # Step 2.3: Calculate the PSEs based on the estimation results
    if(K == 1) {
      PSE[[K]] <- est$est2$`NDE(1,0)`
      PSE[[K + 1]] <- est$est2$`NIE(1,0)`
      ATE <- est$est2$`ATE(1,0)`
      names(PSE) <- c("D->Y","D->M1->Y")
      if(!minimal){
        models_lst_Y[[K]] <- est$models_Y
        models_lst_D[[K]] <- est$models_D
      }
    }
    else if(k == 1 & K >= 2) {
      PSE[[K - k + 1]] <- est$est2$`NDE(1,0)` - prev_NDE
      PSE[[K + 1]] <- est$est2$`NIE(1,0)`
      ATE <- est$est2$`ATE(1,0)`
      names(PSE)[K + 1] <- "D->M1~>Y"
      if(!minimal) {
        models_lst_Y[[k]] <- est$models_Y
        names(models_lst_Y)[k] <- paste0("M_", k)
        models_lst_D[[k]] <- est$models_D
        names(models_lst_D)[k] <- paste0("M_", k)
      }
      if(k + 1 == K) {
        names(PSE)[K - k + 1] <- paste0("D->M",k+1,"->Y")
      }
      else {
        names(PSE)[K - k + 1] <- paste0("D->M",k+1,"~>Y")
      }
    }
  else if(k == K & K >= 2) {
      PSE[[K - k + 1]] <- est$est2$`NDE(1,0)`
      names(PSE)[K - k + 1] <- "D->Y"
      prev_NDE <- est$est2$`NDE(1,0)`
      if(!minimal) {
        models_lst_Y[[k]] <- est$models_Y
        names(models_lst_Y)[k] <- paste0("M_", k)
        models_lst_D[[k]] <- est$models_D
        names(models_lst_D)[k] <- paste0("M_", k)
      }
    }
  else {
      PSE[[K - k + 1]] <- est$est2$`NDE(1,0)` - prev_NDE
      if(!minimal) {
        models_lst_Y[[k]] <- est$models_Y
        names(models_lst_Y)[k] <- paste0("M_", k)
        models_lst_D[[k]] <- est$models_D
        names(models_lst_D)[k] <- paste0("M_", k)
      }
      prev_NDE <- est$est2$`NDE(1,0)`
      if(k + 1 == K) {
      names(PSE)[[K - k + 1]] <- paste0("D->M",k+1,"->Y")
      }
      else {
      names(PSE)[[K - k + 1]] <- paste0("D->M",k+1,"~>Y")
      }
    }
  }

  if(minimal) {
    out <-
      list(
        ATE = ATE,
        PSE = PSE
      )
  } else {
    out <-
      list(
        ATE = ATE,
        PSE = PSE,
        models_lst_D = models_lst_D,
        models_lst_Y = models_lst_Y
      )
  }
  return(out)
}

#' Multiply Robust (MR) Estimation for Path-specific Effects
#'
#' @description
#' `mrpath()` uses a multiply robust (MR) approach to estimate path-specific effects,
#' and it computes inferential statistics using the nonparametric bootstrap.
#' If there are K causally ordered mediators, mrpath provides estimates for a direct effect
#' of the exposure on the outcome that does not operate through any of the mediators, and then
#' K path-specific effects, with each of these effects operating through one mediator, net of the
#' mediators preceding it in causal order. If only one mediator is specified, `mrpath()` computes
#' conventional natural direct and indirect effects.
#'
#' @details
#' `mrpath()` estimates path specific effects using multiply robust estimation,
#' and computes inferential statistics using the nonparametric bootstrap. It will
#' construct the multiply robust estimator based on the Type 2 estimator used in
#' the `mrmed()` function.
#'
#' To compute path-specific effects with K causally ordered mediators, `mrpath()`
#' recursively estimates a series of natural direct effects (NDEs).
#' Specifically:
#'
#' - The path-specific direct effect from treatment to outcome, bypassing all mediators,
#'   is given by \eqn{\text{PSE}_{D \to Y}(d, d^*) = \text{NDE}_{M}(d, d^*)}.
#'
#' - For each intermediate mediator \eqn{M_k}, the path-specific effect is computed
#'   as the difference between successive NDEs:
#'   \eqn{
#'   \text{PSE}_{D \to M_k \rightsquigarrow Y}(d, d^*) = \text{NDE}_{M_{k-1}}(d, d^*) - \text{NDE}_{M_k}(d, d^*).
#'   }
#'
#' - The first path-specific effect (from \eqn{ D \to M_1 \rightsquigarrow Y }) is equal
#'   to the natural indirect effect for \eqn{ M_1 }, i.e.:
#'   \eqn{
#'   \text{PSE}_{D \to M_1 \rightsquigarrow Y}(d, d^*) = \text{NIE}_{M_1}(d, d^*).
#'   }
#'
#' Specifying the `M` Argument:
#'
#' The `M` argument is a list of character vectors identifying the names of the
#' mediator variables. This argument is purposely a list of vectors rather than
#' simply a vector because it accommodates both univariate and multivariate
#' mediators treated as a set. To explain, let's start with a simple example.
#'
#' Suppose you have two single mediators, named `ever_unemp_age3539` and
#' `log_faminc_adj_age3539`, where `ever_unemp_age3539` causally precedes
#' `log_faminc_adj_age3539`. In this case, you would use the following syntax:
#' `M = list("ever_unemp_age3539", "log_faminc_adj_age3539")`.
#'
#' Now, let's say you have a third mediator, named `m3`. You believe that
#' `ever_unemp_age3539` causally precedes both `log_faminc_adj_age3539` and
#' `m3`, but you are unwilling to make an assumption about the relative causal
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
#' `M = list("ever_unemp_age3539", c("log_faminc_adj_age3539", "m3"), c("level2","level3","level4"))`.
#'
#' @param data A data frame.
#' @param D A character scalar identifying the name of the exposure variable in
#'   `data`. `D` is a character string, but the exposure variable it identifies
#'   must be numeric and binary, with two distinct values.
#' @param M A character list (of one or more elements) identifying the names
#'   of the mediator variables in `data`. If you are estimating univariate
#'   natural effects (with a single mediator), `M` should be a character scalar
#'   (i.e., a vector with only one element). If you are estimating path-specific
#'   effects, `M` should be a list identifying all mediators.
#' @param Y A character scalar identifying the name of the outcome variable in
#'   `data`. `Y` is a character string, but the outcome variable it identifies
#'   must be numeric.
#' @param C A character vector identifying the names of the baseline confounders
#'   in `data` that you wish to include in both the mediator and outcome models.
#'   If there are no such covariates you wish to include, leave `C` as its
#'   default null argument.
#' @param interaction_DM A logical scalar indicating whether the outcome model
#'   should include exposure-mediator interactions (interactions of the exposure
#'   with each mediator if there is more than one mediator in `M`).
#' @param interaction_DC A logical scalar indicating whether the outcome
#'   model should include interactions of the exposure with each covariate in `C`.
#' @param interaction_MC A logical scalar indicating whether the outcome model
#'   should include interactions of each mediator in `M` with each covariate in
#'   `C`.
#' @param d The numeric value of the treatment variable that the user defines as
#'   the treatment status. If not equal to 1, the function will recode it as 1.
#' @param dstar The numeric value of the treatment variable that the user defines
#'   as the control status. If not equal to 0, the function will recode it as 0.
#' @param censor A logical scalar indicating whether the IPW weights constructed by
#'   estimation procedure should be censored. By default, this value is `TRUE`.
#' @param censor_low,censor_high A pair of arguments, each a numeric scalar
#'   denoting a probability in \[0,1\]. If `censor` is TRUE, then IPW weights below
#'   the `censor_low` quantile will be bottom-coded, and weights above the
#'   `censor_high` quantile will be top-coded. For example, if the default values
#'   `censor_low = 0.01` and `censor_high = 0.99` are used, then IPW weights will
#'   be censored at their 1st and 99th percentiles.
#' @param boot A logical scalar indicating whether the function should perform
#'   the nonparametric bootstrap and return two-sided confidence intervals and
#'   p-values.
#' @param boot_reps An integer scalar specifying the number of bootstrap replications
#'   to perform.
#' @param boot_conf_level A numeric scalar specifying the confidence level for the
#'   bootstrap interval.
#' @param boot_seed An integer scalar specifying the random-number seed used in
#'   bootstrap resampling.
#' @param boot_parallel A logical scalar indicating whether the bootstrap should
#'   be performed using a parallelized loop to reduce runtime. Parallel computation,
#'   as implemented in this function, requires that the following R packages are installed:
#'   `doParallel`, `doRNG`, and `foreach`. However, you do not need to explicitly
#'   load these packages using `library()`. Note that the results of the parallelized
#'   bootstrap may differ slightly from those of the non-parallelized bootstrap, even if
#'   the same seed is specified, due to differences in how seeds are processed.
#' @param boot_cores An integer scalar specifying the number of CPU cores to use
#'   for the parallelized bootstrap. This argument only affects computation if both
#'   `boot` and `boot_parallel` are TRUE. By default, `boot_cores` is set to the greater
#'   of two values: (a) one, and (b) the number of available CPU cores minus two.
#'   If `boot_cores` equals one, the bootstrap loop will not be parallelized,
#'   regardless of the value of `boot_parallel`.
#'
#' @returns Based on the user's specification, `mrpath()` returns the
#' following elements:
#'
#' \item{ATE:}{A numeric scalar with the estimated average total effect
#'   for the exposure contrast `d - dstar`: ATE(`d`,`dstar`).}
#' \item{PSE:}{A numeric list of length `length(M)+1` with the estimated
#'   path-specific effects for the exposure contrast `d - dstar`. The vector is
#'   named with the path each effect describes.}
#' \item{model_lst_D:}{A list of fitted treatment models, where each element corresponds
#' to a mediator variable. If multiple mediators are included, the list stores
#' separate models for each.}
#' \item{model_lst_Y:}{A list of models regressing the outcome on the treatment,
#' controls, and an increasing sequence of mediators. The first model (`M_1`)
#' includes only the first mediator,the second (`M_2`) includes the first two,
#' the third (`M_3`) includes the first three, and so on.}
#'
#' If you request the bootstrap (by setting the `boot` argument to TRUE), then
#' the function returns all of the elements listed above, as well as the
#' following additional elements:
#'  \describe{
#'   \item{ATE}{List of length 2 with bootstrap results for the Average Treatment Effect.}
#'   \item{PSE: D -> Mk ~>Y or D -> Mk -> Y for the last mediator in the causal chain}{
#'   List of length 2 with bootstrap results for the path-specific effect from D to Y (through Mk).}
#'   \item{PSE:D->Y}{List of length 2 with bootstrap results for the direct effect of D on Y.}
#' }
#'
#' Each sub-list contains:
#' \describe{
#'   \item{stats}{A \code{tibble} with 1 row and 4 columns:
#'     \describe{
#'       \item{\code{ci}}{Character. Confidence interval as a formatted string (e.g., \code{"[-0.25, -0.031]"})}
#'       \item{\code{pvalue}}{Numeric. P-value for the estimate.}
#'       \item{\code{sd}}{Numeric. Standard deviation (bootstrap standard error).}
#'       \item{\code{method_type}}{Character. A label for the estimand (e.g., \code{"ATE"}).}
#'     }
#'   }
#'   \item{org_val}{A \code{tibble} with 3 columns:
#'     \describe{
#'       \item{\code{value}}{Numeric. The bootstrapped estimate in each iteration.}
#'       \item{\code{method_type}}{Character. Estimand repeated for each row.}
#'       \item{\code{boot_id}}{Integer. Bootstrap replicate ID.}
#'     }
#'   }
#' }
#'
#' @export
#'
#' @examples
#' # ----------------------------- #
#' #     Data and shared setup     #
#' # ----------------------------- #
#' data(nlsy)
#'
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
#'
#' key_variables <- c(
#'   "cesd_age40",
#'   "ever_unemp_age3539",
#'   "log_faminc_adj_age3539",
#'   "att22",
#'   covariates
#' )
#'
#' # For convenience in examples, use complete cases
#' nlsy_ex <- nlsy[complete.cases(nlsy[, key_variables]), ]
#' nlsy_ex$std_cesd_age40 <-
#'   (nlsy_ex$cesd_age40 - mean(nlsy_ex$cesd_age40)) /
#'   sd(nlsy_ex$cesd_age40)
#'
#' # ----------------------------------------- #
#' # Example 1: Two mediators, additive models #
#' # ----------------------------------------- #
#' mrpath(
#'   data = nlsy_ex,
#'   D = "att22",
#'   M = list("ever_unemp_age3539", "log_faminc_adj_age3539"),
#'   Y = "std_cesd_age40",
#'   C = covariates,
#'   d = 1,
#'   dstar = 0,
#'   boot = FALSE
#' )
#'
#' # -------------------------------------------------- #
#' # Example 2: Two mediators, models with interactions #
#' # -------------------------------------------------- #
#' mrpath(
#'   data = nlsy_ex,
#'   D = "att22",
#'   M = list("ever_unemp_age3539", "log_faminc_adj_age3539"),
#'   Y = "std_cesd_age40",
#'   C = covariates,
#'   d = 1,
#'   dstar = 0,
#'   interaction_DM = TRUE,
#'   interaction_DC = TRUE,
#'   interaction_MC = TRUE
#' )
#'
#' # --------------------------- #
#' # Example 3: Single mediator  #
#' # Returns ATE, NDE, and NIE   #
#' # --------------------------- #
#' mrpath(
#'   data = nlsy_ex,
#'   D = "att22",
#'   M = "ever_unemp_age3539",
#'   Y = "std_cesd_age40",
#'   C = covariates,
#'   d = 1,
#'   dstar = 0
#' )
#'
#' # ----------------------------------------- #
#' # Example 4: Bootstrap with parallelization #
#' # ----------------------------------------- #
#' mrpath(
#'   data = nlsy_ex,
#'   D = "att22",
#'   M = list("ever_unemp_age3539", "log_faminc_adj_age3539"),
#'   Y = "std_cesd_age40",
#'   C = covariates,
#'   d = 1,
#'   dstar = 0,
#'   boot = TRUE,
#'   boot_reps = 200,
#'   boot_seed = 1234,
#'   boot_parallel = FALSE
#' )
#'
#' # ------------------------------------------- #
#' # Example 5: Three mediators, additive models #
#' # ------------------------------------------- #
#'
#' key_variables5 <- c(
#'   "cesd_age40","cesd_1992","ever_unemp_age3539",
#'   "log_faminc_adj_age3539","att22", covariates
#' )
#' nlsy_ex5 <- nlsy[complete.cases(nlsy[, key_variables5]), ]
#' nlsy_ex5$std_cesd_age40 <-
#'   (nlsy_ex5$cesd_age40 - mean(nlsy_ex5$cesd_age40)) /
#'   sd(nlsy_ex5$cesd_age40)
#'
#' mrpath(
#'   data = nlsy_ex5,
#'   D = "att22",
#'   M = list("cesd_1992","ever_unemp_age3539","log_faminc_adj_age3539"),
#'   Y = "std_cesd_age40",
#'   C = covariates,
#'   d = 1,
#'   dstar = 0
#' )
#'
mrpath <- function(
    D,
    Y,
    M,
    C,
    data,
    d,
    dstar,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99,
    interaction_DM = FALSE,
    interaction_DC = FALSE,
    interaction_MC = FALSE,
    boot = FALSE,
    boot_reps = 200,
    boot_conf_level = 0.95,
    boot_seed = NULL,
    boot_parallel = FALSE,
    boot_cores = max(c(parallel::detectCores()-2,1))
){
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

  # compute point estimates
  est <-
    mrpath_inner(
      D = D,
      Y = Y,
      M = M,
      C = C,
      data = data,
      d = d,
      dstar = dstar,
      minimal = FALSE,
      censor = censor,
      interaction_DM = interaction_DM,
      interaction_DC = interaction_DC,
      interaction_MC = interaction_MC,
      censor_low = censor_low,
      censor_high = censor_high
    )

  # predeclare boot_res
  boot_res <- NULL
  boot_res_lst <- NULL

  # bootstrap, if requested
  if (boot) {
    boot_fnc <- function() {

      boot_data <- data_outer %>% dplyr::sample_frac(replace = TRUE)

      mrpath_inner(
        D = D,
        Y = Y,
        M = M,
        C = C,
        data = boot_data,
        d = d,
        dstar = dstar,
        minimal = TRUE,
        censor = TRUE,
        interaction_DM = interaction_DM,
        interaction_DC = interaction_DC,
        interaction_MC = interaction_MC,
        censor_low = censor_low,
        censor_high = censor_high
      )
    }

    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type = "PSOCK")
      doParallel::registerDoParallel(cl = x_cluster)
      parallel::clusterExport(
        cl = x_cluster,
        varlist = c("mrpath_inner", "mrmed_inner", "trimQ", "censor", "boot_fnc"),
        envir = environment()
      )
      `%dopar%` <- foreach::`%dopar%`
    }

    if (!is.null(boot_seed)) {
      set.seed(boot_seed)
      if (boot_parallel_rev) {
        doRNG::registerDoRNG(boot_seed)
      }
    }

    if (boot_parallel_rev) {
      boot_res <- foreach::foreach(
        i = 1:boot_reps,
        .combine = dplyr::bind_rows,
        .packages = c("dplyr", "rlang", "tidyr", "purrr", "Hmisc", "tibble")
        ) %dopar% {
          out <- boot_fnc()
          out_filtered <- out[!sapply(out, is.null)]
          purrr::imap_dfr(out_filtered, function(.x, .y) {
          if (is.list(.x)) {
            purrr::imap_dfr(.x, function(.xx, .yy) {
              tibble::tibble(value = .xx, method_type = paste(.y, .yy, sep = ":"), boot_id = i)
            })
          }
            else {
              tibble::tibble(value = .x, method_type = .y, boot_id = i)
              }
            })
          }
    }
    else {
      boot_res <- dplyr::bind_rows(
        lapply(seq_len(boot_reps), function(i) {
          out <- boot_fnc()
          out_filtered <- out[!sapply(out, is.null)]
          purrr::imap_dfr(out_filtered, function(.x, .y) {
            if (is.list(.x)) {
              purrr::imap_dfr(.x, function(.xx, .yy) {
                tibble::tibble(value = .xx, method_type = paste(.y, .yy, sep = ":"), boot_id = i)
              })
            } else {
              tibble::tibble(value = .x, method_type = .y, boot_id = i)
            }
          })
        })
      )
    }

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

  boot_pval <- function(x) {
    2 * min(
      mean(x < 0),
      mean(x > 0)
    )
  }

  boot_res_lst <-
    lapply(
      split(boot_res, boot_res$method_type),
      function(rst){
        stat_df <-
          tibble::tibble(
            ci = paste(
              "[",
              round(boot_ci(rst$value)[1],3),",",round(boot_ci(rst$value)[2],3),
              "]"),
            pvalue = boot_pval(rst$value),
            sd = sd(rst$value, na.rm = TRUE),
            method_type = unique(rst$method_type)
          )
        return(
          list(
            stats = stat_df,
            org_val = rst
          )
        )
      }
    )
  }

  if (boot) {
    out <- list(point_rst = est, boot_rst = boot_res_lst)
  } else {
    out <- est
  }
  return(out)
}
