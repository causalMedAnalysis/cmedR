#' Regression-Imputation estimator for path-specific effects
#'
#' @description
#' `pathimp()` is a wrapper for two functions from the `paths` package. It
#' implements the pure imputation estimator and the imputation-based weighting
#' estimator (when a propensity score model is provided) as detailed in Zhou
#' and Yamamoto (2020). You may install the `paths` package from CRAN or GitHub:
#' `devtools::install_github("xiangzhou09/paths")`
#'
#' @details
#' `pathimp()` estimates path-specific effects using pure regression imputation and (optionally)
#' an imputation-based weighting estimator, and it computes inferential statistics using the
#' nonparametric bootstrap.
#'
#' With K causally ordered mediators, the implementation proceeds as follows:
#' (i) it fits a model for the mean of the outcome conditional on the exposure and baseline confounders;
#' (ii) it imputes conventional potential outcomes under using model from (i); (iii) for each mediator
#' k = 1, 2, ..., K, it then fits (iiia) a model for the mean of the outcome conditional on the exposure,
#' baseline confounders, and the mediators Mk = \{M1, ..., Mk\}; (iv) it uses the models from (iii) to
#' impute cross-world potential outcomes; and (v) and finally, it uses the imputed outcomes from all
#' the previous steps to calculate estimates for the path-specific effects.
#'
#' `pathimp()` provides estimates for the total effect and K+1 path-specific effects: the direct effect
#' of the exposure on the outcome that does not operate through any of the mediators, and separate
#' path-specific effects operating through each of the K mediators, net of the mediators that precede
#' them in causal order.
#'
#' If only a single mediator is specified, `pathimp()` reverts to estimates of conventional natural
#' direct and indirect effects through a univariate mediator.
#'
#' @param data A data frame.
#'
#' @param D A character string indicating the name of the treatment variable in
#' `data`. The treatment should be a binary variable taking either 0 or 1.
#'
#' @param Y A character string indicating the name of the outcome variable.
#'
#' @param M A list of \eqn{K} character vectors indicating the names of \eqn{K}
#' causally ordered mediators \eqn{M_1,\ldots, M_K}.
#'
#' @param Y_models A list of \eqn{K+1} fitted model objects describing how the
#' outcome depends on treatment, pretreatment confounders, and varying sets of
#' mediators, where \eqn{K} is the number of mediators.
#' \itemize{
#'   \item the first element is a baseline model of the outcome conditional on
#'   treatment and pretreatment confounders.
#'   \item the \eqn{k}th element is an outcome model conditional on treatment,
#'   pretreatment confounders, and mediators \eqn{M_1,\ldots, M_{k-1}}.
#'   \item the last element is an outcome model conditional on treatment,
#'   pretreatment confounders, and all of the mediators, i.e.,
#'   \eqn{M_1,\ldots, M_K}.
#'   }
#'  The fitted model objects can be of type \code{\link{lm}}, \code{\link{glm}},
#'  \code{\link[gbm]{gbm}}, \code{\link[BART]{wbart}}, or \code{\link[BART]{pbart}}.
#'
#' @param D_model An optional propensity score model for treatment. It can be
#' of type \code{\link{glm}},\code{\link[gbm]{gbm}}, \code{\link[twang]{ps}}, or
#' \code{\link[BART]{pbart}}. When it is provided, the imputation-based weighting
#' estimator is also used to compute path-specific causal effects. Defaults to
#' \code{NULL}. Must be explicitly specified when \code{out_ipw = TRUE} to compute the
#' imputation-based weighting estimator.
#' @param boot A logical scalar indicating whether the function will perform the
#'   nonparametric bootstrap and return two-sided confidence intervals and
#'   p-values. Defaults to FALSE, in which case only point estimates are returned.
#' @param boot_reps An integer scalar for the number of bootstrap replications to
#'   perform. Only consulted when \code{boot = TRUE}.
#' @param boot_conf_level A numeric scalar for the confidence level of the
#'   bootstrap interval.
#' @param boot_seed An integer scalar specifying the random-number seed used in
#'   bootstrap resampling. Defaults to \code{NULL}, indicating that no seed is set.
#' @param boot_parallel A logical scalar indicating whether the bootstrap will
#' be performed with a parallelized loop, with the goal of reducing runtime.
#' Defaults to FALSE. Parallelization uses the `"snow"` backend of
#' \code{\link[boot]{boot}} (via the `paths` package), which works across all
#' operating systems.
#' @param round_decimal The number of decimal digits to which results are
#' rounded and displayed.
#' @param boot_cores An integer scalar specifying the number of CPU cores on
#' which the parallelized bootstrap will run, passed as `ncpus` to the `paths`
#' bootstrap. This argument only has an effect if you requested a parallelized
#' bootstrap (i.e., only if `boot` is TRUE and `boot_parallel` is TRUE). By
#' default, `boot_cores` is equal to the greater of two values: (a) one and (b)
#' the number of available CPU cores minus two. If `boot_cores` equals one, then
#' the bootstrap loop will not be parallelized (regardless of the setting of
#' `boot_parallel`).
#' @param out_ipw A logical value indicating whether to report the
#' imputation-based weighting estimator. If set to \code{TRUE}, the user must
#' specify the propensity score model to calculate the Imputation-based Weighting
#' Estimator . If \code{FALSE}, only the Pure Imputation Estimator will be
#' returned. Defaults to \code{FALSE}.
#'
#' @returns When \code{out_ipw = TRUE}, `pathimp()` returns a data frame with the
#' following information: \describe{
#' \item{Pure Imputation Estimator}{Estimates the direct (ATE) and path-specific
#' effects (PSE)  through mediators \eqn{M_1, \ldots, M_K} using the pure
#' imputation estimator, along with corresponding bootstrap confidence intervals.}
#' \item{Imputation-based Weighting Estimator}{Estimates of direct and
#' path-specific effects via \eqn{M_1, \ldots, M_K} based on the imputation-based
#' weighting estimator,along with corresponding bootstrap confidence intervals.}
#' } When \code{out_ipw = FALSE}, only the pure imputation estimator will be
#' returned.
#'
#' @importFrom paths paths
#' @importFrom methods is
#' @importFrom utils capture.output
#'
#' @export
#'
#' @examples
#' # Example 1: Pure imputation with two mediators
#' ## Prepare data:
#' data(nlsy)
#' covariates <- c(
#' "female",
#' "black",
#' "hispan",
#' "paredu",
#' "parprof",
#' "parinc_prank",
#' "famsize",
#' "afqt3"
#' )
#'
#' key_variables <- c(
#' "cesd_age40",
#' "ever_unemp_age3539",
#' "log_faminc_adj_age3539",
#' "att22",
#' covariates
#' )
#'
#' df <-
#' nlsy[complete.cases(nlsy[,key_variables]),] |>
#' dplyr::mutate(
#' std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
#' )
#'
#' glm_m0 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22,
#'   data = df
#' )
#' glm_m1 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 + ever_unemp_age3539,
#'   data = df
#' )
#' glm_m2 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 + ever_unemp_age3539 +
#'     log_faminc_adj_age3539,
#'   data = df
#' )
#' glm_ymodels <- list(glm_m0, glm_m1, glm_m2)
#'
#' # Fit the paths model:
#' fit_ex1 <- pathimp(
#'   D = "att22",
#'   Y = "std_cesd_age40",
#'   M = list("ever_unemp_age3539","log_faminc_adj_age3539"),
#'   Y_models = glm_ymodels,
#'   data = df,
#'   boot = TRUE,
#'   boot_reps = 250,
#'   boot_conf_level = 0.95,
#'   boot_seed = 2138,
#'   out_ipw = FALSE
#' )
#' print(fit_ex1)
#'
#' # Example 2: Adding imputation-based weighting
#' glm_ps <- glm(
#'   att22 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3,
#'   family = binomial("logit"),
#'   data = df
#' )
#' # Fit the paths model:
#' fit_ex2 <- pathimp(
#'   D = "att22",
#'   Y = "std_cesd_age40",
#'   M = list("ever_unemp_age3539","log_faminc_adj_age3539"),
#'   Y_models = glm_ymodels,
#'   D_model = glm_ps,
#'   data = df,
#'   boot = TRUE,
#'   boot_reps = 250,
#'   boot_conf_level = 0.95,
#'   boot_seed = 2138,
#'   out_ipw = TRUE
#' )
#' print(fit_ex2)
#'
#' # Example 3: Pure imputation with three mediators
#'
#' # Prepare data
#' key_variables3 <- c(
#'   "cesd_age40","cesd_1992","ever_unemp_age3539","log_faminc_adj_age3539",
#'   "att22", covariates
#' )
#' df3 <- nlsy[complete.cases(nlsy[,key_variables3]),] |>
#'   dplyr::mutate(
#'     std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
#'   )
#'
#' glm_m0 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22,
#'   data = df3
#' )
#' glm_m1 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 + cesd_1992,
#'   data = df3
#' )
#' glm_m2 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 + cesd_1992 +
#'     ever_unemp_age3539,
#'   data = df3
#' )
#' glm_m3 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 + cesd_1992 +
#'     ever_unemp_age3539 + log_faminc_adj_age3539,
#'   data = df3
#' )
#' glm_ymodels3 <- list(glm_m0, glm_m1, glm_m2, glm_m3)
#'
#' # Fit the paths model:
#' \dontrun{
#' fit_ex3 <- pathimp(
#'   D = "att22",
#'   Y = "std_cesd_age40",
#'   M = list("cesd_1992","ever_unemp_age3539","log_faminc_adj_age3539"),
#'   Y_models = glm_ymodels3,
#'   data = df3,
#'   boot = TRUE,
#'   boot_reps = 250,
#'   boot_conf_level = 0.95,
#'   boot_seed = 2138,
#'   boot_parallel = TRUE, # enable parallel bootstrap if boot = TRUE
#'   out_ipw = FALSE
#' )
#' print(fit_ex3)
#'}
#'
#' # Example 4: Pure imputation with two mediators including treatment–mediator
#' #            interactions
#'
#' glm_m0 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22,
#'   data = df
#' )
#'
#' # Outcome models including treatment-mediator interactions
#' glm_m1 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 * ever_unemp_age3539,
#'   data = df
#' )
#' glm_m2 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 *
#'     (ever_unemp_age3539 + log_faminc_adj_age3539) ,
#'   data = df
#' )
#' glm_ymodels_intM <- list(glm_m0, glm_m1, glm_m2)
#'
#' # Fit the paths model:
#' fit_ex4 <- pathimp(
#'   D = "att22",
#'   Y = "std_cesd_age40",
#'   M = list("ever_unemp_age3539","log_faminc_adj_age3539"),
#'   Y_models = glm_ymodels_intM,
#'   data = df,
#'   boot = TRUE,
#'   boot_reps = 250,
#'   boot_conf_level = 0.95,
#'   boot_seed = 2138,
#'   out_ipw = FALSE
#' )
#' print(fit_ex4)
#'
#'
#' # Example 5: Pure imputation with two mediators including treatment–mediator
#' #            and treatment–confounder interactions
#'
#' glm_m0 <- glm(
#'   std_cesd_age40 ~ (female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3) * att22,
#'   data = df
#' )
#'
#' # Outcome models including treatment-mediator and treatment-confounder
#' # interactions
#' glm_m1 <- glm(
#'   std_cesd_age40 ~ (female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + ever_unemp_age3539) * att22,
#'   data = df
#' )
#' glm_m2 <- glm(
#'   std_cesd_age40 ~ (female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + ever_unemp_age3539 +
#'     log_faminc_adj_age3539) * att22,
#'   data = df
#' )
#' glm_ymodels_intMC <- list(glm_m0, glm_m1, glm_m2)
#'
#' # Fit the paths model:
#' fit_ex5 <- pathimp(
#'   D = "att22",
#'   Y = "std_cesd_age40",
#'   M = list("ever_unemp_age3539","log_faminc_adj_age3539"),
#'   Y_models = glm_ymodels_intMC,
#'   data = df,
#'   boot = TRUE,
#'   boot_reps = 250,
#'   boot_conf_level = 0.95,
#'   boot_seed = 2138,
#'   out_ipw = FALSE
#' )
#' print(fit_ex5)
#'
#' # Example 6: Pure imputation with 1 mediator
#'
#' # Prepare data
#' key_variables1 <- c("cesd_age40","ever_unemp_age3539","att22", covariates)
#'
#' df1 <- nlsy[complete.cases(nlsy[,key_variables1]),] |>
#'   dplyr::mutate(
#'     std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
#'   )
#'
#' glm_m0 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22,
#'   data = df1
#' )
#' glm_m1 <- glm(
#'   std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 + ever_unemp_age3539,
#'   data = df1
#' )
#' glm_ymodels1 <- list(glm_m0, glm_m1)
#'
#' # Fit the paths model:
#' fit_ex6 <- pathimp(
#'   D = "att22",
#'   Y = "std_cesd_age40",
#'   M = list("ever_unemp_age3539"),
#'   Y_models = glm_ymodels1,
#'   data = df1,
#'   boot = TRUE,
#'   boot_reps = 250,
#'   boot_conf_level = 0.95,
#'   boot_seed = 2138,
#'   out_ipw = FALSE
#' )
#' print(fit_ex6)
#'
#' # Note: With a single mediator, the effect labeled as the path-specific
#' # effect (D -> M -> Y) is the natural indirect effect (NIE), and the
#' # effect of the direct path (D -> Y) is the natural direct effect (NDE).
#'

pathimp <- function(
    D,
    Y,
    M,
    Y_models,
    D_model = NULL,
    data,
    boot = FALSE,
    boot_reps = 200,
		boot_conf_level = 0.95,
		boot_seed = NULL,
		boot_parallel = FALSE,
		round_decimal = 3,
		boot_cores = max(c(parallel::detectCores() - 2L, 1L), na.rm = TRUE),
    out_ipw = FALSE){
  # For Y_models, check model type and model arguments:
  for(i in seq_len(length(Y_models))){
    # Check model type:
    model_type <- is(Y_models[[i]])[1]
    model <- Y_models[[i]]
    if(!model_type %in% c("lm","glm","gbm","wbart","pbart")){
      stop(paste(
        "The model type must be lm, glm, gbm, pbart or wbart"))
    }
  # Grab the regressors:
    if(model_type %in% c("pbart","wbart")){
      regressors <- colnames(model$varprob)
    }else if(model_type %in%c("lm","glm")){
      regressors <- as.character(attr(model$terms,"variables"))[-c(1:2)]
    }else{
      regressors <- as.character(attr(model$Terms,"variables"))[-c(1:2)]
    }

    # Check the arguments of the M models:
    if(i == 1){
      if(sum(unlist(M) %in% regressors) > 0){
        stop(paste(
        "The first model should regress the outcome variable only on controls."
        ))}
      regressors_DC <- regressors
    }else{
      Mk <- unlist(M[c(1:(i-1))])
      if(!setequal(setdiff(regressors, regressors_DC),Mk)){
        stop(
          "Please double-check your model specification; the order of mediators
          in the list should match the order of the specified models."
        )
      }
    }
  }

  # Check the model type of the D_model:
  if(out_ipw == TRUE){
    if(is.null(D_model)){
      stop("Please specify your propensity score model for treatment.")
    }
    model_type_ps <- is(D_model)[1]
    if(!model_type_ps %in% c("glm","gbm","ps","pbart")){
      stop(paste(
        "The model type must be glm, gbm, ps or pbart"))
    }
  }

  # Translate the logical boot_parallel toggle into the backend paths() forwards
  # to boot::boot. "snow" is the only cross-platform backend ("multicore" fails
  # on Windows). Fall back to sequential when only one core is available.
  boot_parallel_use <- if(isTRUE(boot_parallel) && boot_cores > 1L){
    "snow"
  }else{
    "no"
  }

  # Set seed. Snow workers run in separate memory, so a plain set.seed() in the
  # master does not reach them. Switching to the L'Ecuyer-CMRG generator lets
  # boot::boot (inside paths()) distribute reproducible RNG streams to the snow
  # workers via clusterSetRNGStream(), making same-seed runs reproducible.
  if(!is.null(boot_seed)){
    if(boot_parallel_use == "snow"){
      old_rng_kind <- RNGkind("L'Ecuyer-CMRG")[[1L]]
      on.exit(RNGkind(old_rng_kind), add = TRUE)
    }
    set.seed(boot_seed)
  }

  # Under snow, paths() refits each model in a fresh worker process by
  # re-evaluating the model's stored call. If a model was fit by passing the
  # formula as a global object (e.g. lm(yform, data = ...)), the worker cannot
  # find that name and errors with "object '<name>' not found". Inline the
  # resolved formula into each model's call so the refit is self-contained.
  if(boot_parallel_use == "snow"){
    resolve_model_formula <- function(model){
      if(is.null(model)) return(model)
      model_call <- model[["call"]]
      if(is.null(model_call) || is.null(model_call[["formula"]])) return(model)
      resolved <- tryCatch(stats::formula(model), error = function(e) NULL)
      if(!is.null(resolved)) model[["call"]][["formula"]] <- resolved
      model
    }
    Y_models <- lapply(Y_models, resolve_model_formula)
    if(!is.null(D_model)) D_model <- resolve_model_formula(D_model)
  }

  # Fit the paths model. paths() always bootstraps, so for point estimates only
  # (boot = FALSE) we force the minimum nboot and discard the inference. The
  # point estimates (boot_out$t0) do not depend on nboot, so they are unchanged.
  fit_paths_model <- function(nboot, parallel = boot_parallel_use, ncpus = boot_cores){
    result <- NULL
    capture.output(
      result <- paths(
        a = D,
        y = Y,
        m = M,
        models = Y_models,
        ps_model = D_model,
        data = data,
        nboot = nboot,
        conf_level = boot_conf_level,
        ncpus = ncpus,
        parallel = parallel
      )
    )
    result
  }

  if(boot){
    paths_model <- fit_paths_model(boot_reps)
  }else{
    # try nboot = 1, fall back to 2 if the install/backend rejects it.
    paths_model <- tryCatch(
      fit_paths_model(1L, parallel = "no", ncpus = 1L),
      error = function(e) fit_paths_model(2L, parallel = "no", ncpus = 1L)
    )
  }

  # Clean the model output:
  result <- list(paths_model$pure, paths_model$hybrid)
  processed_result <-
    lapply(
      result,
      function(rst_df){
        rst_df <-
          rst_df %>%
          dplyr::filter(
            .data$decomposition == "Type I"
          ) %>%
          dplyr::mutate(
            estimand = dplyr::case_when(
              .data$estimand == "direct" ~ "PSE(D -> Y)",
              .data$estimand == "total" ~ "ATE(1,0)",
              stringr::str_detect(.data$estimand, "via M\\d+") ~
                dplyr::if_else(
                  stringr::str_detect(.data$estimand, paste0("M", length(M))),
                  stringr::str_replace(.data$estimand, "via (M\\d+)", "PSE(D -> \\1 -> Y)"),
                  stringr::str_replace(.data$estimand, "via (M\\d+)", "PSE(D -> \\1 ~> Y)")
                ),
              TRUE ~ .data$estimand
            ),
            est = round(.data$estimate, round_decimal),
            intv = paste0("[", round(.data$lower, round_decimal), ", ", round(.data$upper, round_decimal), "]"),
            out = if(boot) paste(.data$est, .data$intv) else as.character(.data$est)
          ) %>%
          dplyr::mutate(
            order_key = dplyr::case_when(
              stringr::str_detect(.data$estimand, "ATE") ~ 1,
              estimand == "PSE(D -> Y)" ~ 2,
              TRUE ~ 3
            )
          ) %>%
          dplyr::arrange(.data$order_key) %>%
          dplyr::select(-.data$order_key) %>%
          dplyr::select(
            .data$estimator,
            .data$estimand,
            .data$out
          ) %>%
          dplyr::mutate(
            estimator = dplyr::case_when(
              estimator == "pure" ~ "Pure Imputation Estimator",
              estimator == "hybrid" ~ "Imputation-based Weighting Estimator"
            )
          )
      }
    )
  if(out_ipw == TRUE){
    summary_df <- rbind(processed_result[[1]], processed_result[[2]])
  }else{
    summary_df <- processed_result[[1]]
  }

  if(boot){
    org_obj <- paths_model
  }else{
    # Slimmer, inference-free object (cf. impmed when boot = FALSE): drop the
    # se/CI/p columns, the bootstrap draws, and the bootstrap settings.
    org_obj <- paths_model
    for(component in c("pure", "hybrid")){
      comp_df <- org_obj[[component]]
      if(is.data.frame(comp_df) && nrow(comp_df) > 0){
        keep_cols <- setdiff(names(comp_df), c("se", "lower", "upper", "p"))
        org_obj[[component]] <- comp_df[, keep_cols, drop = FALSE]
      }
    }
    org_obj$boot_out <- NULL
    org_obj$nboot <- NA_integer_
    org_obj$conf_level <- NA_real_
  }

  return(list(summary_df = summary_df, org_obj = org_obj))
}
