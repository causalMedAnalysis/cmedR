#' Regression-Imputation estimator for path-specific effects
#'
#' @description
#' `pathimp()` is a wrapper for two functions from the `paths` package. It
#' implements the pure imputation estimator and the imputation-based weighting
#' estimator (when a propensity score model is provided) as detailed in Zhou
#' and Yamamoto (2020). You may install the `paths` package from CRAN or Github:
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
#' @param boot_reps An integer scalar for the number of bootstrap replications to perform.
#' @param boot_conf_level A numeric scalar for the confidence level of the
#'   bootstrap interval.
#' @param boot_seed An integer scalar specifying the random-number seed used in
#'   bootstrap resampling. Defaults to \code{NULL}, indicating that no seed is set.
#' @param boot_parallel Type of parallel operation to be used. Follows the same
#' specification as the \code{\link[boot]{boot}} package. Defaults to \code{"no"}.
#' @param round_decimal The number of decimal digits to which results are
#' rounded and displayed.
#' @param boot_cores An integer scalar specifying the number of CPU cores on
#' which the parallelized bootstrap will run. This argument only has an effect
#' if you requested a parallelized bootstrap (i.e., only if `boot_parallel` is
#' not `no`). By default, `boot_cores` is equal to the greater of two values:
#' (a) one and (b) the number of available CPU cores minus two. If `boot_cores`
#' equals one, then the bootstrap loop will not be parallelized (regardless of
#' the setting of `boot_parallel`).
#' @param out_ipw A logical value indicating whether to report the
#' imputation-based weighting estimator. If set to \code{TRUE}, the user must
#' specify the propensity score model to calculate the Imputation-based Weighting
#' Estimator . If \code{FALSE}, only the Pure Imputation Estimator will be
#' returned.
#'
#' @returns When \code{out_ipw = TRUE}, `pathimp()` returns a dataframe with the
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
#'
#' @export
#'
#' @examples
#' # Example 1: With one exposure-induced confounder
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
#' key_variables <- c(
#' "cesd_age40",
#' "ever_unemp_age3539",
#' "log_faminc_adj_age3539",
#' "att22",
#' covariates
#' )
#' df <-
#' nlsy[complete.cases(nlsy[,key_variables]),] |>
#' dplyr::mutate(
#' std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
#' )
#' ## Set up the Y_models:
#'   glm_m0 <- glm(
#'     std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'      parinc_prank + famsize + afqt3 + att22,
#'     data = df
#'     )
#'
#'   glm_m1 <- glm(
#'     std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 +  ever_unemp_age3539,
#'     data = df
#'     )
#'
#'   glm_m2 <- glm(
#'     std_cesd_age40 ~ female + black + hispan + paredu + parprof +
#'     parinc_prank + famsize + afqt3 + att22 + ever_unemp_age3539 +
#'     log_faminc_adj_age3539,
#'     data = df
#'     )
#'
#' glm_ymodels <- list(glm_m0, glm_m1, glm_m2)
#'
#'
#' ## Set up the D_model:
#' glm_ps <- glm(
#' att22 ~ female + black + hispan + paredu + parprof + parinc_prank +
#' famsize + afqt3,
#' family = binomial("logit"),
#' data = df
#' )
#'
#' # Example 1: Including Imputation-based Weighting Estimator:
#'
#' ## Fit the paths model:
#' glm_paths <-
#' pathimp(
#' D = "att22",
#' Y = "std_cesd_age40",
#' M = list("ever_unemp_age3539","log_faminc_adj_age3539"),
#' Y_models = glm_ymodels,
#' D_model = glm_ps,
#' data = df,
#' boot_reps = 250,
#' boot_conf_level = 0.95,
#' boot_seed = 2138,
#' boot_parallel = "multicore", # Parallel the bootstrapping
#' boot_cores = 5,
#' out_ipw = TRUE
#'  )
#' print(glm_paths)
#'
#' \dontrun{
#' # Example 2: Only Calculate Pure Imputation Estimator:
#' glm_paths <-
#' pathimp(
#' D = "att22",
#' Y = "std_cesd_age40",
#' M = list("ever_unemp_age3539","log_faminc_adj_age3539"),
#' Y_models = glm_ymodels,
#' data = df,
#' boot_reps = 250,
#' boot_conf_level = 0.95,
#' boot_seed = 2138,
#' boot_parallel = "multicore", # Parallel the bootstrapping
#' boot_cores = 5,
#' out_ipw = FALSE
#'  )
#' print(glm_paths)
#' }

pathimp <- function(
    D,
    Y,
    M,
    Y_models,
    D_model = NULL,
    data,
    boot_reps,
		boot_conf_level = 0.95,
		boot_seed = NULL,
		boot_parallel = "no",
		round_decimal = 3,
		boot_cores = max(c(parallel::detectCores()-2,1)),
    out_ipw){
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
      if(sum(unlist(M) %in% regressors) > 1){
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

  # Set Seed:
  if(!is.null(boot_seed)){
  set.seed(boot_seed)
  }

  # Fit the paths model:
  paths_model <-
    paths(
      a = D,
      y = Y,
      m = M,
      models = Y_models,
      ps_model = D_model,
      data = data,
      nboot = boot_reps,
      conf_level = boot_conf_level,
      ncpus = boot_cores,
      parallel = boot_parallel
    )

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
            out = paste(.data$est,.data$intv)
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
    final_rst <- rbind(processed_result[[1]], processed_result[[2]])
    return(list(summary_df = final_rst, org_obj = paths_model))
  }else{
    return(list(summary_df = processed_result[[1]], org_obj = paths_model))
  }
}
