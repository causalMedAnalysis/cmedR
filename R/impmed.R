#' Regression-Imputation Estimator for Natural Direct and Indirect Effects
#'
#' @description
#' `impmed()` estimates the total effect and the natural direct and indirect
#' effects by pure regression imputation using a simplified version of the
#' `pathimp()` workflow. The user supplies two fitted outcome models:
#'
#' 1. A reduced outcome model for `E(Y | D, C)`.
#' 2. A full outcome model for `E(Y | D, M, C)`, where `M` may contain one or
#'    more mediators.
#'
#' @details
#' Relative to `pathimp()`, `impmed()` treats all mediators as a single block.
#' It does not estimate mediator-specific path-specific effects and it does not
#' implement the imputation-based weighting estimator. Instead, it reproduces
#' the simpler regression-imputation logic:
#'
#' 1. Use the reduced model to impute `E{Y(d, M(d)) | C}` and
#'    `E{Y(dstar, M(dstar)) | C}`.
#' 2. Use the full model to predict `Y` under `D = d` while leaving the
#'    observed mediators unchanged.
#' 3. Regress those predictions on the reduced-model regressors to recover
#'    `E{Y(d, M(dstar)) | C}`.
#' 4. Average the three potential-outcome means and take contrasts to obtain the
#'    total effect, natural direct effect, and natural indirect effect.
#'
#' If more than one mediator is supplied, `impmed()` reports the joint natural
#' direct and indirect effects, labeled `MNDE` and `MNIE` to mirror the Stata
#' implementation.
#'
#' Treatment-mediator, treatment-confounder, and confounder-mediator
#' interactions should be specified directly in the model formulas when
#' needed. Interaction terms are handled automatically by `predict()`
#' based on the fitted model.
#'
#' Confidence intervals are computed with a nonparametric percentile bootstrap.
#' In each bootstrap replication, both user-supplied outcome models are refit on
#' the resampled analysis data, the bridge regression is refit, and all effects
#' are recalculated.
#'
#' @param D A character scalar naming the treatment variable in `data`. The
#' treatment variable should be numeric because the estimator intervenes on it
#' using the numeric values `d` and `dstar`.
#' @param Y A character scalar naming the outcome variable in `data`.
#' @param M A character vector, or a list of character vectors, naming one or
#' more mediators. All mediators are treated jointly rather than as causally
#' ordered mediator blocks.
#' @param Y_models A list of exactly two fitted outcome models. The first
#' element must be the reduced outcome model, excluding the mediators in `M`.
#' The second element must be the full outcome model, including all mediators in
#' `M`. Currently, `impmed()` supports either two fitted `lm` objects or two
#' fitted `glm` objects. For `glm`, both models must use the same family and
#' link, and the supported specifications are `gaussian(identity)` and
#' `binomial(logit)`.
#' @param data A data frame used to fit the supplied models. This argument is
#' retained for API consistency with the rest of the package. Internally,
#' `impmed()` recovers the common analysis sample from the fitted models and
#' requires that both models were fit on the same rows.
#' @param d A numeric scalar giving the treatment value for the first
#' intervention.
#' @param dstar A numeric scalar giving the treatment value for the comparison
#' intervention.
#' @param boot A logical scalar indicating whether the function will perform the
#' nonparametric bootstrap and return two-sided confidence intervals and
#' p-values. Defaults to FALSE, in which case only point estimates are returned.
#' @param boot_reps An integer scalar giving the number of bootstrap
#' replications. In practice, we recommend a minimum of 1000 replications. Only
#' consulted when `boot = TRUE`.
#' @param boot_conf_level A numeric scalar giving the bootstrap confidence
#' level. Defaults to `0.95`.
#' @param boot_seed An optional integer scalar specifying the random-number seed
#' used for bootstrap resampling. Defaults to `NULL`, meaning that no seed is
#' set.
#' @param boot_parallel A logical scalar indicating whether the bootstrap will
#' be performed with a parallelized loop, with the goal of reducing runtime.
#' Defaults to FALSE. Parallelization uses the `"snow"` backend of
#' [boot::boot()], which works across all operating systems.
#' @param round_decimal An integer scalar giving the number of decimal places
#' used in the printed summary table.
#' @param boot_cores An integer scalar specifying the number of CPU cores on
#' which the parallelized bootstrap will run, passed as `ncpus` to
#' [boot::boot()]. This argument only has an effect if you requested a
#' parallelized bootstrap (i.e., only if `boot` is TRUE and `boot_parallel` is
#' TRUE). By default, `boot_cores` is equal to the greater of two values: (a)
#' one and (b) the number of available CPU cores minus two. If `boot_cores`
#' equals one, then the bootstrap loop will not be parallelized (regardless of
#' whether `boot_parallel` is TRUE).
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{summary_df}{A data frame with the formatted estimates for the ATE,
#'   NDE or MNDE, and NIE or MNIE. When `boot = TRUE`, the estimate column also
#'   carries percentile bootstrap confidence intervals and the data frame gains
#'   a column of two-sided p-values.}
#'   \item{org_obj}{A list containing the numeric estimates, the implied
#'   potential-outcome means, and the intervention values. When `boot = TRUE`,
#'   it additionally contains the confidence intervals, p-values, bootstrap
#'   draws, and bootstrap settings.}
#' }
#'
#' @examples
#' \dontrun{
#' data(Brader)
#'
#' ## Prepare data
#' Y <- "std_immigr"
#' D <- "tone_eth"
#' C <- c("ppage", "female", "hs", "sc", "ba", "ppincimp")
#'
#' # number of bootstrap replications
#' nboot <- 1000
#'
#' # set seed
#' boot_seed <- 1234
#'
#' # Example 1: Single mediator
#' M1 <- "p_harm"
#' key_vars <- c("immigr", D, M, C)
#'
#' Brader1 <- Brader[complete.cases(Brader[, key_vars]), ]
#' Brader1$std_immigr <-
#'   (Brader1$immigr - mean(Brader1$immigr)) / sd(Brader1$immigr)
#'
#' glm_reduced <- glm(
#'   as.formula(paste(Y, "~", D, "+", paste(C, collapse = " + "))),
#'   data = Brader1
#' )
#'
#' glm_full <- glm(
#'   as.formula(paste(Y, "~", D, "+", M1, "+", paste(C, collapse = " + "))),
#'   data = Brader1
#')
#'
#' impmed(
#'   D = D,
#'   Y = Y,
#'   M = M,
#'   Y_models = list(glm_reduced, glm_full),
#'   data = Brader1,
#'   d = 1,
#'   dstar = 0,
#'   boot_reps = nboot,
#'   boot_conf_level = 0.95,
#'   boot_seed = boot_seed,
#'   boot_parallel = "multicore"
#' )$summary_df
#'
#' # Example 2: Two mediators
#' M <- c("p_harm", "emo")
#' key_vars <- c("immigr", D, M, C)
#'
#' Brader2 <- Brader[complete.cases(Brader[, key_vars]), ]
#' Brader2$std_immigr <-
#'   (Brader2$immigr - mean(Brader2$immigr)) / sd(Brader2$immigr)
#'
#' glm_reduced2 <- glm(
#'   as.formula(paste(Y, "~", D, "+", paste(C, collapse = " + "))),
#'   data = Brader2
#' )
#'
#' glm_full2 <- glm(
#'   as.formula(paste(Y, "~", D, "+", M, "+", paste(C, collapse = " + "))),
#'   data = Brader2
#')
#' impmed(
#'   D = D,
#'   Y = Y,
#'   M = M,
#'   Y_models = list(glm_reduced2, glm_full2),
#'   data = Brader2,
#'   d = 1,
#'   dstar = 0,
#'   boot_reps = nboot,
#'   boot_conf_level = 0.95,
#'   boot_seed = boot_seed,
#'   boot_parallel = "multicore"
#' )$summary_df
#' }
#'
#' @importFrom boot boot
#' @export
impmed <- function(
    D,
    Y,
    M,
    Y_models,
    data,
    d = 1,
    dstar = 0,
    boot = FALSE,
    boot_reps = 200,
    boot_conf_level = 0.95,
    boot_seed = NULL,
    boot_parallel = FALSE,
    round_decimal = 3,
    boot_cores = max(c(parallel::detectCores() - 2L, 1L), na.rm = TRUE)) {

  mediator_names <- .impmed_normalize_mediators(M = M)

  .impmed_validate_inputs(
    D = D,
    Y = Y,
    M = mediator_names,
    Y_models = Y_models,
    data = data,
    d = d,
    dstar = dstar,
    boot_reps = boot_reps,
    boot_conf_level = boot_conf_level,
    boot_parallel = boot_parallel,
    round_decimal = round_decimal,
    boot_cores = boot_cores
  )

  # recover the common analysis sample used by both supplied outcome models
  analysis_data <- .impmed_get_analysis_data(
    Y_models = Y_models,
    data = data,
    Y = Y,
    D = D,
    M = mediator_names
  )

  # compute the point estimates using the user-supplied reduced and full
  # outcome models plus an internal bridge regression for the cross-world mean
  point_out <- .impmed_point_estimates(
    reduced_model = Y_models[[1L]],
    full_model = Y_models[[2L]],
    data = analysis_data,
    D = D,
    d = d,
    dstar = dstar,
    effect_names = .impmed_effect_names(length(mediator_names))
  )

  effect_names <- names(point_out$effects)

  estimand_labels <- paste0(
    effect_names,
    "(",
    .impmed_format_value(d, round_decimal),
    ", ",
    .impmed_format_value(dstar, round_decimal),
    ")"
  )

  if (boot) {
    # "snow" is the only cross-platform boot backend ("multicore" fails on Windows)
    boot_parallel_use <- if (isTRUE(boot_parallel) && boot_cores > 1L) {
      "snow"
    } else {
      "no"
    }

    if (!is.null(boot_seed)) {
      # L'Ecuyer-CMRG lets boot distribute reproducible RNG streams to snow workers
      if (boot_parallel_use == "snow") {
        old_rng_kind <- RNGkind("L'Ecuyer-CMRG")[[1L]]
        on.exit(RNGkind(old_rng_kind), add = TRUE)
      }
      set.seed(boot_seed)
    }

    # snow workers are bare R processes, so build the cluster ourselves and
    # export the internal helpers the bootstrap statistic depends on
    boot_cluster <- NULL
    if (boot_parallel_use == "snow") {
      boot_cluster <- parallel::makePSOCKcluster(boot_cores)
      on.exit(parallel::stopCluster(boot_cluster), add = TRUE)
      helper_env <- topenv(environment())
      parallel::clusterExport(
        boot_cluster,
        varlist = ls(helper_env, all.names = TRUE, pattern = "^\\.impmed_"),
        envir = helper_env
      )
    }

    boot_error_env <- new.env(parent = emptyenv())
    boot_error_env$messages <- character(0)
    boot_error_env$counter <- 0L

    # refit both outcome models on each bootstrap sample and recompute the full
    # set of estimands
    boot_out <- boot::boot(
      data = analysis_data,
      statistic = function(data, indices) {
        sampled_data <- .impmed_prepare_boot_sample(data = data, indices = indices)
        rep_id <- .impmed_next_boot_rep_id(
          boot_error_env = boot_error_env,
          parallel = boot_parallel_use
        )

        refit_models <- tryCatch(
          .impmed_refit_models(Y_models = Y_models, data = sampled_data),
          error = function(e) {
            .impmed_record_boot_issue(
              boot_error_env = boot_error_env,
              stage = "model refit",
              message_text = conditionMessage(e),
              rep_id = rep_id,
              parallel = boot_parallel_use
            )
            NULL
          }
        )

        if (is.null(refit_models)) {
          return(rep(NA_real_, 3L))
        }

        boot_point_out <- tryCatch(
          .impmed_point_estimates(
            reduced_model = refit_models[[1L]],
            full_model = refit_models[[2L]],
            data = sampled_data,
            D = D,
            d = d,
            dstar = dstar,
            effect_names = .impmed_effect_names(length(mediator_names))
          ),
          error = function(e) {
            .impmed_record_boot_issue(
              boot_error_env = boot_error_env,
              stage = "point estimation",
              message_text = conditionMessage(e),
              rep_id = rep_id,
              parallel = boot_parallel_use
            )
            NULL
          }
        )

        if (is.null(boot_point_out)) {
          return(rep(NA_real_, 3L))
        }

        boot_estimates <- as.numeric(boot_point_out$effects)
        if (length(boot_estimates) != 3L || any(!is.finite(boot_estimates))) {
          .impmed_record_boot_issue(
            boot_error_env = boot_error_env,
            stage = "point estimation",
            message_text = "returned non-finite effect estimates.",
            rep_id = rep_id,
            parallel = boot_parallel_use
          )
          return(rep(NA_real_, 3L))
        }

        return(as.numeric(boot_estimates))
      },
      R = boot_reps,
      parallel = boot_parallel_use,
      ncpus = if (boot_parallel_use == "no") 1L else boot_cores,
      cl = boot_cluster
    )

    boot_draws <- as.matrix(boot_out$t)
    colnames(boot_draws) <- effect_names
    successful_boot <- apply(boot_draws, 1L, function(x) all(is.finite(x)))
    n_boot_fail <- sum(!successful_boot)
    boot_draws <- boot_draws[successful_boot, , drop = FALSE]

    if (nrow(boot_draws) < 2L) {
      stop(
        "Too few successful bootstrap replications were available to compute the ",
        "confidence intervals.",
        .impmed_format_boot_issue_summary(
          messages = boot_error_env$messages,
          parallel = boot_parallel_use
        )
      )
    }

    alpha <- 1 - boot_conf_level
    conf_int <- t(apply(
      boot_draws,
      2L,
      stats::quantile,
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE,
      names = FALSE
    ))
    colnames(conf_int) <- c("lower", "upper")

    # compute two-sided bootstrap p-values (CI-inversion method)
    pvalues <- apply(boot_draws, 2L, function(x) 2 * min(mean(x < 0), mean(x > 0)))

    if (n_boot_fail > 0L) {
      warning(
        n_boot_fail,
        " bootstrap replication(s) failed and were omitted from the percentile ",
        "confidence intervals.",
        .impmed_format_boot_issue_summary(
          messages = boot_error_env$messages,
          parallel = boot_parallel_use
        )
      )
    }

    summary_df <- data.frame(
      estimator = rep("Pure Imputation Estimator", length(effect_names)),
      estimand = estimand_labels,
      out = paste0(
        round(point_out$effects, round_decimal),
        " [",
        round(conf_int[, "lower"], round_decimal),
        ", ",
        round(conf_int[, "upper"], round_decimal),
        "]"
      ),
      pvalue = round(pvalues, round_decimal),
      stringsAsFactors = FALSE
    )

    org_obj <- list(
      call = match.call(),
      estimates = point_out$effects,
      conf_int = conf_int,
      pvalues = pvalues,
      potential_outcome_means = point_out$means,
      boot_draws = boot_draws,
      boot_failures = n_boot_fail,
      boot_conf_level = boot_conf_level,
      boot_reps = boot_reps,
      intervention_values = list(d = d, dstar = dstar),
      mediator_names = mediator_names,
      analysis_n = nrow(analysis_data),
      model_class = .impmed_get_model_class(Y_models[[1L]]),
      model_family = .impmed_get_model_family(Y_models[[1L]]),
      boot_object = boot_out
    )
  } else {
    # point estimates only: mirror the slimmer return shape used by the other
    # estimators (e.g. linmed, ipwmed, mrmed) when boot = FALSE, omitting the
    # confidence-interval and p-value fields entirely rather than emitting blank
    # placeholders.
    summary_df <- data.frame(
      estimator = rep("Pure Imputation Estimator", length(effect_names)),
      estimand = estimand_labels,
      out = as.character(round(point_out$effects, round_decimal)),
      stringsAsFactors = FALSE
    )

    org_obj <- list(
      call = match.call(),
      estimates = point_out$effects,
      potential_outcome_means = point_out$means,
      intervention_values = list(d = d, dstar = dstar),
      mediator_names = mediator_names,
      analysis_n = nrow(analysis_data),
      model_class = .impmed_get_model_class(Y_models[[1L]]),
      model_family = .impmed_get_model_family(Y_models[[1L]])
    )
  }

  return(list(summary_df = summary_df, org_obj = org_obj))
}


# convert the mediator specification to a flat character vector
.impmed_normalize_mediators <- function(M) {
  if (is.character(M)) {
    mediator_names <- M
  } else if (is.list(M) && all(vapply(M, is.character, logical(1)))) {
    mediator_names <- unlist(M, use.names = FALSE)
  } else {
    stop("`M` must be a character vector or a list of character vectors.")
  }

  mediator_names <- unname(mediator_names)
  if (length(mediator_names) < 1L) {
    stop("`M` must contain at least one mediator name.")
  }

  if (anyNA(mediator_names) || any(mediator_names == "")) {
    stop("`M` must contain non-missing mediator names.")
  }

  if (anyDuplicated(mediator_names)) {
    stop("`M` must not contain duplicate mediator names.")
  }

  mediator_names
}


# validate the public arguments and both supplied fitted outcome models
.impmed_validate_inputs <- function(
    D,
    Y,
    M,
    Y_models,
    data,
    d,
    dstar,
    boot_reps,
    boot_conf_level,
    boot_parallel,
    round_decimal,
    boot_cores) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  for (var_name in c("D", "Y")) {
    value <- get(var_name)
    if (!is.character(value) || length(value) != 1L) {
      stop("`", var_name, "` must be a character scalar naming a variable.")
    }
    if (!value %in% names(data)) {
      stop("`", var_name, "` must name a column in `data`.")
    }
  }

  if (!all(M %in% names(data))) {
    stop("Every mediator named in `M` must be present in `data`.")
  }

  if (!is.numeric(data[[D]])) {
    stop("`D` must refer to a numeric treatment variable.")
  }

  for (arg_name in c("d", "dstar")) {
    value <- get(arg_name)
    if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
      stop("`", arg_name, "` must be a non-missing numeric scalar.")
    }
  }

  if (!is.list(Y_models) || length(Y_models) != 2L) {
    stop("`Y_models` must be a list of exactly two fitted outcome models.")
  }

  model_classes <- vapply(Y_models, .impmed_get_model_class, character(1))
  if (!all(model_classes %in% c("lm", "glm"))) {
    stop("Each element of `Y_models` must be a fitted `lm` or `glm` object.")
  }

  if (length(unique(model_classes)) != 1L) {
    stop("Both elements of `Y_models` must have the same model class.")
  }

  model_families <- vapply(Y_models, .impmed_get_model_family, character(1))
  model_links <- vapply(Y_models, .impmed_get_model_link, character(1))

  if (model_classes[[1L]] == "glm") {
    if (length(unique(model_families)) != 1L || length(unique(model_links)) != 1L) {
      stop("Both `glm` objects in `Y_models` must use the same family and link.")
    }

    family_link <- paste0(model_families[[1L]], "(", model_links[[1L]], ")")
    if (!family_link %in% c("gaussian(identity)", "binomial(logit)")) {
      stop(
        "Supported `glm` specifications are `gaussian(identity)` and ",
        "`binomial(logit)`."
      )
    }
  }

  if (!is.numeric(boot_reps) || length(boot_reps) != 1L ||
      boot_reps < 2L || boot_reps != as.integer(boot_reps)) {
    stop("`boot_reps` must be a single integer greater than or equal to 2.")
  }

  if (!is.numeric(boot_conf_level) || length(boot_conf_level) != 1L ||
      boot_conf_level <= 0 || boot_conf_level >= 1) {
    stop("`boot_conf_level` must be a single number between 0 and 1.")
  }

  if (!is.logical(boot_parallel) || length(boot_parallel) != 1L ||
      is.na(boot_parallel)) {
    stop("`boot_parallel` must be a single logical value (TRUE or FALSE).")
  }

  if (!is.numeric(round_decimal) || length(round_decimal) != 1L ||
      round_decimal < 0 || round_decimal != as.integer(round_decimal)) {
    stop("`round_decimal` must be a single non-negative integer.")
  }

  if (!is.numeric(boot_cores) || length(boot_cores) != 1L || boot_cores < 1) {
    stop("`boot_cores` must be a single positive number.")
  }

  reduced_terms <- all.vars(stats::formula(Y_models[[1L]]))
  full_terms <- all.vars(stats::formula(Y_models[[2L]]))

  if (!Y %in% reduced_terms || !Y %in% full_terms) {
    stop("Both supplied outcome models must use `Y` as the outcome variable.")
  }

  if (!D %in% reduced_terms || !D %in% full_terms) {
    stop("Both supplied outcome models must include the treatment variable `D`.")
  }

  if (any(M %in% reduced_terms)) {
    stop(
      "The first element of `Y_models` must be the reduced model and must not ",
      "include any mediators named in `M`."
    )
  }

  if (!all(M %in% full_terms)) {
    stop(
      "The second element of `Y_models` must be the full model and must include ",
      "all mediators named in `M`."
    )
  }
}


# recover the common analysis sample from the two fitted models and return the
# corresponding rows from the original data frame so raw variables remain
# available during bootstrap refits
.impmed_get_analysis_data <- function(Y_models, data, Y, D, M) {
  model_frames <- lapply(
    Y_models,
    function(model) {
      tryCatch(
        stats::model.frame(model),
        error = function(e) NULL
      )
    }
  )

  if (any(vapply(model_frames, is.null, logical(1)))) {
    stop(
      "Unable to recover the estimation sample from one or both outcome models. ",
      "Please fit the models with a standard `data =` argument so that the ",
      "model frames can be reconstructed."
    )
  }

  sample_rows <- lapply(model_frames, rownames)
  if (!identical(sample_rows[[1L]], sample_rows[[2L]])) {
    stop(
      "The two outcome models in `Y_models` must be fit on the same analysis ",
      "sample. In practice this usually means fitting both models on the same ",
      "mediator-complete data set before calling `impmed()`."
    )
  }

  analysis_rows <- sample_rows[[1L]]
  if (length(analysis_rows) < 1L) {
    stop("The recovered analysis sample is empty.")
  }

  if (is.null(rownames(data)) || any(!analysis_rows %in% rownames(data))) {
    stop(
      "Unable to align the fitted model samples with `data`. Please supply the ",
      "same data frame that was used to fit both outcome models."
    )
  }

  analysis_data <- data[analysis_rows, , drop = FALSE]
  needed_vars <- c(Y, D, M)
  if (!all(needed_vars %in% names(analysis_data))) {
    stop(
      "The recovered analysis data do not contain `Y`, `D`, and every mediator ",
      "named in `M` as raw columns."
    )
  }

  if (nrow(analysis_data) < 1L) {
    stop("The recovered analysis sample is empty.")
  }

  if (inherits(Y_models[[1L]], "glm") &&
      identical(stats::family(Y_models[[1L]])$family, "binomial")) {
    observed_outcomes <- analysis_data[[Y]]
    if (!all(stats::na.omit(observed_outcomes) %in% c(0, 1))) {
      stop(
        "For `binomial(logit)` models, the observed outcome `Y` must be coded ",
        "as 0/1 on the analysis sample."
      )
    }
  }

  analysis_data
}


# refit both user-supplied outcome models on a bootstrap sample. Refit from the
# evaluated formula so the refit needs only `data`, not the caller's variables.
.impmed_refit_models <- function(Y_models, data) {
  lapply(
    Y_models,
    function(model) {
      model_formula <- stats::formula(model)
      environment(model_formula) <- environment()
      if (inherits(model, "glm")) {
        stats::glm(model_formula, family = stats::family(model), data = data)
      } else {
        stats::lm(model_formula, data = data)
      }
    }
  )
}


# Standardize bootstrap samples to plain data frames with clean row names so the
# refit path does not depend on how the original analysis data stored them.
.impmed_prepare_boot_sample <- function(data, indices) {
  sampled_data <- as.data.frame(data[indices, , drop = FALSE])
  rownames(sampled_data) <- NULL
  sampled_data
}


# predict the conditional mean of the outcome on the response scale
.impmed_predict_response <- function(model, newdata) {
  if (inherits(model, "glm")) {
    return(as.numeric(stats::predict(model, newdata = newdata, type = "response")))
  }

  return(as.numeric(stats::predict(model, newdata = newdata)))
}


# create a counterfactual data set in which treatment is set to a common value
# for every observation.
.impmed_set_treatment <- function(data, D, d_value) {
  cf_data <- data
  cf_data[[D]] <- d_value
  cf_data
}


# build and fit the bridge regression used to transport the full-model
# predictions from the observed treatment distribution to `D = dstar`
.impmed_fit_bridge_model <- function(reduced_model, data, pseudo_outcome) {
  bridge_data <- data
  bridge_data[[".impmed_pseudo_outcome"]] <- pseudo_outcome

  bridge_formula <- stats::formula(reduced_model)
  bridge_formula[[2L]] <- as.name(".impmed_pseudo_outcome")

  if (inherits(reduced_model, "glm")) {
    model_family <- stats::family(reduced_model)

    # R's binomial glm warns when the pseudo-outcome is fractional. Using
    # quasibinomial preserves the same mean-model estimating equations while
    # avoiding spurious bootstrap warnings
    if (identical(model_family$family, "binomial")) {
      model_family <- stats::quasibinomial(link = model_family$link)
    }

    return(stats::glm(
      formula = bridge_formula,
      family = model_family,
      data = bridge_data
    ))
  }

  stats::lm(formula = bridge_formula, data = bridge_data)
}


# compute the three imputed potential-outcome means and convert them into the
# total, natural direct, and natural indirect effects
.impmed_point_estimates <- function(
    reduced_model,
    full_model,
    data,
    D,
    d,
    dstar,
    effect_names) {
  data_d <- .impmed_set_treatment(data = data, D = D, d_value = d)
  data_dstar <- .impmed_set_treatment(data = data, D = D, d_value = dstar)

  pred_d_md <- .impmed_predict_response(model = reduced_model, newdata = data_d)
  .impmed_assert_finite_vector(
    x = pred_d_md,
    label = "Reduced-model predictions under D = d"
  )
  mu_d_md <- mean(pred_d_md)

  pred_dstar_mdstar <- .impmed_predict_response(
    model = reduced_model,
    newdata = data_dstar
  )
  .impmed_assert_finite_vector(
    x = pred_dstar_mdstar,
    label = "Reduced-model predictions under D = dstar"
  )
  mu_dstar_mdstar <- mean(pred_dstar_mdstar)

  pseudo_outcome <- .impmed_predict_response(model = full_model, newdata = data_d)
  .impmed_assert_finite_vector(
    x = pseudo_outcome,
    label = "Full-model predictions under D = d"
  )
  bridge_model <- .impmed_fit_bridge_model(
    reduced_model = reduced_model,
    data = data,
    pseudo_outcome = pseudo_outcome
  )
  pred_d_mdstar <- .impmed_predict_response(model = bridge_model, newdata = data_dstar)
  .impmed_assert_finite_vector(
    x = pred_d_mdstar,
    label = "Bridge-model predictions under D = dstar"
  )
  mu_d_mdstar <- mean(pred_d_mdstar)

  .impmed_assert_finite_vector(
    x = c(mu_d_md, mu_dstar_mdstar, mu_d_mdstar),
    label = "Imputed potential-outcome means"
  )

  effects <- c(
    ATE = mu_d_md - mu_dstar_mdstar,
    NDE = mu_d_mdstar - mu_dstar_mdstar,
    NIE = mu_d_md - mu_d_mdstar
  )
  names(effects) <- effect_names
  .impmed_assert_finite_vector(
    x = effects,
    label = "Implied causal effect estimates"
  )

  return(list(
    effects = effects,
    means = c(
      mu_d_md = mu_d_md,
      mu_dstar_mdstar = mu_dstar_mdstar,
      mu_d_mdstar = mu_d_mdstar
    )
  ))
}


#identify the broad model class while respecting that `glm` also inherits from lm.
.impmed_get_model_class <- function(model) {
  if (inherits(model, "glm")) {
    return("glm")
  }

  if (inherits(model, "lm")) {
    return("lm")
  }

  return(NA_character_)
}


# name the direct and indirect effects differently depending on whether the
# mediators are univariate or multivariate
.impmed_effect_names <- function(n_mediators) {
  c(
    "ATE",
    if (n_mediators == 1L) "NDE" else "MNDE",
    if (n_mediators == 1L) "NIE" else "MNIE"
  )
}


# extract the model family in a way that works for both `lm` and `glm`
.impmed_get_model_family <- function(model) {
  if (inherits(model, "glm")) {
    return(stats::family(model)$family)
  }

  "gaussian"
}


# extract the link function in a way that works for both `lm` and `glm`
.impmed_get_model_link <- function(model) {
  if (inherits(model, "glm")) {
    return(stats::family(model)$link)
  }

  "identity"
}


# Fail loudly when a prediction or effect vector contains NA, NaN, or Inf so
# the bootstrap branch can report the real issue instead of silently returning
# all-NA replications.
.impmed_assert_finite_vector <- function(x, label) {
  if (length(x) < 1L || any(!is.finite(x))) {
    stop(label, " contained non-finite values.")
  }
}


# Track bootstrap replication order when bootstrapping runs in the main R
# process. Parallel backends may not share this state, so the identifier is
# optional in that case.
.impmed_next_boot_rep_id <- function(boot_error_env, parallel) {
  if (identical(parallel, "no")) {
    boot_error_env$counter <- boot_error_env$counter + 1L
    return(boot_error_env$counter)
  }

  NA_integer_
}


# Store and echo bootstrap failures so the first real errors are visible during
# debugging instead of being swallowed into NA return values.
.impmed_record_boot_issue <- function(
    boot_error_env,
    stage,
    message_text,
    rep_id,
    parallel) {
  prefix <- if (is.na(rep_id)) {
    "Bootstrap replication"
  } else {
    paste("Bootstrap replication", rep_id)
  }

  issue_text <- paste0(prefix, " failed during ", stage, ": ", message_text)
  boot_error_env$messages <- unique(c(boot_error_env$messages, issue_text))

  if (length(boot_error_env$messages) <= 10L || identical(parallel, "no")) {
    message(issue_text)
  }
}


# Summarize the first few captured bootstrap failures in the final warning or
# error message.
.impmed_format_boot_issue_summary <- function(messages, parallel) {
  if (length(messages) < 1L) {
    if (!identical(parallel, "no")) {
      return(
        paste0(
          " No bootstrap error details were captured in the main R process. ",
          "If you are debugging a parallel bootstrap, rerun with ",
          "`boot_parallel = FALSE` to surface the per-replication error ",
          "messages directly."
        )
      )
    }

    return("")
  }

  n_show <- min(length(messages), 3L)
  paste0(
    " First bootstrap issue(s): ",
    paste(messages[seq_len(n_show)], collapse = " | ")
  )
}


# format intervention values cleanly in the summary table without changing the
# underlying numeric values stored in `org_obj`
.impmed_format_value <- function(x, digits) {
  format(round(x, digits), trim = TRUE, scientific = FALSE)
}
