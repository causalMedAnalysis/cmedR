#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG %dorng%
#' @importFrom foreach foreach
#' @importFrom stats as.formula predict quantile rbinom rmultinom rnorm rpois runif sigma
#' @importFrom MASS polr
#' @importFrom nnet multinom
NULL

#' Print method for medsim objects (point.est)
#'
#' Displays only the point estimates and suppresses the full model output.
#'
#' @param x   An object of class \code{medsim}.
#' @param ...   Other arguments (ignored).
#'
#' @noRd
#' @export
#' @method print medsim
print.medsim <- function(x, ...) {
  est_names <- setdiff(names(x), c("Mmodels", "Ymodel"))
  cat("\u2192 medsim point estimates:\n")
  print(unlist(x[est_names]))
  invisible(x)
}

#' Print method for medsim objects (bootstrap)
#'
#' Displays only the bootstrap results table (point estimates, p-values, and confidence intervals)
#' and suppresses the full model output.
#'
#' @param x   An object of class \code{medsim_boot}, as returned by \code{\link{medsim}(..., boot = TRUE)}.
#' @param ... Other arguments (ignored).
#'
#' @noRd
#' @export
#' @method print medsim_boot
print.medsim_boot <- function(x, ...) {
  cat("\u2192 medsim (bootstrap) results:\n")
  # 'results' is a data.frame with columns: point.est, p.value, ll.*, ul.*
  print(x$results, row.names = TRUE)
  invisible(x)
}

#' Simulate Mediated Effects (Core)
#'
#' This function simulates mediation effects in a causal inference setting. It handles
#' various types of models, including GLMs, multinomial and ordered logistic regression, and allows
#' for the inclusion of survey weights.
#'
#' @param data A data frame containing the variables specified in the model specifications.
#' @param num_sim An integer specifying the number of simulations to run. Default is 2000.
#' @param cat_list A character vector specifying the categories for the treatment variable.
#' @param treatment A string specifying the treatment variable in the dataset.
#' @param intv_med A list specifying the intervention on the mediators. Default is NULL.
#' @param model_spec A list of lists, where each inner list specifies a model for the mediator
#'   or outcome. Each model specification must include the `func` (function name), `formula`,
#'   and `args` (optional arguments).
#' @param weights Optional. A string specifying the name of the column in the data that contains the
#'   weights to be used in weighted regression models and for calculating weighted means. If `NULL`,
#'   weights are not used.
#'
#' @return A list of point estimates for the mediation effects.
#' @noRd
medsim_core <- function(data, num_sim = 2000, cat_list = c("0", "1"), treatment,
                        intv_med = NULL, model_spec, weights = NULL, minimal = FALSE) {

  # Initialize point estimates list
  point.est <- list()

  # Read the dataset
  df <- data

  # Number of mediator models
  num_mediators <- length(model_spec) - 1

  # Initialize list to store mediator models
  Mmodels <- list()

  mediator_mappings <- list()

  # Construct mediator models and extract mediator names
  mediators <- c()

  controlled_value <- list()

  # Define the modify_formula function for intv_med
  modify_formula <- function(formula, mediators) {
    # get the character vector of individual terms
    expanded_terms <- attr(terms(formula), "term.labels")

    # keep only those terms whose all.vars() does _not_ include any mediator name
    keep <- vapply(expanded_terms, function(term) {
      vars_in_term <- all.vars(as.formula(paste("~", term)))
      !any(vars_in_term %in% mediators)
    }, logical(1))

    # rebuild
    lhs <- as.character(formula)[2]
    new_rhs <- paste(expanded_terms[keep], collapse = " + ")
    as.formula(paste(lhs, "~", new_rhs))
  }

  parse_intv_med <- function(intv_med) {
    intv_med_list <- list()
    for (x in intv_med) {
      if (grepl("=", x)) {
        parts <- strsplit(sub("=", "=&=", x), "=&=")[[1]] # splitting on the first equal sign
        intv_med_list[[parts[1]]] <- eval(parse(text=parts[2]))
      } else {
        intv_med_list[[x]] <- NA
      }
    }
    return(intv_med_list)
  }

  intv_med_values <- parse_intv_med(intv_med)

  # Function to convert to factor with integer levels and store mapping
  convert_to_integer_factor <- function(vec) {
    factor_vec <- as.factor(vec)
    levels_map <- levels(factor_vec)

    if (length(levels_map) == 2) {
      if (is.factor(vec) && all(levels_map == c("0", "1"))) {
        levels_map <- NULL
      } else {
        levels(factor_vec) <- c(0, 1)
      }
    } else {
      if (is.factor(vec) && all(levels_map == as.character(1:length(levels_map)))) {
        levels_map <- NULL
      } else {
        levels(factor_vec) <- seq_along(levels_map)
      }
    }

    integer_factor <- as.integer(as.character(factor_vec))
    return(list(factor = integer_factor, levels_map = levels_map))
  }

  # Function to map back from integer factor to original levels
  map_back_to_original_levels <- function(integer_factor, levels_map) {
    if (length(levels_map) == 2) {
      factor(integer_factor, levels = c(0, 1), labels = levels_map)
    } else {
      factor(integer_factor, levels = 1:length(levels_map), labels = levels_map)
    }
  }

  # Process each mediator model
  for (i in 1:num_mediators) {
    mediators[i] <- all.vars(model_spec[[i]]$formula)[1]

    # Check if the family is not NULL and evaluate if it's "binomial"
    fam_arg_mediator <- model_spec[[i]]$args$family

    # pull out a single family name, whether it's a character or a family object
    fam_name_mediator <- NULL
    if (!is.null(fam_arg_mediator)) {
      if (is.character(fam_arg_mediator)) {
        fam_name_mediator <- fam_arg_mediator
      } else if (inherits(fam_arg_mediator, "family")) {
        fam_name_mediator <- fam_arg_mediator$family
      }
    }

    # now test
    is_binomial <- !is.null(fam_name_mediator) && fam_name_mediator == "binomial"

    # Check if the function is multinom or polr
    is_multinom_or_polr <- model_spec[[i]]$func == "multinom" || model_spec[[i]]$func == "polr"

    if (is_binomial || is_multinom_or_polr) {
      # Convert mediator to factor with integer levels and store mapping
      mediator_conversion <- convert_to_integer_factor(df[[mediators[i]]])
      df[[mediators[i]]] <- factor(mediator_conversion$factor)
      mediator_mappings[[mediators[i]]] <- mediator_conversion$levels_map
    }

    if (mediators[i] %in% names(intv_med_values)) {
      if (is.na(intv_med_values[[mediators[i]]])) {
        # Modify formula for intervened mediator
        new_formula <- modify_formula(model_spec[[i]]$formula, mediators)
        Mmodels[[i]] <- do.call(model_spec[[i]]$func, c(list(formula = new_formula, data = df, weights = if (!is.null(weights)) df[[weights]]), model_spec[[i]]$args))
      } else {
        # Handle controlled value by converting it to the new scale

        controlled_value[[i]] <- intv_med_values[[mediators[i]]]

        if (!is.null(mediator_mappings[[mediators[i]]])) {
          # Convert the controlled_value to the new scale using the match function
          new_scale_value <- match(as.character(controlled_value[[i]]), mediator_mappings[[mediators[i]]])

          # Adjust the new_scale_value if there are only two levels (binary case)
          if (length(mediator_mappings[[mediators[i]]]) == 2) {
            new_scale_value <- new_scale_value - 1
          }

          Mmodels[[i]] <- as.factor(new_scale_value)
        } else {
          Mmodels[[i]] <- controlled_value[[i]]
        }
      }
    } else {
      # Fit the model using the original formula, include weights if specified
      Mmodels[[i]] <- do.call(model_spec[[i]]$func, c(list(formula = model_spec[[i]]$formula, data = df, weights = if (!is.null(weights)) df[[weights]]), model_spec[[i]]$args))
    }
  }

  # Convert the outcome variable to factor with integer levels if necessary
  outcome <- all.vars(model_spec[[length(model_spec)]]$formula)[1]

  # Check if the family is not NULL and evaluate if it's "binomial"
  fam_arg_outcome <- model_spec[[length(model_spec)]]$args$family

  # pull out a single family name, whether it's a character or a family object
  fam_name_outcome <- NULL
  if (!is.null(fam_arg_outcome)) {
    if (is.character(fam_arg_outcome)) {
      fam_name_outcome <- fam_arg_outcome
    } else if (inherits(fam_arg_outcome, "family")) {
      fam_name_outcome <- fam_arg_outcome$family
    }
  }

  # now test
  is_binomial  <- !is.null(fam_name_outcome) && fam_name_outcome == "binomial"

  # Check if the function is multinom or polr
  is_multinom_or_polr <- model_spec[[length(model_spec)]]$func == "multinom" || model_spec[[length(model_spec)]]$func == "polr"

  if (is_binomial || is_multinom_or_polr) {
    outcome_conversion <- convert_to_integer_factor(df[[outcome]])
    df[[outcome]] <- factor(outcome_conversion$factor)
    outcome_mappings <- outcome_conversion$levels_map
  }

  # Fit the outcome model, include weights if specified
  Ymodel <- do.call(model_spec[[length(model_spec)]]$func, c(list(formula = model_spec[[length(model_spec)]]$formula, data = df, weights = if (!is.null(weights)) df[[weights]]), model_spec[[length(model_spec)]]$args))

  simulate_draw <- function(model, newdata, model_spec) {
    if (is.numeric(model) || is.factor(model)) {
      return(rep(model, nrow(newdata)))
    } else {
      if (inherits(model, "glm")) {
        phat <- predict(model, newdata = newdata, type = "response")  # Use type = "response" for glm models

        fam <- family(model)$family
        if (fam == "binomial") {
          sim_outcomes <- as.factor(as.character(rbinom(nrow(newdata), size = 1, prob = phat)))
          if (all.vars(model_spec$formula)[1] == outcome && !is.null(outcome_mappings)) {
            mapped_outcomes <- map_back_to_original_levels(as.integer(as.character(sim_outcomes)), outcome_mappings)
            return(mapped_outcomes)
          } else {
            return(sim_outcomes)
          }
        } else if (fam == "quasibinomial") {
          trials <- 100
          return(rbinom(nrow(newdata), size = trials, prob = phat) / trials)
        } else if (fam == "poisson") {
          return(rpois(nrow(newdata), lambda = phat))
        } else if (fam == "gaussian") {
          return(rnorm(nrow(newdata), mean = phat, sd = sigma(model)))
        } else {
          stop("simulate_draw: unsupported glm family '", fam, "'")
        }
      } else if ("lm" %in% class(model) && !("glm" %in% class(model))) {
        phat <- predict(model, newdata = newdata, type = "response")

        return(rnorm(nrow(newdata), mean = phat, sd = sigma(model)))
      } else if (inherits(model, "multinom") || inherits(model, "polr")) {
        phat <- predict(model, newdata = newdata, type = "probs")  # Use type = "probs" for multinom and polr model

        sim_outcomes <- apply(phat, 1, function(prob) {
          rmultinom(n = 1, size = 1, prob = prob)
        })
        sim_outcomes <- as.factor(as.character(apply(sim_outcomes, 2, function(x) which(x == 1))))
        if (all.vars(model_spec$formula)[1] == outcome && !is.null(outcome_mappings)) {
          mapped_outcomes <- map_back_to_original_levels(as.integer(as.character(sim_outcomes)), outcome_mappings)
          return(mapped_outcomes)
        } else {
          return(sim_outcomes)
        }
      } else {
        stop("Unsupported model type")
      }
    }
  }

  if (!is.null(intv_med)) {

    # Initialize vectors to store overall simulation results
    Y_all_treated_intv <- numeric(num_sim)
    Y_all_controlled_intv <- numeric(num_sim)
    Y_intv_controlled <- numeric(num_sim)
    Y_intv_treated <- numeric(num_sim)

    for (i in 1:num_sim) {

      # Initialize dataframes
      idata_intv_treated <- idata_intv_controlled <- idata_all_treated_intv <- idata_all_controlled_intv <- df

      # Check if the treatment variable in the original dataframe is a factor
      if(is.factor(df[[treatment]])) {
        # Convert cat_list values to the same levels as the treatment factor
        idata_all_controlled_intv[[treatment]] <- as.factor(as.character(cat_list[1]))
        idata_all_treated_intv[[treatment]] <- as.factor(as.character(cat_list[2]))
        idata_intv_controlled[[treatment]] <- as.factor(as.character(cat_list[2]))
        idata_intv_treated[[treatment]] <- as.factor(as.character(cat_list[1]))
      } else {
        # If it is numeric, keep cat_list values as numeric
        idata_all_controlled_intv[[treatment]] <- as.numeric(cat_list[1])
        idata_all_treated_intv[[treatment]] <- as.numeric(cat_list[2])
        idata_intv_controlled[[treatment]] <- as.numeric(cat_list[2])
        idata_intv_treated[[treatment]] <- as.numeric(cat_list[1])
      }

      # Iterate over mediators
      for (m in 1:num_mediators) {

        # Simulate for current mediator for OE
        M_simulated_all_treated <- simulate_draw(Mmodels[[m]], idata_all_treated_intv, model_spec[[m]])
        M_simulated_all_controlled <- simulate_draw(Mmodels[[m]], idata_all_controlled_intv, model_spec[[m]])

        if (mediators[m] %in% names(intv_med_values)) {
          # Simulate for current mediator for IE
          M_simulated_intv_treated <- simulate_draw(Mmodels[[m]], idata_intv_controlled, model_spec[[m]])
          M_simulated_intv_controlled <- simulate_draw(Mmodels[[m]], idata_intv_treated, model_spec[[m]])
        } else {
          # Simulate for current mediator for IE
          M_simulated_intv_treated <- simulate_draw(Mmodels[[m]], idata_intv_treated, model_spec[[m]])
          M_simulated_intv_controlled <- simulate_draw(Mmodels[[m]], idata_intv_controlled, model_spec[[m]])
        }

        # Assign simulated values for IE
        idata_intv_treated[[mediators[m]]] <- M_simulated_intv_treated
        idata_intv_controlled[[mediators[m]]] <- M_simulated_intv_controlled

        # Assign simulated values for OE
        idata_all_treated_intv[[mediators[m]]] <- M_simulated_all_treated
        idata_all_controlled_intv[[mediators[m]]] <- M_simulated_all_controlled
      }
      # OE
      Y_all_treated_intv[i] <- if (!is.null(weights)) {
        weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_treated_intv, model_spec[[length(model_spec)]]))), df[[weights]])
      } else {
        mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_treated_intv, model_spec[[length(model_spec)]]))))
      }

      Y_all_controlled_intv[i] <- if (!is.null(weights)) {
        weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_controlled_intv, model_spec[[length(model_spec)]]))), df[[weights]])
      } else {
        mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_controlled_intv, model_spec[[length(model_spec)]]))))
      }

      # IE
      Y_intv_treated[i] <- if (!is.null(weights)) {
        weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_intv_treated, model_spec[[length(model_spec)]]))), df[[weights]])
      } else {
        mean(as.numeric(as.character(simulate_draw(Ymodel, idata_intv_treated, model_spec[[length(model_spec)]]))))
      }

      Y_intv_controlled[i] <- if (!is.null(weights)) {
        weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_intv_controlled, model_spec[[length(model_spec)]]))), df[[weights]])
      } else {
        mean(as.numeric(as.character(simulate_draw(Ymodel, idata_intv_controlled, model_spec[[length(model_spec)]]))))
      }
    }

    # Check if there is at least one numeric value in Mmodels
    if (any(sapply(Mmodels, is.numeric)) || any(sapply(Mmodels, is.factor))) {
      CDE_val <- mean(Y_intv_controlled) - mean(Y_intv_treated)
      # Initialize the base label for CDE
      CDE_label_base <- paste("CDE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1])

      # Gather parts for mediators with numeric values
      numeric_mediators_parts <- c()
      for (m in 1:num_mediators) {
        if (is.numeric(Mmodels[[m]]) || is.factor(Mmodels[[m]])) {
          numeric_mediators_parts <- c(numeric_mediators_parts, paste(",", mediators[m], "=", controlled_value[[m]]))
        }
      }

      # Finalize the label by concatenating all parts
      CDE_label <- paste(CDE_label_base, paste(numeric_mediators_parts, collapse = ","), ")")

      # Assign the CDE value with the constructed label
      point.est[[CDE_label]] <- CDE_val
    } else {
      # Calculate point estimates
      IDE_val <- mean(Y_intv_controlled) - mean(Y_all_controlled_intv)
      IIE_val <- mean(Y_all_treated_intv) - mean(Y_intv_controlled)
      OE_val <- mean(Y_all_treated_intv) - mean(Y_all_controlled_intv)

      point.est[[paste("IDE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")]] <- IDE_val
      point.est[[paste("IIE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")]] <- IIE_val
      point.est[[paste("OE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")]] <- OE_val

      if(is.factor(df[[outcome]])) {
        point.est[paste("E[", outcome, "(", treatment, "=", cat_list[1], ")]")] <- mean(Y_all_controlled_intv)
        point.est[paste("E[", outcome, "(", treatment, "=", cat_list[2], ")]")] <- mean(Y_all_treated_intv)
        point.est[paste("E[(", outcome, "(", treatment, "=", cat_list[2], ",", intv_med, "(", treatment, "=", cat_list[1], "))]")] <- mean(Y_intv_controlled)
        point.est[paste("E[(", outcome, "(", treatment, "=", cat_list[1], ",", intv_med, "(", treatment, "=", cat_list[2], "))]")] <- mean(Y_intv_treated)
      }
    }

  } else {

    # Initialize vectors to store overall simulation results
    Y_all_treated_pse <- numeric(num_sim)
    Y_all_controlled_pse <- numeric(num_sim)

    # Initialize lists of lists for PSE results
    Y_first_treated <- vector("list", num_sim)
    Y_first_controlled <- vector("list", num_sim)
    for (i in 1:num_sim) {
      Y_first_treated[[i]] <- vector("numeric", num_mediators)
      Y_first_controlled[[i]] <- vector("numeric", num_mediators)
    }

    for (i in 1:num_sim) {

      # Initialize dataframes
      idata_all_treated_pse <- idata_all_controlled_pse <- df

      # Check if the treatment variable in the original dataframe is a factor
      if(is.factor(df[[treatment]])) {
        # Convert cat_list values to the same levels as the treatment factor
        idata_all_controlled_pse[[treatment]] <- as.factor(as.character(cat_list[1]))
        idata_all_treated_pse[[treatment]] <- as.factor(as.character(cat_list[2]))
      } else {
        # If it is numeric, keep cat_list values as numeric
        idata_all_controlled_pse[[treatment]] <- as.numeric(cat_list[1])
        idata_all_treated_pse[[treatment]] <- as.numeric(cat_list[2])
      }

      # Initialize lists to store dataframes for each mediator and scenario
      idata_first_controlled <- vector("list", num_mediators)
      idata_first_treated <- vector("list", num_mediators)
      M_simulated_first_treated <- vector("list", num_mediators)
      M_simulated_first_controlled <- vector("list", num_mediators)

      # Iterate over mediators to initialize inner lists
      for (m in 1:num_mediators) {
        M_simulated_first_treated[[m]] <- vector("list", m)
        M_simulated_first_controlled[[m]] <- vector("list", m)
        idata_first_controlled[[m]] <- vector("list", m)
        idata_first_treated[[m]] <- vector("list", m)

        # Initialize vectors in inner lists
        for (j in 1:m) {
          M_simulated_first_treated[[m]][[j]] <- numeric()
          M_simulated_first_controlled[[m]][[j]] <- numeric()
          idata_first_controlled[[m]][[j]] <- numeric()
          idata_first_treated[[m]][[j]] <- numeric()
        }
      }

      # Iterate over mediators
      for (m in 1:num_mediators) {

        # Simulate for current mediator for PSE
        M_simulated_treated_pse <- simulate_draw(Mmodels[[m]], idata_all_treated_pse, model_spec[[m]])
        M_simulated_controlled_pse <- simulate_draw(Mmodels[[m]], idata_all_controlled_pse, model_spec[[m]])

        idata_all_treated_pse[[mediators[m]]] <- M_simulated_treated_pse
        idata_all_controlled_pse[[mediators[m]]] <- M_simulated_controlled_pse

        if (m == 1) {
          idata_first_controlled[[m]][[m]] <- idata_all_treated_pse
          idata_first_treated[[m]][[m]] <- idata_all_controlled_pse

          if(is.factor(df[[treatment]])) {
            # Convert cat_list values to the same levels as the treatment factor
            idata_first_controlled[[m]][[m]][[treatment]] <- as.factor(as.character(cat_list[1]))
            idata_first_treated[[m]][[m]][[treatment]] <- as.factor(as.character(cat_list[2]))
          } else {
            idata_first_controlled[[m]][[m]][[treatment]] <- as.numeric(cat_list[1])
            idata_first_treated[[m]][[m]][[treatment]] <- as.numeric(cat_list[2])
          }
        } else {
          for (j in 2:m) {
            idata_first_controlled[[m]][[j]] <- idata_first_controlled[[m]][[j-1]] <- idata_first_controlled[[m-1]][[j-1]]
            idata_first_treated[[m]][[j]] <- idata_first_treated[[m]][[j-1]] <- idata_first_treated[[m-1]][[j-1]]

            M_simulated_first_treated[[m]][[j-1]] <- simulate_draw(Mmodels[[m]], idata_first_treated[[m-1]][[j-1]], model_spec[[m]])
            M_simulated_first_controlled [[m]][[j-1]]<- simulate_draw(Mmodels[[m]], idata_first_controlled[[m-1]][[j-1]], model_spec[[m]])

            idata_first_controlled[[m]][[j-1]][[mediators[m]]] <- M_simulated_first_controlled[[m]][[j-1]]
            idata_first_treated[[m]][[j-1]][[mediators[m]]] <- M_simulated_first_treated[[m]][[j-1]]
            idata_first_controlled[[m]][[j]][[mediators[m]]] <- M_simulated_treated_pse
            idata_first_treated[[m]][[j]][[mediators[m]]] <- M_simulated_controlled_pse
          }
        }
      }

      # TE
      Y_all_treated_pse[i] <- if (!is.null(weights)) {
        weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_treated_pse, model_spec[[length(model_spec)]]))), df[[weights]])
      } else {
        mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_treated_pse, model_spec[[length(model_spec)]]))))
      }

      Y_all_controlled_pse[i] <- if (!is.null(weights)) {
        weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_controlled_pse, model_spec[[length(model_spec)]]))), df[[weights]])
      } else {
        mean(as.numeric(as.character(simulate_draw(Ymodel, idata_all_controlled_pse, model_spec[[length(model_spec)]]))))
      }

      # PSE
      for (m in 1:num_mediators) {
        Y_first_treated[[i]][[m]] <- if (!is.null(weights)) {
          weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_first_treated[[num_mediators]][[m]], model_spec[[length(model_spec)]]))), df[[weights]])
        } else {
          mean(as.numeric(as.character(simulate_draw(Ymodel, idata_first_treated[[num_mediators]][[m]], model_spec[[length(model_spec)]]))))
        }

        Y_first_controlled[[i]][[m]] <- if (!is.null(weights)) {
          weighted.mean(as.numeric(as.character(simulate_draw(Ymodel, idata_first_controlled[[num_mediators]][[m]], model_spec[[length(model_spec)]]))), df[[weights]])
        } else {
          mean(as.numeric(as.character(simulate_draw(Ymodel, idata_first_controlled[[num_mediators]][[m]], model_spec[[length(model_spec)]]))))
        }
      }
    }


    TE_val <- mean(Y_all_treated_pse) - mean(Y_all_controlled_pse)
    point.est[paste("TE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")] <- TE_val

    if (num_mediators == 1) {
      # PSE for the direct path from treatment to outcome
      point.est[paste("NDE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")] <- mean(sapply(Y_first_treated, function(x) x[[num_mediators]])) - mean(Y_all_controlled_pse)
      # PSE for the all indirect path from treatment to outcome
      point.est[paste("NIE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")] <- mean(Y_all_treated_pse) - mean(sapply(Y_first_treated, function(x) x[[num_mediators]]))
    } else {
      # PSE for the direct path from treatment to outcome
      point.est[paste("PSE(", treatment, "->", outcome, ") or MNDE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")] <- mean(sapply(Y_first_treated, function(x) x[[num_mediators]])) - mean(Y_all_controlled_pse)

      # PSE for the all indirect path from treatment to outcome
      point.est[paste("MNIE(", treatment, "=", cat_list[2], ",", treatment, "*=", cat_list[1], ")")] <- mean(Y_all_treated_pse) - mean(sapply(Y_first_treated, function(x) x[[num_mediators]]))

      # PSE for paths involving mediators
      for (m in 1:num_mediators) {
        if (m == 1) {
          point.est[paste("PSE(", treatment, "->", mediators[[m]], "~>", outcome, ")")] <- mean(Y_all_treated_pse) - mean(sapply(Y_first_treated, function(x) x[[m]]))
        } else if (m < num_mediators){
          point.est[paste("PSE(", treatment, "->", mediators[[m]], "~>", outcome, ")")] <- mean(sapply(Y_first_treated, function(x) x[[m-1]])) - mean(sapply(Y_first_treated, function(x) x[[m]]))
        } else {
          point.est[paste("PSE(", treatment, "->", mediators[[m]], "->", outcome, ")")] <- mean(sapply(Y_first_treated, function(x) x[[m-1]])) - mean(sapply(Y_first_treated, function(x) x[[m]]))
        }
      }
    }

    if(is.factor(df[[outcome]])) {
      point.est[paste("E[", outcome, "(", treatment, "=", cat_list[1], ")]")] <- mean(Y_all_controlled_pse)
      point.est[paste("E[", outcome, "(", treatment, "=", cat_list[2], ")]")] <- mean(Y_all_treated_pse)

      for (m in 1:num_mediators) {
        first_treated_mediator_terms <- paste(lapply(1:m, function(i) paste(mediators[i], "(", treatment, "=",  cat_list[1], ")")), collapse = ", ")
        first_controlled_mediator_terms <- paste(lapply(1:m, function(i) paste(mediators[i], "(", treatment, "=",  cat_list[2], ")")), collapse = ", ")
        point.est[paste("E[(", outcome, "(", treatment, "=", cat_list[2], ",", first_treated_mediator_terms, "))]")] <- mean(sapply(Y_first_treated, function(x) x[[m]]))
        point.est[paste("E[(", outcome, "(", treatment, "=", cat_list[1], ",", first_controlled_mediator_terms, "))]")] <- mean(sapply(Y_first_controlled, function(x) x[[m]]))
      }
    }
  }

  if (minimal) {
    # return only the numeric scalars
    est_names <- setdiff(names(point.est), c("Mmodels","Ymodel"))
    return(unname(unlist(point.est[est_names])))
  }

  # attach full model objects if not minimal
  point.est$Mmodels <- Mmodels
  point.est$Ymodel  <- Ymodel

  class(point.est) <- c("medsim","list")

  return(point.est)
}

#' Simulate Mediated Effects with Bootstrapping
#'
#' This function simulates mediated effects and optionally performs bootstrapping to estimate confidence intervals and p-values.
#'
#' @param data Data frame containing the variables.
#' @param num_sim Number of simulations to run. Default is `2000`.
#' @param cat_list Vector of levels for the treatment. Default is `c("0", "1")`.
#' @param treatment Name of the treatment variable.
#' @param intv_med Intervened mediator. Default is `NULL`.
#' @param model_spec List of model specifications, where each model specification includes the `func`, `formula`, and `args`.
#' @param weights Optional. A string specifying the name of the column in the data that contains the weights to be used in weighted regression models and for calculating weighted means. Default is `NULL`.
#' @param seed Seed for reproducibility. Default is `NULL`.
#' @param boot Logical indicating whether to perform bootstrapping. Default is `FALSE`.
#' @param boot_reps Number of bootstrap replications. Default is 100.
#' @param boot_cores Number of CPU cores for parallel bootstrap. Defaults to available cores minus 2.
#' @param boot_conf_level Confidence level for bootstrap intervals. Default is `0.95`.
#'
#' @return A data frame of point estimates, confidence intervals, and p-values if bootstrapping is performed; otherwise a list of point estimates
#' @export
medsim <- function(data, num_sim = 2000, cat_list = c("0", "1"), treatment,
                   intv_med = NULL, model_spec, weights = NULL, seed = NULL,
                   boot = FALSE, boot_reps = 100, boot_cores = NULL, boot_conf_level = 0.95) {

  # Read the dataset
  df <- data

  # Set a global seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check if bootstrapping is requested
  if (boot) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' needed for this function to work. Please install it.")
    }
    if (!requireNamespace("doRNG", quietly = TRUE)) {
      stop("Package 'doRNG' needed for this function to work. Please install it.")
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'foreach' needed for this function to work. Please install it.")
    }

    # Use boot_cores if provided; otherwise, use available cores minus 2.
    available_cores <- parallel::detectCores()
    if (!is.null(boot_cores)) {
      if (boot_cores > available_cores) {
        stop(paste0("Error: boot_cores (", boot_cores,
                    ") is greater than the available cores (", available_cores, ")."))
      }
      no_cores <- boot_cores
    } else {
      no_cores <- available_cores - 2
    }

    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)

    # Export the medsim_core function and any other necessary objects to each worker
    parallel::clusterExport(cl, varlist = c("medsim_core"), envir=environment())

    # Use doRNG::registerDoRNG to ensure reproducibility
    doRNG::registerDoRNG(seed)

    # Initialize a list to store bootstrap results
    `%dorng%` <- doRNG::`%dorng%`

    bootstrap_results <- foreach::foreach(n = 1:boot_reps, .combine = cbind, .packages = c("dplyr", "MASS", "nnet")) %dorng% {
      boot.df <- df[sample(nrow(df), nrow(df), replace = TRUE), ]
      medsim_core(data = boot.df, num_sim = num_sim, cat_list = cat_list, treatment = treatment,
                  intv_med = intv_med, model_spec = model_spec, weights = weights, minimal = TRUE)
    }

    parallel::stopCluster(cl)

    # Restructure the results
    medsim.boot <- matrix(unlist(bootstrap_results), ncol=nrow(bootstrap_results), byrow=TRUE)

    # Calculate bootstrap confidence intervals based on boot_conf_level
    lower_prob <- (1 - boot_conf_level) / 2
    upper_prob <- 1 - lower_prob
    ci_limits <- apply(medsim.boot, 2, function(x) stats::quantile(x, probs = c(lower_prob, upper_prob)))

    # Calculate p-values
    p_values <- apply(medsim.boot, 2, function(x) {
      ll <- mean(x > 0)
      ul <- mean(x < 0)
      p_val <- 2 * min(ll, ul)
      return(p_val)
    })

    # Using medsim_core for point estimate
    core_out <- medsim_core(data = df, num_sim = num_sim, cat_list = cat_list, treatment = treatment,
                                   intv_med = intv_med, model_spec = model_spec, weights = weights)

    # Extract just the numeric estimates
    est_names <- setdiff(names(core_out), c("Mmodels", "Ymodel"))
    point.est <- unlist(core_out[est_names])

    # Round the estimates, p-values and confidence intervals to 3 digits
    point.est <- round(point.est, 3)
    p_values  <- round(p_values, 3)
    ci_limits <- round(ci_limits, 3)

    # Dynamically create CI labels
    ll_label <- paste0("ll.", boot_conf_level * 100, "ci")
    ul_label <- paste0("ul.", boot_conf_level * 100, "ci")

    results <- data.frame(point.est = point.est,
                    p.value = p_values)
    results[[ll_label]] <- ci_limits[1, ]
    results[[ul_label]] <- ci_limits[2, ]

    row.names(results) <- est_names

    out <- list(
      results = results,
      Mmodels = core_out$Mmodels,
      Ymodel  = core_out$Ymodel
    )

    class(out) <- c("medsim_boot","list")

    return(out)
  } else {
    return(medsim_core(data = df, num_sim = num_sim, cat_list = cat_list, treatment = treatment,
                       intv_med = intv_med, model_spec = model_spec, weights = weights))
  }
}
