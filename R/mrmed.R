#' Multiply robust (MR) estimator for total, natural direct and natural
#' indirect effects - inner function
#'
#' @description
#' Internal function used within `mrmed()`. See the `mrmed()` function
#' documentation for a description of shared function arguments. Here, we will
#' only document the one argument that is not shared by `mrmed_inner()` and
#' `mrmed()`: the `minimal` argument.
#'
#' @param minimal A logical scalar indicating whether the function should
#'   return only a minimal set of output. The `mrmed()` function uses the
#'   default of FALSE when calling `mrmed_inner()` to generate the point
#'   point estimates and sets the argument to TRUE when calling `mrmed_inner()`
#'   to perform the bootstrap.
#'
#' @noRd
#---------------------- The Inner Function ----------------------#
mrmed_inner <- function(
    D,
    Y,
    M,
    C,
    D_C_model, # D ~ C
    D_MC_model = NULL, # D ~ M,C
    Y_DC_model = NULL, # Y ~ D,C
    Y_DMC_model, # Y ~ D,M,C
    M_DC_model = NULL, # M ~ D,C
    data,
    d = 1,
    dstar = 0,
    minimal = FALSE,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99
){
  # ---------------Section 1: Check the arguments -----------------#
  # 1. Make sure the treatment is dummy and numeric:
  if (!is.numeric(data[[D]])) {
    stop(paste0("Variable '", D, "' must be a numeric variable."))
  }

  if (!all(unique(data[[D]][!is.na(data[[D]])]) %in% c(d, dstar))) {
    stop(
      paste0("Variable '", D, "' must be a dummy variable (d/d* only).
           Found values: ",
             paste(unique(data[[D]][!is.na(data[[D]])]), collapse = ", ")))
  }
  # 2. Make sure there is no missing values in the data:
  key_vars <- c(D, M, Y, C)
  if (!minimal) {
    miss_summary <- sapply(
      key_vars,
      FUN = function(v) c(
        nmiss = sum(!is.na(data[[v]])),
        miss = sum(is.na(data[[v]]))
      )
    ) |>
      t() |>
      as.data.frame()
    if (max(miss_summary$miss)>0) {
      warning("Warning: There is missing data in at least one of the variables
               specified for D, M, Y, and C. See the miss_summary data frame
              in the output.")
    }
  }
  # 3. Make sure the required nuisance models are properly specified:
  # Helper function to validate model formula
  check_formula <- function(model, expected_var, var_label) {
    model <- as.formula(model)
    actual_var <- as.character(attr(terms(model), "variables"))[2]
    if (actual_var != expected_var) {
      stop(
        paste0("The ", var_label, " variable in the model is ", actual_var,
               ", but you expected it to be ", expected_var, ".")
      )
    }
    return(model)
  }
  # Check  nuisance function specification:
  if (is.null(M_DC_model) && (is.null(D_MC_model) || is.null(Y_DC_model))) {
    warning(
      "Please specify models for the required nuisance function(s):\n",
      "- Specify a model for P(M|C,D) to implement the MR estimator in equation (6.17);\n",
      "- Specify models for P(D|C,M) and E(Y|D,C) to implement the MR estimator in equation (6.20);\n",
      "- Specify all the above models to implement both MR estimators."
    )
  } else {
    method_type <- c()
    if (!is.null(M_DC_model)) {
      method_type <- c(method_type, 1)
      # P(M|D,C):
      M_DC_model <- check_formula(M_DC_model, M, "outcome")
      if (!all(unique(data[[M]][!is.na(data[[M]])]) %in% c(0, 1))) {
        stop(
          paste0("Variable '", M, "' must be a dummy variable (0/1 only).
          Found values: ",
                 paste(unique(data[[M]][!is.na(data[[M]])]), collapse = ", ")))
      }
    }

    if (!is.null(D_MC_model) && !is.null(Y_DC_model)) {
      method_type <- c(method_type, 2)
      # π(D|M,C):
      D_MC_model <- check_formula(D_MC_model, D, "exposure")
      # E(Y|D,C):
      Y_DC_model <- check_formula(Y_DC_model, Y, "outcome")
      if(any(c(unlist(M)) %in% attr(terms(Y_DC_model), "term.labels"))){
        stop(
          paste0("The outcome model should only include baseline covariates
              and treatment variables; mediators are incorrectly included."
          )
        )
      }
    }
  }

  # 4. Make sure the formula specification here is right:
  # 4.1 D_C_model: E(D|C)
  D_C_model <- check_formula(D_C_model, D, "exposure")
  if(any(c(D,unlist(M)) %in% attr(terms(D_C_model), "term.labels"))){
    stop(
      paste0("The exposure model should only include baseline covariates
              as regressors; treatment and mediators are incorrectly included."
      )
    )
  }

  # 4.2 Y_DMC_model: E(Y|D,M,C)
  Y_DMC_model <- check_formula(Y_DMC_model, Y, "outcome")

  # Code the treatment variable to match the specified d and dstar:
  data <-
    data %>%
    mutate(
      !!sym(D) :=
        if_else(!!sym(D) == d, 1, 0)
    )

  # ------------Section 2: Initialize some common specifications --------------#
  # Prediction Matrix for μD(C,M):
  # For E[Y|C,M,d*|C,d*], E[Y|C,M,d*], P(M|C,d*)
  pred_dstar <- data %>% mutate(!!sym(D) := dstar)
  # For E[Y|C,M,d|C,d], E[Y|C,M,d], P(M|C,d),:
  pred_d <- data %>% mutate(!!sym(D) := d)

  # -----Fit the exposure model πD(C), πD(C,M),get the predicted weights------#
  # πD(C):P(D|C):
  pi_DC <- glm(D_C_model, family = binomial("logit"), data = data)
  data$pi_hat_DC <- pi_DC$fitted.values

  # ------Fit the outcome model: μD(C,M), get the predicted values ----------#

  mu_DMC <- lm(Y_DMC_model, data = data)
  data$mu_hat_DMC_d <- predict(mu_DMC, newdata = pred_d)
  data$mu_hat_DMC_dstar <- predict(mu_DMC, newdata = pred_dstar)

# ==============================================================================#
#   Section 3: Initialize additional specifications, and estimate separately
# ==============================================================================#

  if(any(method_type %in% list(c(2),c(1,2)))){
    #=======================================#
    # Additional Specification for mrmed2:
    #=======================================#
    # πD(C,M):P(D|C,M)
    pi_DCM <- glm(D_MC_model, family = binomial("logit"), data = data)
    data$pi_hat_DCM <- pi_DCM$fitted.values

    # ----Fit the outcome model: νD(C,M), get the predicted values --------------#
    nu_DMC_d_hat_model <-
      reformulate(
        attr(terms(Y_DC_model), "term.labels"),
        response = "mu_hat_DMC_d")
    nu_DMC_dstar_hat_model <-
      reformulate(
        attr(terms(Y_DC_model), "term.labels"),
        response = "mu_hat_DMC_dstar")
    nu_DMC_d_hat <- lm(nu_DMC_d_hat_model, data = data)
    nu_DMC_dstar_hat <- lm(nu_DMC_dstar_hat_model, data = data)
    data$nu_d_d <- predict(nu_DMC_d_hat, newdata = pred_d)
    data$nu_dstar_d <- predict(nu_DMC_d_hat, newdata = pred_dstar)
    data$nu_d_dstar <- predict(nu_DMC_dstar_hat, newdata = pred_d)
    data$nu_dstar_dstar <- predict(nu_DMC_dstar_hat, newdata = pred_dstar)

    #=======================================#
    # Generate the final results for mrmed2: #
    #=======================================#

    stat_df2 <-
      data %>%
      select(
        Y,
        D,
        .data$pi_hat_DC,
        .data$pi_hat_DCM,
        .data$mu_hat_DMC_d,
        .data$mu_hat_DMC_dstar,
        .data$nu_d_dstar,
        .data$nu_d_d,
        .data$nu_dstar_d,
        .data$nu_dstar_dstar
      ) %>%
      rename(
        pi_hat_d_DCM = .data$pi_hat_DCM,
        pi_hat_d_DC = .data$pi_hat_DC
      ) %>%
      mutate(
        pi_hat_dstar_DCM = 1 - .data$pi_hat_d_DCM,
        pi_hat_dstar_DC = 1 - .data$pi_hat_d_DC
      )

     val_map <- list("d" = d, "dstar" = dstar)

     final_calculation2 <-
      lapply(
        list(
          c("dstar", "d"),
          c("d", "d"),
          c("d", "dstar"),
          c("dstar", "dstar")
        ),
        function(comb){
          d1 <- comb[1]
          d1_val <- val_map[[d1]]
          d2 <- comb[2]
          d2_val <- val_map[[d2]]
          rst_df <-
            stat_df2 %>%
            mutate(
              # Weight in the first Line:
              !!sym(paste0("W1_", d1 ,"_", d2)) :=
                ((as.double(.data[[D]] == d2_val))/ !!sym(paste0("pi_hat_",d1,"_DC"))) *
                ((!!sym(paste0("pi_hat_",d1,"_DCM")))/ (!!sym(paste0("pi_hat_",d2,"_DCM")))),
              # Weight in the second Line:
              !!sym(paste0("W2_", d1)) :=
                ((as.double(.data[[D]] == d1_val))/ !!sym(paste0("pi_hat_",d1,"_DC")))
            )

          if(censor == TRUE){
            # censor the Weight accordingly:
            rst_df[[paste0("W1_", d1, "_", d2)]][rst_df[[D]] == d2_val] <-
              trimQ(
                rst_df[[paste0("W1_", d1, "_", d2)]][rst_df[[D]] == d2_val],
                low = censor_low,
                high = censor_high)

            rst_df[[paste0("W2_", d1)]][rst_df[[D]] == d1_val] <-
              trimQ(
                rst_df[[paste0("W2_", d1)]][rst_df[[D]] == d1_val],
                low = censor_low,
                high = censor_high)
          }

          rst_trm_df <-
            rst_df %>%
            mutate(
              !!sym(paste0("S_",d1,"_",d2)) :=
                # First Line in equation (6.17)
                !!sym(paste0("W1_", d1 ,"_", d2)) * (!!sym(Y) - !!sym(paste0("mu_hat_DMC_",d2))) +
                !!sym(paste0("W2_", d1)) * (!!sym(paste0("mu_hat_DMC_",d2)) - !!sym(paste0("nu_",d1,"_",d2))) +
                !!sym(paste0("nu_",d1,"_",d2))
            ) %>%
            mutate(
              .row_id = row_number()
            )
          return(rst_trm_df)
        }
      )

    final_df2 <-
      purrr::reduce(
        final_calculation2,
        left_join,
        by =
          purrr::reduce(
            purrr::map(
              final_calculation2,
              colnames
            ),
            intersect)
      ) %>%
      select(
        -.data$W2_d.x,
        -.data$W2_dstar.x
      ) %>%
      rename(
        W2_dstar = .data$W2_dstar.y,
        W2_d = .data$W2_d.y
      ) %>%
      mutate(
        ATE = .data$S_d_d - .data$S_dstar_dstar,
        NDE = .data$S_dstar_d - .data$S_dstar_dstar,
        NIE = .data$S_d_d - .data$S_dstar_d
      ) %>%
      summarise(
        `ATE(1,0)` = wtd.mean(.data$ATE),
        `NDE(1,0)` = wtd.mean(.data$NDE),
        `NIE(1,0)` = wtd.mean(.data$NIE)
      )

    model2_rst <-
      list(
        est2 = final_df2,
        models_D = list(
          pi_DC = pi_DC,
          pi_DCM = pi_DCM
        ),
        models_Y = list(
          mu_DMC = mu_DMC,
          nu_DMC_d_hat = nu_DMC_d_hat,
          nu_DMC_dstar_hat = nu_DMC_dstar_hat
        )
      )
  }


  if(any(method_type %in% list(c(1),c(1,2)))){

    #=======================================#
    # Additional Specification for mrmed1:
    #=======================================#
    # For E[Y|C,d,m]:
    pred_d_m <- data %>% mutate(!!sym(D) := d) %>% mutate(!!sym(M) := 1)
    # For E[Y|C,d,m*]:
    pred_d_mstar <- data %>% mutate(!!sym(D) := d) %>% mutate(!!sym(M) := 0)
    # For E[Y|C,d*,m]:
    pred_dstar_m <- data %>% mutate(!!sym(D) := dstar) %>% mutate(!!sym(M) := 1)
    # For E[Y|C,d*,m*]:
    pred_dstar_mstar <- data %>% mutate(!!sym(D) := dstar) %>% mutate(!!sym(M) := 0)
    # ---------Fit the mediator model P(M|C,D),get the predicted weights--------#
    M_DC <- glm(M_DC_model, family = binomial("logit"), data = data)
    data$M_hat_d <- predict(M_DC, newdata = pred_d, type = "response")
    data$M_hat_dstar <- predict(M_DC, newdata = pred_dstar, type = "response")
    # ------- For μd(C,m), μd(C,m*), μd*(C,m), μd*(C,m*)--------#
    data$mu_hat_DMC_d_m <- predict(mu_DMC, newdata = pred_d_m)
    data$mu_hat_DMC_d_mstar <- predict(mu_DMC, newdata = pred_d_mstar)
    data$mu_hat_DMC_dstar_m <- predict(mu_DMC, newdata = pred_dstar_m)
    data$mu_hat_DMC_dstar_mstar <- predict(mu_DMC, newdata = pred_dstar_mstar)

    #=======================================#
    # Generate the final results for mrmed1: #
    #=======================================#
    stat_df1 <-
      data %>%
      select(
        Y,
        D,
        M,
        .data$pi_hat_DC,
        starts_with("M_hat_"),
        starts_with("mu_hat_")
      ) %>%
      rename(
        pi_hat_d_DC = .data$pi_hat_DC
      ) %>%
      mutate(
        M_hat_d_m = .data$M_hat_d,
        M_hat_dstar_m = .data$M_hat_dstar,
        pi_hat_dstar_DC = 1 - .data$pi_hat_d_DC,
        M_hat_d_mstar = 1 - .data$M_hat_d_m,
        M_hat_dstar_mstar = 1 - .data$M_hat_dstar_m
      ) %>%
      mutate(
        M_hat_d = ifelse(
          !!sym(M) == 1,
          .data$M_hat_d_m,
          1 - .data$M_hat_d_m
        ),
        M_hat_dstar = ifelse(
          !!sym(M) == 1,
          .data$M_hat_dstar_m,
          1 - .data$M_hat_dstar_m
        )
      )

    val_map <- list("d" = d, "dstar" = dstar)

    final_calculation1 <-
      lapply(
        list(
          c("dstar", "d"),
          c("d", "d"),
          c("d", "dstar"),
          c("dstar", "dstar")
        ),
        function(comb){
          d1 <- comb[1]
          d1_val <- val_map[[d1]]
          d2 <- comb[2]
          d2_val <- val_map[[d2]]
          rst_df <-
            stat_df1 %>%
            mutate(
              # Weight in the first Line:
              !!sym(paste0("W1_", d1 ,"_", d2)) :=
                ((as.double(.data[[D]] == d2_val))/ !!sym(paste0("pi_hat_",d2,"_DC"))) *
                ((!!sym(paste0("M_hat_",d1)))/ (!!sym(paste0("M_hat_",d2)))),
              # Weight in the second Line:
              !!sym(paste0("W2_", d1)) :=
                ((as.double(.data[[D]] == d1_val))/ !!sym(paste0("pi_hat_",d1,"_DC")))
            )

          if(censor == TRUE){
            # Trim the Weight accordingly:
            rst_df[[paste0("W1_", d1, "_", d2)]][rst_df[[D]] == d2_val] <-
              trimQ(
                rst_df[[paste0("W1_", d1, "_", d2)]][rst_df[[D]] == d2_val],
                low = censor_low,
                high = censor_high)

            rst_df[[paste0("W2_", d1)]][rst_df[[D]] == d1_val] <-
              trimQ(
                rst_df[[paste0("W2_", d1)]][rst_df[[D]] == d1_val],
                low = censor_low,
                high = censor_high)
          }

          rst_trm_df <-
            rst_df %>%
            mutate(
              !!sym(paste0("S_",d1,"_",d2)) :=
                # First Line in equation (6.14):
                !!sym(paste0("W1_", d1 ,"_", d2)) * (!!sym(Y) - !!sym(paste0("mu_hat_DMC_",d2))) +
                # Second Line in equation (6.14):
                !!sym(paste0("W2_", d1)) * (
                  !!sym(paste0("mu_hat_DMC_",d2)) -
                    (
                      !!sym(paste0("mu_hat_DMC_",d2,"_","m")) * !!sym(paste0("M_hat_",d1,"_","m"))  +
                        !!sym(paste0("mu_hat_DMC_",d2,"_","mstar")) * !!sym(paste0("M_hat_",d1,"_","mstar"))
                    )
                ) +
                # Third Line in equation (6.14):
                (
                  !!sym(paste0("mu_hat_DMC_",d2,"_","m")) * !!sym(paste0("M_hat_",d1,"_","m"))  +
                    !!sym(paste0("mu_hat_DMC_",d2,"_","mstar")) * !!sym(paste0("M_hat_",d1,"_","mstar"))
                )
            ) %>%
            mutate(
              .row_id = row_number()
            )
          return(rst_trm_df)
        }
      )

    final_df1 <-
      purrr::reduce(
        final_calculation1,
        left_join,
        by =
          purrr::reduce(
            purrr::map(
              final_calculation1,
              colnames
            ),
            intersect)
      ) %>%
      select(
        -.data$W2_d.x,
        -.data$W2_dstar.x
      ) %>%
      rename(
        W2_dstar = .data$W2_dstar.y,
        W2_d = .data$W2_d.y
      ) %>%
      mutate(
        ATE = .data$S_d_d - .data$S_dstar_dstar,
        NDE = .data$S_dstar_d - .data$S_dstar_dstar,
        NIE = .data$S_d_d - .data$S_dstar_d
      ) %>%
      summarise(
        `ATE(1,0)` = wtd.mean(.data$ATE),
        `NDE(1,0)` = wtd.mean(.data$NDE),
        `NIE(1,0)` = wtd.mean(.data$NIE)
      )
    model1_rst <-
      list(
        est1 = final_df1,
        model_D = list(
          pi_DC = pi_DC
        ),
        model_M = list(
          M_DC = M_DC
        ),
        models_Y = list(
          mu_DMC = mu_DMC
        )
      )
  }

  # Compile and Output:
  if (minimal){
    if(all(method_type == c(2))){
      return(
        list(est2  = model2_rst$est2)
        )
    }else if(all(method_type == c(1))){
      return(
       list(est1 = model1_rst$est1)
        )
    }else(
      return(
        list(
          est1 = model1_rst$est1,
          est2 = model2_rst$est2
        )
      )
    )
  }else{
    if(all(method_type == c(2))){
      return(
        model2_rst
      )
    }else if(all(method_type == c(1))){
      return(
        model1_rst
      )
    }else(
      return(
        list(
          model1 = model1_rst,
          model2 = model2_rst
        )
      )
    )
  }
}

#' Multiply robust (MR) estimator for total, natural direct and natural
#' indirect effects
#'
#' @description
#' `mrmed()` uses a multiply robust (MR) approach to estimate
#' the total effect (ATE), natural direct effect (NDE), and natural indirect
#' effect (NIE). The function supports estimation of both single mediator and
#' multiple mediator effects.
#'
#' @details
#' `mrmed()` performs causal mediation analysis using multiply robust estimation,
#' and computes inferential statistics using the nonparametric bootstrap. It will
#' estimate the effects using different formulas depending on the nuisance function
#' specifications.
#'
#' For the Type 1 estimator (equation 6.17 in Wodtke and Zhou), `mrmed()` requires a logit model for P(M|D,C),
#' a logit model for P(D|C), and a linear model for E(Y|C,D,M).
#'
#' For the Type 2 estimator (equation 6.20 in Wodtke and Zhou), it requires a logit model for P(D|C), another logit
#' model for P(D|C,M), a linear model for E(Y|C,M,D), and another linear model for E(E(Y|C,M,D=d)|C,D).
#'
#' When analyzing multiple mediators simultaneously, users must employ the
#' Type 2 estimator, in which case `mrmed()` estimates multivariate natural
#' effects operating through the entire set of mediators.
#'
#' @param data A data frame.
#' @param D A character scalar identifying the name of the exposure variable in
#'   `data`. `D` is a character string, but the exposure variable it identifies
#'   must be numeric and binary, with two distinct values.
#' @param M A character vector (of one or more elements) identifying the names
#'   of the mediator variables in `data`. If you are estimating univariate
#'   natural effects (with a single mediator), `M` should be a character scalar
#'   (i.e., a vector with only one element)—e.g., `M = "ever_unemp_age3539"`. If you
#'   are estimating multivariate natural effects (with multiple mediators), `M`
#'   should be a list identifying all mediators—e.g.,
#'   `M = list("ever_unemp_age3539", "log_faminc_adj_age3539")`.
#' @param Y A character scalar identifying the name of the outcome variable in
#'   `data`. `Y` is a character string, but the outcome variable it identifies
#'   must be numeric.
#' @param C A character vector (of one or more elements) identifying the names
#'   of the covariate variables in `data` that you wish to include in both the
#'   mediator and outcome models. If there are no such covariates you wish to
#'   include, leave `C` as its default null argument.
#' @param D_C_model A character scalar specifying the formula to be fitted for a
#'   logit model of the exposure given baseline covariates (denoted in the book as
#'   π(D|C)). This specification is required for both types of formula. E.g.,
#'   `D_C_model = "att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3"`.
#' @param D_MC_model A character scalar specifying the formula to be fitted for a
#'   logit model of the exposure given baseline covariates and the mediator(s)
#'   (denoted in the book as π(D|C,M)). This specification is required only for
#'   the Type 2 estimator. When this input is not NULL, the function will implement the Type 2
#'   estimator by default and will throw an error if other models required for Type 2 estimation are
#'   NULL. E.g.,
#'   `D_MC_model = "att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + ever_unemp_age3539"`.
#' @param Y_DMC_model A character scalar specifying the formula to be fitted for a
#'   linear model of the outcome given baseline covariates, mediator(s), and the
#'   treatment variable (denoted in the book as μ(Y|C,M,D)). This specification is
#'   required for both types of estimator. E.g.,
#'   `Y_DMC_model = "std_cesd_age40 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + att22 + ever_unemp_age3539"`.
#' @param Y_DC_model A character scalar specifying the formula to be fitted for a
#'   linear model of the conditional mean of μ(Y|C,M,D) given baseline covariates
#'   and the treatment variable (denoted in the book as ν_D(C)). This specification
#'   allows the user to specify interactions between D and C. In implementation,
#'   the outcome variable is substituted with the estimated conditional mean from
#'   the `Y_DMC_model`. This specification is required only for the Type 2 estimator.
#'   When this input is not NULL, the function will implement the Type 2 estimator by
#'   default and will throw an error if other models required for Type 2 estimation are
#'   NULL. E.g.,
#'   `Y_DC_model = "std_cesd_age40 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + att22"`.
#' @param M_DC_model A character scalar specifying the formula to be fitted for a
#'   logit model of the conditional mean of P(M|C,D) given baseline covariates
#'   and the treatment variable. This specification allows the user to specify
#'   interactions between D and C, and is required only for Type 1 estimation.
#'    When this input is not NULL, the function will implement the Type 1 estimator by default. E.g.,
#'   `M_DC_model = "ever_unemp_age3539 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + att22"`.
#' @param d The numeric value of the treatment variable that the user defines as
#'   the treatment status. If not equal to 1, the function will recode it as 1.
#' @param dstar The numeric value of the treatment variable that the user defines
#'   as the control status. If not equal to 0, the function will recode it as 0.
#' @param censor A logical scalar indicating whether the IPW weights constructed by
#'   estimation should be censored. By default, this value is `TRUE`.
#' @param censor_low,censor_high A pair of arguments, each a numeric scalar
#'   denoting a probability in \[0,1\]. If `censor` is TRUE, then IPW weights below
#'   the `censor_low` quantile will be bottom-coded, and weights above the
#'   `censor_high` quantile will be top-coded. E.g., if the default values
#'   `censor_low = 0.01` and `censor_high = 0.99` are used, then IPW weights will
#'   be censored at their 1st and 99th percentiles. By default, weights are censored
#'   to the \[1st, 99th\] percentile range.
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
#'   `doParallel`, `doRNG`, and `foreach`. (However, you do not need to explicitly
#'   load these packages using `library()`.) Note that the results of the parallelized
#'   bootstrap may differ slightly from those of the non-parallelized bootstrap, even if
#'   the same seed is specified, due to differences in how seeds are processed.
#' @param boot_cores An integer scalar specifying the number of CPU cores to use
#'   for the parallelized bootstrap. This argument only affects computation if both
#'   `boot` and `boot_parallel` are TRUE. By default, `boot_cores` is set to the greater
#'   of two values: (a) one, and (b) the number of available CPU cores minus two.
#'   If `boot_cores` equals one, the bootstrap loop will not be parallelized,
#'   regardless of the value of `boot_parallel`.
#'
#' @returns Based on the user's specification, `mrmed()` returns a list with the
#' following elements:
#'
#' \item{est1, est2}{A tibble containing the point estimates of \eqn{ATE(1,0)}, \eqn{NDE(1,0)},
#' and \eqn{NIE(1,0)}. When the user specifies all the nuisance functions for both the
#' Type 1 and Type 2 estimators, the function will return a list of estimated results for both;
#' otherwise, it will return a tibble named `est1` for Type 1 or `est2` for Type 2.}
#'
#' \item{models_D}{A list of model objects from the fitted exposure models. For Type 1
#' estimation, this corresponds to the π(D|C) model (`D_C_model`); for Type 2,
#' this includes both the π(D|C,M) model (`D_MC_model`) and π(D|C) model (`D_C_model`).}
#'
#' \item{models_M}{A model object from the fitted mediator model. This will only be returned
#' if the user specifies a model for the nuisance function `M_DC_model` for Type 1 estimation.}
#'
#' \item{models_Y}{A list of model objects from the fitted outcome models.
#'  For the Type 1 estimator, this includes μ(Y|C,M,D), corresponding to
#'  `Y_DMC_model`. For the Type 2 estimator, it additionally includes
#'  the ν_D(C) model under D = d and D = d*.}
#'
#' If you request the bootstrap (by setting the `boot` argument to TRUE), the
#' function returns all of the elements listed above, as well as the
#' following additional elements:
#'
#' \item{ci_ATE}{A numeric vector with the bootstrap confidence interval for the
#'   total average treatment effect (ATE).}
#' \item{ci_NDE}{A numeric vector with the bootstrap confidence interval for the
#'   natural direct effect (NDE).}
#' \item{ci_NIE}{A numeric vector with the bootstrap confidence interval for the
#'   natural indirect effect (NIE).}
#' \item{pvalue_ATE}{A numeric scalar with the p-value from a two-sided test of
#'   whether the ATE differs from zero, as computed from the bootstrap.}
#' \item{pvalue_NDE}{A numeric scalar with the p-value from a two-sided test of
#'   whether the NDE differs from zero, as computed from the bootstrap.}
#' \item{pvalue_NIE}{A numeric scalar with the p-value from a two-sided test of
#'   whether the NIE differs from zero, as computed from the bootstrap.}
#' \item{sd_ATE}{A numeric scalar with the standard error of the bootstrapped
#'   ATE estimates.}
#' \item{sd_NDE}{A numeric scalar with the standard error of the bootstrapped
#'   NDE estimates.}
#' \item{sd_NIE}{A numeric scalar with the standard error of the bootstrapped
#'   NIE estimates.}
#' \item{method_type}{A vector indicating the method type returned. Can be `c(1)`,
#'   `c(2)`, or `c(1, 2)`.}
#' \item{boot_ATE}{A numeric vector of length `boot_reps` containing the ATE
#'   estimates from all bootstrap replicate samples.}
#' \item{boot_NDE}{A numeric vector of length `boot_reps` containing the NDE
#'   estimates from all bootstrap replicate samples.}
#' \item{boot_NIE}{A numeric vector of length `boot_reps` containing the NIE
#'   estimates from all bootstrap replicate samples.}
#'
#' @import rlang
#' @importFrom dplyr select mutate bind_rows if_else sample_frac summarise left_join
#' @importFrom dplyr as_tibble starts_with `%>%` case_when rename row_number n
#' @importFrom Hmisc wtd.mean
#' @export
#'
#' @examples
#' #-----------------------------------------#
#' #  Initial Specification and Clean the Data:
#' #----------------------------------------#
#' data(nlsy)
#' # outcome
#' Y <- "std_cesd_age40"
#'
#' # exposure
#' D <- "att22"
#'
#' # mediators
#' M <- list(
#' "ever_unemp_age3539",
#' "log_faminc_adj_age3539"
#' )
#' # baseline confounders:
#' C <- c(
#' "female",
#' "black",
#' "hispan",
#' "paredu",
#' "parprof",
#' "parinc_prank",
#' "famsize",
#' "afqt3"
#' )
#' # key variables
#' key_vars <- c(
#'  "cesd_age40", # unstandardized version of Y
#'   D,
#'   unlist(M),
#'   C
#' )
#'
#' # Clean the Data:
#' df <-
#' nlsy[complete.cases(nlsy[,key_vars]),] |>
#' dplyr::mutate(std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40))
#'
#' #-----------------------------------------#
#' #  Specify the models:
#' #----------------------------------------#
#' # D Models:
#' D_C_model <- as.formula(paste(D, " ~ ", paste(C, collapse= "+")))
#' D_MC_model <- as.formula(paste(D, " ~ ", paste(c(C, M[1]), collapse= "+")))
#' # Y Models:
#' Y_DC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D), collapse= "+")))
#' Y_DMC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D, M[1]), collapse= "+")))
#'
#' # Example 1: Single Mediator:
#' #-----------------------------------------#
#' # Replicate Table 6-2 Col 3:
#' #----------------------------------------#
#' \dontrun{
#' mrmed_rst1 <-
#'   mrmed(
#'     D,
#'     Y,
#'     M[[1]],
#'     C,
#'     D_C_model,
#'     D_MC_model,
#'     Y_DC_model,
#'     Y_DMC_model,
#'     M_DC_model = NULL,
#'     data = df,
#'     d = 1,
#'     dstar = 0,
#'     censor = TRUE,
#'     censor_low = 0.01,
#'     censor_high = 0.99,
#'     boot = TRUE,
#'     boot_reps = 2000,
#'     boot_conf_level = 0.95,
#'     boot_seed = 02138,
#'     boot_parallel = FALSE,
#'   )
#'  }


#---------------------- The Outer Function ----------------------#

mrmed <- function(
    D,
    Y,
    M,
    C,
    D_C_model, # D ~ C
    D_MC_model = NULL, # D ~ M,C
    Y_DC_model = NULL, # Y ~ D,C
    Y_DMC_model, # Y ~ D,M,C
    M_DC_model = NULL, # M ~ D,C
    data,
    d = 1,
    dstar = 0,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99,
    boot = TRUE,
    boot_reps = 2,
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
  }

  # compute point estimates
  est <-
    mrmed_inner(
      D,
      Y,
      M,
      C,
      D_C_model,
      D_MC_model,
      Y_DC_model,
      Y_DMC_model,
      M_DC_model,
      data = data_outer,
      d,
      dstar,
      minimal = FALSE,
      censor,
      censor_low,
      censor_high
    )

  # bootstrap, if requested
  if (boot) {
    # bootstrap function
    boot_fnc <- function() {
      # sample from the data with replacement
      boot_data <- data_outer %>% sample_frac(replace = TRUE)

      # compute point estimates in the replicate sample
      mrmed_inner(
        D,
        Y,
        M,
        C,
        D_C_model,
        D_MC_model,
        Y_DC_model,
        Y_DMC_model,
        M_DC_model,
        data = boot_data,
        d,
        dstar,
        minimal = TRUE,
        censor,
        censor_low,
        censor_high
      )
    }


    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type = "PSOCK")
      doParallel::registerDoParallel(cl = x_cluster)
      parallel::clusterExport(
        cl = x_cluster,
        varlist = c("mrmed_inner", "censor", "boot_fnc"),  # add other needed functions/vars
        envir = environment()
      )
      `%dopar%` <- foreach::`%dopar%`
    }

    # set seed
    if (!is.null(boot_seed)) {
      set.seed(boot_seed)
      if (boot_parallel_rev) {
        doRNG::registerDoRNG(boot_seed)
      }
    }

    # compute estimates for each replicate sample
    if (boot_parallel_rev) {
      boot_res <- foreach::foreach(
        i = 1:boot_reps,
        .combine = bind_rows,
        .packages = c("dplyr", "rlang", "tidyr", "purrr", "Hmisc", "tibble")
      ) %dopar% {
        out <- boot_fnc()
        out_filtered <- out[!sapply(out, is.null)]
        purrr::imap_dfr(out_filtered, ~
                          .x %>%
                          mutate(method_type = .y, boot_id = i)
        )
      }
    } else {
      boot_res <- bind_rows(
        lapply(seq_len(boot_reps), function(i) {
          out <- boot_fnc()
          out_filtered <- out[!sapply(out, is.null)]
          purrr::imap_dfr(out_filtered, ~
                            .x %>%
                            mutate(method_type = .y, boot_id = i)
          )
        })
      )
    }
    # clean up
    if (boot_parallel_rev) {
      parallel::stopCluster(x_cluster)
      rm(x_cluster)
    }
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
        list(
        ci_ATE = boot_ci(rst$`ATE(1,0)`),
        ci_NDE = boot_ci(rst$`NDE(1,0)`),
        ci_NIE = boot_ci(rst$`NIE(1,0)`),
        pvalue_ATE = boot_pval(rst$`ATE(1,0)`),
        pvalue_NDE = boot_pval(rst$`NDE(1,0)`),
        pvalue_NIE = boot_pval(rst$`NIE(1,0)`),
        sd_ATE = sd(rst$`ATE(1,0)`, na.rm = TRUE),
        sd_NDE = sd(rst$`NDE(1,0)`, na.rm = TRUE),
        sd_NIE = sd(rst$`NIE(1,0)`, na.rm = TRUE),
        method_type = unique(rst$method_type),
        boot_ATE = rst$`ATE(1,0)`,
        boot_NDE = rst$`NDE(1,0)`,
        boot_NIE = rst$`NIE(1,0)`
        )
      }
    )

    # final output
    out <- est
    if (boot) {
     out[["boot_rst_lst"]] <- boot_res_lst
    }
    return(out)
}







