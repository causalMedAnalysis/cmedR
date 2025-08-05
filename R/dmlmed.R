#' Debiased Machine Learning Estimation for Natural Effects
#'
#' @description
#' `dmlmed()` uses debiased machine learning to estimate the total effect (ATE)
#' natural direct effect (NDE), and natural indirect effect (NIE). The function
#' supports estimation of both single mediator and multiple mediator effects.
#'
#' @details
#' The estimator implemented in `dmlmed()` follows the debiased machine learning (DML) framework
#' of Chernozhukov et al. (2018), combining machine learning with repeated cross-fitting
#' to enable valid inference for causal effects. The method addresses bias that can
#' arise when the same data are used to both estimate high-dimensional nuisance functions and
#' evaluate treatment effects.
#'
#' By separating the data used to estimate nuisance functions from the data used to estimate the
#' target effect, and by repeating this process across multiple folds, DML reduces overfitting bias
#' and allows for robust inference. The resulting estimator is \eqn{\sqrt{n}}-consistent and normal,
#' even if the nuisance models are not themselves \eqn{\sqrt{n}}-consistent. This is made possible
#' by the multiplicative form of the bias: convergence is guaranteed if the product of the nuisance
#' model convergence rates exceeds \eqn{n^{-1/2}}, which can be achieved using many modern machine
#' learning methods such as LASSO, random forests, or neural networks.
#'
#' The function accommodates two types of identification strategies based on user-specified
#' nuisance functions:
#'
#' - **Type 1 Specification**: Requires an exposure model (π(D|C)), a mediator model (P(M|D,C)),
#'   and an outcome model (μ(Y|C,M,D)).
#' - **Type 2 Specification**: Requires two exposure models (π(D|C) and π(D|C,M)), and two outcome
#'   models (μ(Y|C,M,D) and its projected version ν_D(C)). This setup is particularly well-suited
#'   for analyses of multivariate mediators.
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
#'   model of the exposure given baseline covariates (denoted in the book as
#'   π(D|C)). This specification is required for Type 1 and Type 2 estimators. E.g.,
#'   `D_C_model = "att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3"`.
#' @param D_MC_model A character scalar specifying the formula to be fitted for a
#'   model of the exposure given baseline covariates and the mediator(s)
#'   (denoted in the book as π(D|C,M)). This specification is required only for
#'   Type 2 estimation. When this input is not NULL, the model will compute the Type 2
#'   estimator by default. E.g.,
#'   `D_MC_model = "att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + ever_unemp_age3539"`.
#' @param Y_DMC_model A character scalar specifying the formula to be fitted for a
#'   model of the outcome given baseline covariates, mediator(s), and the
#'   treatment variable (denoted in the book as μ(Y|C,M,D)). This specification is
#'   required for both Type 1 and Type 2 estimators. E.g.,
#'   `Y_DMC_model = "std_cesd_age40 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + att22 + ever_unemp_age3539"`.
#' @param Y_DC_model A character scalar specifying the formula to be fitted for a
#'   model of the conditional mean of μ(Y|C,M,D) given baseline covariates
#'   and the treatment variable (denoted in the book as ν_D(C)). This specification
#'   allows the user to specify interactions between D and C. In implementation,
#'   the outcome variable is substituted with the estimated conditional mean from
#'   the `Y_DMC_model`. This specification is required only for Type 2 estimation.
#'   When this input is not NULL, the model will compute the Type 2 estimator by
#'   default. E.g.,
#'   `Y_DC_model = "std_cesd_age40 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + att22"`.
#' @param M_DC_model A character scalar specifying the formula to be fitted for a
#'   model of P(M|C,D), the conditional probability of the mediator given baseline covariates
#'   and the treatment variable. This specification allows the user to specify
#'   interactions between D and C, and is required only for Type 1 estimation.
#'   When this input is not NULL, the model will compute the Type 1 estimator by default. E.g.,
#'   `M_DC_model = "ever_unemp_age3539 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + att22"`.
#' @param d The numeric value of the treatment variable that the user defines as
#'   the treatment status. If not equal to 1, the function will recode it as 1.
#' @param dstar The numeric value of the treatment variable that the user defines
#'   as the control status. If not equal to 0, the function will recode it as 0.
#' @param K Integer. The number of folds (partitions) used for repeated cross-fitting.
#'   The data is randomly divided into \code{K} approximately equal-sized subsets.
#'   Models are trained on \code{K - 1} folds and evaluated on the held-out fold.
#'   The procedure is repeated across all partitions. Typical values range from 4 to 10.
#'   The default number is \code{5}.
#' @param V Integer. The number of folds used in the Super Learner's internal cross-validation.
#'   procedure for training and evaluating candidate algorithms. Must be an explicit
#'   integer (e.g., \code{5L}) to ensure proper handling by certain functions. Smaller sample
#'   sizes may benefit from a higher \code{V} to better balance bias and variance
#'   in model evaluation. The default number is \code{5L}.
#' @param seed Seed value for reproducibility. Controls the randomization in cross-validation
#'   and other stochastic components of the procedure.
#' @param SL.library Character vector. Specifies the set of candidate algorithms to be
#'   used in the Super Learner ensemble. Each element should be the name of a valid learner
#'   (e.g., \code{"SL.mean"}, \code{"SL.glmnet"}, \code{"SL.ranger"}). Learners can
#'   be chosen from base learners available in the \pkg{SuperLearner} package or user-defined wrappers.
#'   The default setting is \code{c("SL.mean", "SL.glmnet")}.
#' @param stratifyCV Logical. If \code{TRUE}, stratified sampling is used when creating cross-validation
#'   folds within each Super Learner, preserving the distribution of the outcome variable
#'   across folds. This is especially useful for binary or imbalanced outcomes. The default setting
#'   is \code{TRUE}.
#' @param censor A logical scalar indicating whether the IPW weights used in the DML estimators
#'   should be censored. By default, this value is `TRUE`.
#' @param censor_low,censor_high A pair of arguments, each a numeric scalar
#'   denoting a probability in \[0,1\]. If `censor` is TRUE, then IPW weights below
#'   the `censor_low` quantile will be bottom-coded, and weights above the
#'   `censor_high` quantile will be top-coded. E.g., if the default values
#'   `censor_low = 0.01` and `censor_high = 0.99` are used, then IPW weights will
#'   be censored at their 1st and 99th percentiles. By default, weights are censored
#'   to the \[1st, 99th\] percentile range.
#' @param minimal A logical scalar indicating whether the function should
#'   return only a minimal set of output.
#'
#' @returns Based on the user's specification, `dmlmed()` returns a list with the
#' following elements:
#'
#' \item{est1, est2}{A tibble containing the point estimates, standard errors, and
#' 95% confidence intervals for \eqn{ATE(1,0)}, \eqn{NDE(1,0)}, and \eqn{NIE(1,0)}}.
#'
#' If `minimal` is set to `FALSE`, the function will return the following additional items,
#' along with a summary of missingness for the input data:
#'
#' \item{df1, df2}{Data frames containing the calculated \eqn{S_{d(\ast), d(\ast)}} objects,
#' using method 1 or method 2.}
#'
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @import SuperLearner
#' @importFrom caret createFolds
#' @export
#'
#' @examples
#' #------------------------------------------#
#' # Initial Specification and Clean the Data:
#' #------------------------------------------#
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
#'
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
#'
#' # key variables
#' key_vars <- c(
#'  "cesd_age40",
#'   D,
#'   unlist(M),
#'   C
#' )
#'
#' # Clean the Data:
#' df <- nlsy[complete.cases(nlsy[,key_vars]),] |>
#' dplyr::mutate(std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40))
#'
#' #-----------------------------------------#
#' #  Specify the models:
#' #----------------------------------------#
#' # D Models:
#' D_C_model <- as.formula(paste(D, " ~ ", paste(C, collapse= "+")))
#' D_MC_model <- as.formula(paste(D, " ~ ", paste(c(C, M[1]), collapse= "+")))
#' # M model:
#' M_DC_model <- as.formula(paste(M[1], " ~ ", paste(c(C, D), collapse= "+")))
#' # Y Models:
#' Y_DC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D), collapse= "+")))
#' Y_DMC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D, M[1]), collapse= "+")))
#'
#' # Example 1: Type 1 Estimation
#' \dontrun{
#'  dmlmed1_rst <-
#'  dmlmed(
#'   D,
#'   Y,
#'   M[[1]],
#'   C,
#'   D_C_model, # D ~ C
#'   D_MC_model = NULL, # D ~ M,C
#'   Y_DC_model = NULL, # Y ~ D,C
#'   Y_DMC_model, # Y ~ D,M,C
#'   M_DC_model, # M ~ D,C
#'   data = df,
#'   d = 1,
#'   dstar = 0,
#'   K = 5,
#'   seed = 1238,
#'   SL.library = c("SL.mean", "SL.glmnet","SL.ranger"),
#'   stratifyCV = TRUE,
#'   minimal = TRUE,
#'   censor = TRUE,
#'   censor_low = 0.01,
#'   censor_high = 0.99
#' )
#' }
#'
#' # Example 2: Type 2 Estimation
#' \dontrun{
#'  dmlmed2_rst <-
#'  dmlmed(
#'   D,
#'   Y,
#'   M[[1]],
#'   C,
#'   D_C_model, # D ~ C
#'   D_MC_model, # D ~ M,C
#'   Y_DC_model, # Y ~ D,C
#'   Y_DMC_model, # Y ~ D,M,C
#'   M_DC_model = NULL, # M ~ D,C
#'   data = df,
#'   d = 1,
#'   dstar = 0,
#'   K = 5,
#'   seed = 1238,
#'   SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
#'   stratifyCV = TRUE,
#'   minimal = TRUE,
#'   censor = TRUE,
#'   censor_low = 0.01,
#'   censor_high = 0.99
#' )
#' }

dmlmed <- function(
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
    d,
    dstar,
    K = 5,
    V = 5L,
    seed,
    SL.library = c("SL.mean", "SL.glmnet"),
    stratifyCV = TRUE,
    minimal = TRUE,
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

  # 2. Make sure there is no missing value in the data:
  key_vars <- c(D, unlist(M), Y, C)
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

  # 3. Make sure the nuisance functions are completely listed:
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

  if (is.null(M_DC_model) && (is.null(D_MC_model) || is.null(Y_DC_model))) {
    warning(
      "Please specify the nuisance function(s):\n",
      "- Specify P(M|C,D) to estimate equation (6.17);\n",
      "- Specify P(D|C,M) and E(Y|D,C) to estimate equation (6.20);\n",
      "- Specify both to estimate both equations."
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

  # Code the mediator and treatment variable to match the specified d and dstar:
  data <-
    data %>%
    mutate(
      !!sym(D) :=
        if_else(!!sym(D) == d, 1, 0)
    )


  # ------------ Section 2: Initialize DML estimation --------------#

  # Generate the design matrices for model fitting and prediction:
  # For model fitting:
  dm_pi_DC <- model.matrix(D_C_model, data = data)[, -1] %>% as_tibble()
  dm_mu_DMC <- model.matrix(Y_DMC_model, data = data)[, -1] %>% as_tibble()
  # For prediction:
  # μd(M,C):
  dm_mu_DMC_d <-
    model.matrix(Y_DMC_model, data = mutate(data, !!sym(D) := d))[, -1] %>% as_tibble()
  # μd*(M,C):
  dm_mu_DMC_dstar <-
    model.matrix(Y_DMC_model, data = mutate(data, !!sym(D) := dstar))[, -1] %>% as_tibble()

  # ------------ Section 3: Compute DML Estimates -----------------#

  # Section 3.1: Compute Type 2 DML Estimates (Equation 6.20):
  if(any(method_type %in% list(c(2),c(1,2)))){
    # Set the seed and generate the CV_folds:
    set.seed(seed)
    cf_folds <- caret::createFolds(data[,Y, drop = TRUE], K)
    # Generate the design matrices for model fitting and prediction:
    # For model fitting:
    # πD(M,C):P(D|C,M):
    dm_pi_DMC <- model.matrix(D_MC_model, data = data)[, -1] %>% as_tibble()
    # πD(C):P(D|C):
    dm_pi_DC <- model.matrix(D_C_model, data = data)[, -1] %>% as_tibble()
    # μD(C): E(Y|D,C):
    dm_mu_DC <- model.matrix(Y_DC_model, data = data)[, -1] %>% as_tibble()
    # For prediction:
    # μd(C):
    dm_mu_DC_d <-
      model.matrix(Y_DC_model, data = mutate(data, !!sym(D) := d))[, -1] %>% as_tibble()
    # μd*(C):
    dm_mu_DC_dstar <-
      model.matrix(Y_DC_model, data = mutate(data, !!sym(D) := dstar))[, -1] %>% as_tibble()

    pred_k_lst <- list()

    for(k in 1:K){

      #cat(" Type 2 DML Estimator: cross-fitting fold ", k, "\n")

      predict_fold <- cf_folds[[k]]
      train_k <- data[-predict_fold,]
      pred_k <- data[predict_fold,]

      # Fit Treatment Models:
      # Estimate πD(C):
      pi_DC <- SuperLearner(
        Y          = train_k[[D]],
        X          = dm_pi_DC[-predict_fold,],
        newX       = dm_pi_DC,
        family     = binomial(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
        cvControl  = list(V = V, stratifyCV = stratifyCV, shuffle = TRUE, validRows = NULL)
      )
      # Estimate πD(M,C):
      pi_DMC <- SuperLearner(
        Y          = train_k[[D]],
        X          = dm_pi_DMC[-predict_fold,],
        newX       = dm_pi_DMC,
        family     = binomial(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
        cvControl  = list(V = V, stratifyCV = stratifyCV, shuffle = TRUE, validRows = NULL)
      )
      # Fit Outcome Models:
      # Estimate μD(M,C):
      mu_DMC <- SuperLearner(
        Y          = train_k[[Y]],
        X          = dm_mu_DMC[-predict_fold,],
        family     = gaussian(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE),
        cvControl  = list(V = V, shuffle = TRUE, validRows = NULL)
      )

      # Calculate E(Y|d,M,C) and E(Y|d*,M,C) as the outcome for νD(C):
      # Estimate μd(M,C):
      data$mu_hat_DMC_d <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_d)$pred
      # Estimate μd*(M,C):
      data$mu_hat_DMC_dstar <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_dstar)$pred
      # Estimate νd*(C):
      nu_DMC_dstar_hat <- SuperLearner(
        Y          = data$mu_hat_DMC_dstar[-predict_fold],
        X          = dm_mu_DC[-predict_fold,],
        family     = gaussian(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE),
        cvControl  = list(V = V, shuffle = TRUE, validRows = NULL)
      )
      # Estimate νd(C):
      nu_DMC_d_hat <- SuperLearner(
        Y          = data$mu_hat_DMC_d[-predict_fold],
        X          = dm_mu_DC[-predict_fold,],
        family     = gaussian(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE),
        cvControl  = list(V = V, shuffle = TRUE, validRows = NULL)
      )
      # Compute predictions:
      data$pi_hat_DC <- pi_DC$SL.predict
      data$pi_hat_DCM <- pi_DMC$SL.predict
      data$nu_d_d <-
        predict.SuperLearner(nu_DMC_d_hat, newdata = dm_mu_DC_d)$pred
      data$nu_dstar_d <-
        predict.SuperLearner(nu_DMC_d_hat, newdata = dm_mu_DC_dstar)$pred
      data$nu_d_dstar <-
        predict.SuperLearner(nu_DMC_dstar_hat, newdata = dm_mu_DC_d)$pred
      data$nu_dstar_dstar <-
        predict.SuperLearner(nu_DMC_dstar_hat, newdata = dm_mu_DC_dstar)$pred
      pred_k <- data[predict_fold,]
      pred_k_lst <- append(pred_k_lst, list(pred_k))
    }

    # Generate the final results
    main_df <- bind_rows(pred_k_lst)
    stat_df2 <-
      main_df %>%
      dplyr::select(
        dplyr::all_of(Y),
        dplyr::all_of(D),
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
              !!sym(paste0("W1_", d1 ,"_", d2)) :=
                ((as.double(.data[[D]] == d2_val))/ !!sym(paste0("pi_hat_",d1,"_DC"))) *
                ((!!sym(paste0("pi_hat_",d1,"_DCM")))/ (!!sym(paste0("pi_hat_",d2,"_DCM")))),
              !!sym(paste0("W2_", d1)) :=
                ((as.double(.data[[D]] == d1_val))/ !!sym(paste0("pi_hat_",d1,"_DC")))
            )

          if(censor == TRUE){
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
                !!sym(paste0("W1_", d1 ,"_", d2)) * (!!sym(Y) - !!sym(paste0("mu_hat_DMC_",d2))) +
                !!sym(paste0("W2_", d1)) * (!!sym(paste0("mu_hat_DMC_",d2)) - !!sym(paste0("nu_",d1,"_",d2))) +
                !!sym(paste0("nu_",d1,"_",d2))
            ) %>%
            mutate(
              .row_id = dplyr::row_number()
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
      dplyr::select(
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
      )

    final_stat2 <-
      final_df2 %>%
      summarise(
        `ATE(1,0)_Mean` = mean(.data$ATE),
        `NDE(1,0)_Mean` = mean(.data$NDE),
        `NIE(1,0)_Mean` = mean(.data$NIE),
        `ATE(1,0)_SE` = sd(.data$ATE) / sqrt(n()),
        `NDE(1,0)_SE` = sd(.data$NDE) / sqrt(n()),
        `NIE(1,0)_SE` = sd(.data$NIE) / sqrt(n())
      ) %>%
      pivot_longer(
        cols = dplyr::everything(),
        names_to = c("Estimand", "stat"),
        names_sep = "_",
        values_to = "value"
      ) %>%
      pivot_wider(names_from = .data$stat, values_from = .data$value) %>%
      mutate(
        lower = round(.data$Mean - 1.96 * .data$SE, 3),
        upper = round(.data$Mean + 1.96 * .data$SE, 3),
        out = paste0(
          round(.data$Mean, 3)," [",
          .data$lower,", ",
          .data$upper, "]")
      ) %>%
      dplyr::select(
        .data$Estimand,
        .data$Mean,
        .data$SE,
        .data$out
      )
  }

  # Section 3.2: Compute Type 1 DML Estimates (Equation 6.17):

  if(any(method_type %in% list(c(1),c(1,2)))){
    # Set the seed and generate the CV_folds:
    set.seed(seed)
    cf_folds <- caret::createFolds(data[,Y,drop = TRUE], K)
    # Generate the design matrices for model fitting and prediction:
    # For Model Fitting:
    # E(M|D,C):
    dm_M <- model.matrix(M_DC_model, data = data)[, -1] %>% as_tibble()
    # For Model Prediction:
    # μd(C,m):
    dm_mu_DMC_d_m <- dm_mu_DMC_d %>% mutate(!!sym(M) := 1)
    # μd(C,m*):
    dm_mu_DMC_d_mstar <- dm_mu_DMC_d %>% mutate(!!sym(M) := 0)
    # μd*(C,m):
    dm_mu_DMC_dstar_m <- dm_mu_DMC_dstar %>% mutate(!!sym(M) := 1)
    # μd*(C,m*):
    dm_mu_DMC_dstar_mstar <- dm_mu_DMC_dstar %>% mutate(!!sym(M) := 0)
    # μd(C,M):
    dm_M_d <-
      model.matrix(M_DC_model, data = mutate(data, !!sym(D) := d))[, -1] %>% as_tibble()
    # μd*(C,M):
    dm_M_dstar <-
      model.matrix(M_DC_model, data = mutate(data, !!sym(D) := dstar))[, -1] %>% as_tibble()

    pred_k_lst <- list()

    for(k in 1:K){

      #cat(" Type 1 DML Estimator: cross-fitting fold ", k, "\n")

      predict_fold <- cf_folds[[k]]
      train_k <- data[-predict_fold,]
      pred_k <- data[predict_fold,]

      # Fit Treatment Models:
      # πD(C):
      pi_DC <- SuperLearner(
        Y          = train_k[[D]],
        X          = dm_pi_DC[-predict_fold,],
        newX       = dm_pi_DC,
        family     = binomial(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
        cvControl  = list(V = V, stratifyCV = stratifyCV, shuffle = TRUE, validRows = NULL)
      )
      # Fit Mediator Models:
      # E(M|D,C):
      M_DC <- SuperLearner(
        Y          = train_k[[M]],
        X          = dm_M[-predict_fold,],
        newX       = dm_M,
        family     = binomial(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
        cvControl  = list(V = V, stratifyCV = stratifyCV, shuffle = TRUE, validRows = NULL)
      )

      data$M_hat_d <-
        predict.SuperLearner(M_DC, newdata = dm_M_d)$pred
      data$M_hat_dstar <-
        predict.SuperLearner(M_DC, newdata = dm_M_dstar)$pred

      # Fit Outcome Models:
      # μD(M,C):
      mu_DMC <- SuperLearner(
        Y          = train_k[[Y]],
        X          = dm_mu_DMC[-predict_fold,],
        family     = gaussian(),
        SL.library = SL.library,
        control    = list(saveFitLibrary = TRUE),
        cvControl  = list(V = V, shuffle = TRUE, validRows = NULL)
      )

      # Compute predictions:
      data$pi_hat_DC <- pi_DC$SL.predict
      data$mu_hat_DMC_dstar <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_dstar)$pred
      data$mu_hat_DMC_d <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_d)$pred
      data$mu_hat_DMC_dstar_mstar <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_dstar_mstar)$pred
      data$mu_hat_DMC_d_mstar <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_d_mstar)$pred
      data$mu_hat_DMC_dstar_m <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_dstar_m)$pred
      data$mu_hat_DMC_d_m <-
        predict.SuperLearner(mu_DMC, newdata = dm_mu_DMC_d_m)$pred
      pred_k <- data[predict_fold,]
      pred_k_lst <- append(pred_k_lst, list(pred_k))
    }

    main_df <- bind_rows(pred_k_lst)

    # Generate the final results
    stat_df1 <-
      main_df %>%
      dplyr::select(
        dplyr::all_of(Y),
        dplyr::all_of(D),
        dplyr::all_of(M),
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
              !!sym(paste0("W1_", d1 ,"_", d2)) :=
                ((as.double(.data[[D]] == d2_val))/ !!sym(paste0("pi_hat_",d2,"_DC"))) *
                ((!!sym(paste0("M_hat_",d1)))/ (!!sym(paste0("M_hat_",d2)))),
              !!sym(paste0("W2_", d1)) :=
                ((as.double(.data[[D]] == d1_val))/ !!sym(paste0("pi_hat_",d1,"_DC")))
            )

          if(censor == TRUE){
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
                !!sym(paste0("W1_", d1 ,"_", d2)) * (!!sym(Y) - !!sym(paste0("mu_hat_DMC_",d2))) +
                !!sym(paste0("W2_", d1)) * (
                  !!sym(paste0("mu_hat_DMC_",d2)) -
                    (
                      !!sym(paste0("mu_hat_DMC_",d2,"_","m")) * !!sym(paste0("M_hat_",d1,"_","m"))  +
                        !!sym(paste0("mu_hat_DMC_",d2,"_","mstar")) * !!sym(paste0("M_hat_",d1,"_","mstar"))
                    )
                ) +
                (
                  !!sym(paste0("mu_hat_DMC_",d2,"_","m")) * !!sym(paste0("M_hat_",d1,"_","m"))  +
                    !!sym(paste0("mu_hat_DMC_",d2,"_","mstar")) * !!sym(paste0("M_hat_",d1,"_","mstar"))
                )
            ) %>%
            mutate(
              .row_id = dplyr::row_number()
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
      dplyr::select(
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
      )

    final_stat1 <-
      final_df1 %>%
      summarise(
        `ATE(1,0)_Mean` = mean(.data$ATE),
        `NDE(1,0)_Mean` = mean(.data$NDE),
        `NIE(1,0)_Mean` = mean(.data$NIE),
        `ATE(1,0)_SE` = sd(.data$ATE) / sqrt(n()),
        `NDE(1,0)_SE` = sd(.data$NDE) / sqrt(n()),
        `NIE(1,0)_SE` = sd(.data$NIE) / sqrt(n())
      ) %>%
      pivot_longer(
        cols = dplyr::everything(),
        names_to = c("Estimand", "stat"),
        names_sep = "_",
        values_to = "value"
      ) %>%
      pivot_wider(names_from = .data$stat, values_from = .data$value) %>%
      mutate(
        lower = round(.data$Mean - 1.96 * .data$SE, 3),
        upper = round(.data$Mean + 1.96 * .data$SE, 3),
        out = paste0(
          round(.data$Mean, 3)," [",
          .data$lower,", ",
          .data$upper, "]")
      ) %>%
      dplyr::select(
        .data$Estimand,
        .data$Mean,
        .data$SE,
        .data$out
      )
  }

  if (all(method_type == 2)) {
    return(if (minimal) final_stat2 else list(est2 = final_stat2, df2 = final_df2, miss_summary = miss_summary))
  }

  if (all(method_type == 1)) {
    return(if (minimal) final_stat1 else list(est1 = final_stat1, df1 = final_df1, miss_summary = miss_summary))
  }

  if (minimal) {
    return(list(est2 = final_stat2, est1 = final_stat1))
  } else {
    return(list(
      est = list(est2 = final_stat2, est1 = final_stat1),
      df = list(df2 = final_df2, df1 = final_df1),
      miss_summary = miss_summary
    ))
  }
}
