#' Debiased Machine Learning (DML) Estimation for Path-specific Effects
#'
#' @description
#' `dmlpath()` uses a debiased machine learning approach to estimate path-specific effects.
#' If there are K causally ordered mediators, dmlpath provides estimates for: a direct effect
#' of the exposure on the outcome that does not operate through any of the mediators, and then
#' K path-specific effects, with each of these effects operating through one mediator, net of the
#' mediators preceding it in causal order. If only one mediator is specified, `dmlpath()` computes
#' conventional natural direct and indirect effects.
#'
#' @details
#' `dmlpath()` estimates path specific effects using debiased machine learning. It will
#' repeatedly implement the type 2 multiply robust estimator in the `dmlmed()` function
#' to compute natural effects, and then take differences between these effects
#' to generate the PSEs of interest.
#'
#' The Type 2 estimator (Equation 6.20 in Wodtke and Zhou) requires modeling
#' P(D|C), P(D|C,M), E(Y|C,M,D), and E(E(Y|C,M,D=d)|C,D). The `dmlpath()` function uses
#' the SuperLearner algorithm, allowing users to specify multiple machine learning algorithms
#' and to estimate the above models non-parametrically.

#' #' To compute path-specific effects with K causally ordered mediators, `dmlpath()`
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
#' recommend that you dummy code the levels of the factor and treat the dummy
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
#' @param interaction_DM A logical scalar indicating whether the outcome model
#'   should include exposure-mediator interactions (interactions of the exposure
#'   with each mediator if there is more than one mediator in `M`).
#' @param interaction_DC A logical scalar indicating whether both the outcome
#'   model and each of the mediator models should include interactions of the
#'   exposure with each covariate in `C`.
#' @param interaction_MC A logical scalar indicating whether the outcome model
#'   should include interactions of each mediator in `M` with each covariate in
#'   `C`.
#' @param d The numeric value of the treatment variable that the user defines as
#'   the treatment status. If not equal to 1, the function will recode it as 1.
#' @param dstar The numeric value of the treatment variable that the user defines
#'   as the control status. If not equal to 0, the function will recode it as 0.
#' @param censor A logical scalar indicating whether the IPW weights that enter
#'  the estimating equations should be censored. By default, this value is `TRUE`.
#' @param censor_low,censor_high A pair of arguments, each a numeric scalar
#'   denoting a probability in \[0,1\]. If `censor` is TRUE, then IPW weights below
#'   the `censor_low` quantile will be bottom-coded, and weights above the
#'   `censor_high` quantile will be top-coded. E.g., if the default values
#'   `censor_low = 0.01` and `censor_high = 0.99` are used, then IPW weights will
#'   be censored at their 1st and 99th percentiles. By default, weights are censored
#'   to the \[1st, 99th\] percentile range.
#' @param num_folds Integer. The number of folds (partitions) used for repeated cross-fitting.
#'   The data is randomly divided into \code{K} approximately equal-sized subsets.
#'   Models are trained on \code{K - 1} folds and evaluated on the held-out fold.
#'   The procedure is repeated across all partitions. Typical values range from 4 to 10.
#'   The default number is \code{5}.
#' @param V Integer. The number of folds used in the Super Learner's internal cross-validation.
#'   procedure for training and evaluating candidate algorithms. Must be an explicit
#'   integer (e.g., \code{5L}) to ensure proper handling by certain functions.
#'   The default number is \code{5L}.
#' @param seed Seed value for reproducibility. Controls the randomization in cross-validation
#'   and other stochastic components of the estimation procedure.
#' @param SL.library Character vector. Specifies the set of candidate algorithms to be
#'   used in the Super Learner ensemble. Each element should be the name of a valid learner
#'   (e.g., \code{"SL.mean"}, \code{"SL.glmnet"}, \code{"SL.ranger"}). Learners can
#'   be chosen from the base learners available in the \pkg{SuperLearner} package or user-defined wrappers.
#'   The default learners are \code{c("SL.mean", "SL.glmnet")}.
#' @param stratifyCV Logical. If \code{TRUE}, stratified sampling is used when creating cross-validation
#' folds within each SuperLearner for the treatment, preserving the distribution of the treatment variable
#' across folds. The default setting is \code{TRUE}.

#' @returns By default, `dmlpath()` returns a tibble with the following columns:
#'
#' \item{Estimand}{A character vector describing each estimand, including the total effect (ATE)
#' and the path-specific effects (PSEs) for each mediator path. Notation follows the format used
#' in Wodtke and Zhou, e.g., \eqn{ATE(1,0)} or \eqn{PSE: D → M1 → Y(1,0)}.}
#'
#' \item{Mean}{A numeric vector of point estimates for each estimand.}
#'
#' \item{SE}{A numeric vector of standard errors corresponding to each point estimate.}
#'
#' \item{out}{A character vector combining the rounded point estimate and its 95% confidence
#' interval in the format: `"estimate [lower, upper]"`, where the interval is constructed
#' as `Mean ± 1.96 × SE`.}
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
#'
#' data(nlsy)
#'
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
#' # baseline confounders
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
#' df <-
#' nlsy[complete.cases(nlsy[,key_vars]),] |>
#' dplyr::mutate(std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40))
#' \dontrun{
#' table6_4_col2 <-
#'   dmlpath(
#'     D = D,
#'     Y = Y,
#'     M = M,
#'     C = C,
#'     data = df,
#'     d = 1,
#'     dstar = 0,
#'     censor = TRUE,
#'     censor_low = 0.01,
#'     censor_high = 0.99,
#'     interaction_DM = FALSE,
#'     interaction_DC = FALSE,
#'     interaction_MC = FALSE,
#'     num_folds = 5,
#'     V = 5L,
#'     seed = 02138,
#'     SL.library = c("SL.mean","SL.glmnet","SL.ranger"),
#'     stratifyCV = TRUE
#'   )
#'   }

dmlpath <- function(
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
    num_folds = 5,
    V = 5L,
    seed,
    SL.library = c("SL.mean", "SL.glmnet"),
    stratifyCV = TRUE
){

  # Step 1: Get dimensions
  K <- length(M)
  PSE <- vector("list", length = K + 1)
  SE <- vector("list", length = K + 1)

  # Step 2: Calculate the NDE and NIE for each k
  for(k in rev(seq_len(K))){

    # Step 2.1: Construct the estimation models
    # D Models
    D_C_model <- as.formula(paste(D, " ~ ", paste(C, collapse= "+")))

    if(interaction_MC == FALSE){
      D_MC_model <- as.formula(paste(D, " ~ ", paste(c(C, unlist(M[1:k])), collapse= "+")))
    } else {
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
    if(interaction_DC == FALSE){
      Y_DC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D), collapse= "+")))
    } else {
      Y_DC_model <- as.formula(
        paste(
          Y,
          " ~ ",
          paste(c(C, D), collapse= "+"),
          "+",
          paste(D, C, sep = ":", collapse = " + "))
      )
    }

    if(!any(interaction_DC, interaction_DM, interaction_MC)){
      Y_DMC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D, unlist(M[1:k])), collapse= "+")))
    } else {
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
      dmlmed(
        D = D,
        Y = Y,
        M = unlist(M[1:k]),
        C = C,
        D_C_model = D_C_model,
        D_MC_model = D_MC_model,
        Y_DC_model = Y_DC_model,
        Y_DMC_model = Y_DMC_model,
        M_DC_model = NULL,
        data = data,
        d = d,
        dstar = dstar,
        K = num_folds,
        V = V,
        seed = seed,
        SL.library = SL.library,
        stratifyCV = stratifyCV,
        minimal = FALSE,
        censor = censor,
        censor_low = censor_low,
        censor_high = censor_high
      )

    # Step 2.3: Calculate PSEs
    if(K == 1) {
      PSE[[K]] <- est[["est2"]] %>% dplyr::filter(.data$Estimand == "NDE(1,0)") %>% dplyr::pull(.data$Mean)
      PSE[[K + 1]] <- est[["est2"]] %>% dplyr::filter(.data$Estimand == "NIE(1,0)") %>% dplyr::pull(.data$Mean)
      SE[[K]] <- est[["est2"]] %>% dplyr::filter(.data$Estimand == "NDE(1,0)") %>% dplyr::pull(SE)
      SE[[K + 1]] <- est[["est2"]] %>% dplyr::filter(.data$Estimand == "NIE(1,0)") %>% dplyr::pull(SE)
      ATE <- est[["est2"]] %>% dplyr::filter(.data$Estimand == "ATE(1,0)")
      names(PSE) <- c("D->Y","D->M1->Y")
      names(SE) <- c("D->Y","D->M1->Y")
    }

    else if(k == 1 & K >= 2) {
      PSE[[K - k + 1]] <- mean(est$df2$NDE) - mean(prev_NDE)
      PSE[[K + 1]] <- mean(est$df2$NIE)
      SE[[K - k + 1]] <- sd(est$df2$NDE - prev_NDE) / sqrt(sum(!is.na(est$df2$NDE - prev_NDE)))
      SE[[K + 1]] <- sd(est$df2$NIE) / sqrt(sum(!is.na(est$df2$NIE)))
      ATE <- est[["est2"]] %>% dplyr::filter(.data$Estimand == "ATE(1,0)")
      names(PSE)[K + 1] <- "D->M1~>Y"
      names(SE)[K + 1] <- "D->M1~>Y"
      if(k + 1 == K){
        names(PSE)[K - k + 1] <- paste0("D->M",k+1,"->Y")
        names(SE)[K - k + 1] <- paste0("D->M",k+1,"->Y")
      } else{
        names(PSE)[K - k + 1] <- paste0("D->M",k+1,"~>Y")
        names(SE)[K - k + 1] <- paste0("D->M",k+1,"~>Y")
      }
    }
    else if(k == K & K >= 2) {
      PSE[[K - k + 1]] <- mean(est$df2$NDE)
      SE[[K - k + 1]] <- sd(est$df2$NDE) / sqrt(sum(!is.na(est$df2$NDE)))
      names(PSE)[K - k + 1] <- "D->Y"
      names(SE)[K - k + 1] <- "D->Y"
      prev_NDE <- est$df2$NDE
    }
    else {
      PSE[[K - k + 1]] <- mean(est$df2$NDE) - mean(prev_NDE)
      SE[[K - k + 1]] <- sd(est$df2$NDE - prev_NDE)/ sqrt(sum(!is.na(est$df2$NDE - prev_NDE)))
      prev_NDE <- est$df2$NDE
      if(k + 1 == K){
        names(PSE)[[K - k + 1]] <- paste0("D->M",k+1,"->Y")
        names(SE)[[K - k + 1]] <- paste0("D->M",k+1,"->Y")
      } else {
        names(PSE)[[K - k + 1]] <- paste0("D->M",k+1,"~>Y")
        names(SE)[[K - k + 1]] <- paste0("D->M",k+1,"~>Y")
      }
    }
  }

  # Collect the results
  final_df <-
    rbind(
      ATE = ATE,
      tibble::tibble(
        Estimand = paste0("PSE",":", names(PSE),"(1,0)"),
        Mean = unlist(PSE)
       ) %>%
      left_join(
       tibble::tibble(
         Estimand = paste0("PSE",":", names(SE),"(1,0)"),
         SE = unlist(SE)
       ),
       by = "Estimand"
      ) %>%
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
    )
  return(final_df)
}


