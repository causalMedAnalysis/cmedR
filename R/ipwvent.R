# -- Utility helpers ------------------------------------------------------------

# Gumbel CDF functions used by polr() for loglog/cloglog links. Defined locally
# to avoid a hidden dependency on unexported MASS internals.
pgumbel  <- function(q, loc = 0, scale = 1, lower.tail = TRUE) {
  q <- (q - loc) / scale
  p <- exp(-exp(-q))
  if (!lower.tail) 1 - p else p
}
pGumbel <- function(q, ...) 1 - pgumbel(-q, ...)

# Suppress R CMD check NOTE for the .ipwvent_polr_weights column name used via
# NSE in MASS::polr(weights = .ipwvent_polr_weights, ...).
utils::globalVariables(".ipwvent_polr_weights")

# trimQ: winsorise a numeric vector at the low and high quantile boundaries.
trimQ <- function(x, low = 0.01, high = 0.99) {
  lo <- quantile(x, probs = low,  na.rm = TRUE)
  hi <- quantile(x, probs = high, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

# comb_list_vec: foreach .combine helper -- accumulates scalar-list results into
# a list of vectors.
comb_list_vec <- function(a, b) {
  list(
    ATE = c(a$ATE, b$ATE),
    IDE = c(a$IDE, b$IDE),
    IIE = c(a$IIE, b$IIE)
  )
}
# -------------------------------------------------------------------------------

#' Interventional-effects inverse probability weighting (IPW) estimator: inner
#' function
#'
#' @description
#' Internal function used within `ipwvent()`. See the `ipwvent()` function
#' documentation for a description of shared function arguments. Here, we will
#' only document the one argument that is not shared by `ipwvent_inner()` and
#' `ipwvent()`: the `minimal` argument.
#'
#' @param minimal A logical scalar indicating whether the function should
#'   return only a minimal set of output. The `ipwvent()` function uses the
#'   default of FALSE when calling `ipwvent_inner()` to generate the point
#'   estimates and sets the argument to TRUE when calling `ipwvent_inner()`
#'   to perform the bootstrap.
#'
#' @noRd
.ipwvent_build_rhs_formula <- function(formula, response_name) {
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  as.formula(paste(response_name, "~", rhs))
}


.ipwvent_build_polr_start <- function(formula, data) {
  mf <- stats::model.frame(formula = formula, data = data, na.action = stats::na.fail)
  y <- stats::model.response(mf)

  if (!is.ordered(y)) {
    stop(paste(strwrap("Error: The ordered-logit model for the exposure-induced confounder requires an ordered response."), collapse = "\n"))
  }
  if (nlevels(y) < 3L) {
    stop(paste(strwrap("Error: The ordered-logit model for the exposure-induced confounder requires at least three observed categories in the estimation sample."), collapse = "\n"))
  }

  x_mat <- stats::model.matrix(object = attr(mf, "terms"), data = mf)
  intercept_pos <- match("(Intercept)", colnames(x_mat), nomatch = 0L)
  if (intercept_pos > 0L) {
    x_mat <- x_mat[, -intercept_pos, drop = FALSE]
  }

  cum_probs <- cumsum(prop.table(table(y)))
  thresh_probs <- cum_probs[-length(cum_probs)]
  thresh_probs <- pmin(pmax(thresh_probs, 1e-6), 1 - 1e-6)
  zeta_start <- stats::qlogis(thresh_probs)

  c(rep(0, ncol(x_mat)), zeta_start)
}


.ipwvent_fit_polr <- function(formula, data, weight_name) {
  data_polr <- data
  data_polr[[".ipwvent_polr_weights"]] <- data_polr[[weight_name]]

  tryCatch(
    MASS::polr(
      formula = formula,
      data = data_polr,
      weights = .ipwvent_polr_weights,
      Hess = TRUE,
      model = TRUE
    ),
    error = function(e) {
      if (!grepl("attempt to find suitable starting values failed", conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }

      MASS::polr(
        formula = formula,
        data = data_polr,
        weights = .ipwvent_polr_weights,
        start = .ipwvent_build_polr_start(
          formula = formula,
          data = data_polr
        ),
        Hess = TRUE,
        model = TRUE
      )
    }
  )
}


# Refit a pre-fitted glm or polr model on new data (used for bootstrap resamples
# and for the point-estimate fit on the working dataset inside ipwvent_inner).
# response_var: if non-NULL, substitute this column name as the response (needed
# for polr, which requires an ordered-factor response stored under a temp name).
.ipwvent_refit <- function(model, new_data, response_var = NULL, base_weights_rsc = NULL) {
  f <- formula(model)
  if (!is.null(response_var)) {
    f <- .ipwvent_build_rhs_formula(f, response_var)
  }
  environment(f) <- environment()

  if (inherits(model, "polr")) {
    new_data[[".base_weight_col"]] <- if (is.null(base_weights_rsc)) 1 else base_weights_rsc
    return(.ipwvent_fit_polr(f, new_data, weight_name = ".base_weight_col"))
  }

  # glm branch
  fam <- model$family
  if (is.null(base_weights_rsc)) {
    return(glm(f, data = new_data, family = fam))
  }
  glm(f, data = new_data, family = fam, weights = base_weights_rsc)
}


.ipwvent_get_l_values <- function(x) {
  if (is.factor(x) || is.ordered(x)) {
    levels(droplevels(x))
  }
  else {
    sort(unique(x))
  }
}


.ipwvent_prepare_l_response <- function(x, L_model) {
  if (inherits(L_model, "polr")) {
    if (is.ordered(x)) {
      return(droplevels(x))
    }
    if (is.factor(x)) {
      x_drop <- droplevels(x)
      return(ordered(as.character(x_drop), levels = levels(x_drop)))
    }
    return(ordered(x, levels = sort(unique(x))))
  }

  x
}


.ipwvent_set_l_value <- function(x_template, value, n) {
  if (is.ordered(x_template)) {
    return(ordered(rep(as.character(value), n), levels = levels(x_template)))
  }
  if (is.factor(x_template)) {
    return(factor(rep(as.character(value), n), levels = levels(x_template)))
  }
  if (is.character(x_template)) {
    return(rep(as.character(value), n))
  }
  rep(value, n)
}



.ipwvent_as_prob_matrix <- function(probs, n_rows) {
  if (is.null(dim(probs))) {
    probs_vec <- probs
    probs <- matrix(probs_vec, nrow = n_rows, byrow = TRUE)
    if (!is.null(names(probs_vec))) {
      colnames(probs) <- names(probs_vec)
    }
  }
  else {
    probs <- as.matrix(probs)
  }
  probs
}


.ipwvent_predict_l_probs <- function(
  model,
  newdata,
  l_values
) {
  l_level_names <- as.character(l_values)

  if (!inherits(model, "polr")) {
    p1 <- predict(model, newdata = newdata, type = "response")
    out <- cbind(
      `0` = 1 - p1,
      `1` = p1
    )
    return(out[, l_level_names, drop = FALSE])
  }

  if (inherits(model, "polr")) {
    x_terms <- stats::delete.response(stats::terms(model))
    mf <- stats::model.frame(
      formula = x_terms,
      data = newdata,
      na.action = function(x) x,
      xlev = model$xlevels
    )
    x_mat <- stats::model.matrix(
      object = x_terms,
      data = mf,
      contrasts.arg = model$contrasts
    )
    intercept_pos <- match("(Intercept)", colnames(x_mat), nomatch = 0L)
    if (intercept_pos > 0L) {
      x_mat <- x_mat[, -intercept_pos, drop = FALSE]
    }

    beta_hat <- stats::coef(model)
    if (length(beta_hat) == 0L) {
      eta <- rep(0, nrow(x_mat))
    }
    else {
      if (!all(names(beta_hat) %in% colnames(x_mat))) {
        stop(paste(strwrap("Error: The ordered-logit prediction step could not align the new-data design matrix with the coefficients retained in the fitted model for L."), collapse = "\n"))
      }
      x_mat <- x_mat[, names(beta_hat), drop = FALSE]
      eta <- drop(x_mat %*% beta_hat)
    }

    q <- length(model$zeta)
    zeta_mat <- matrix(model$zeta, nrow = nrow(x_mat), ncol = q, byrow = TRUE)
    pfun <- switch(
      model$method,
      logistic = stats::plogis,
      probit = stats::pnorm,
      loglog = pgumbel,
      cloglog = pGumbel,
      cauchit = stats::pcauchy
    )
    cumpr <- matrix(
      pfun(zeta_mat - eta),
      nrow = nrow(x_mat),
      ncol = q
    )
    probs <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
    colnames(probs) <- model$lev

    if (!all(l_level_names %in% colnames(probs))) {
      stop(paste(strwrap("Error: The fitted ordered-logit model for the exposure-induced confounder did not return probabilities for all observed levels of L."), collapse = "\n"))
    }
    return(probs[, l_level_names, drop = FALSE])
  }
}


# BUG 0.16: helper to compute P(M=k | newdata) for each row from an ordinal M
# polr model.  Mirrors the ologit branch of .ipwvent_predict_l_probs.
.ipwvent_predict_m_polr_probs <- function(model, newdata, m_values) {
  m_level_names <- as.character(m_values)

  x_terms <- stats::delete.response(stats::terms(model))
  mf <- stats::model.frame(
    formula = x_terms,
    data = newdata,
    na.action = function(x) x,
    xlev = model$xlevels
  )
  x_mat <- stats::model.matrix(
    object = x_terms,
    data = mf,
    contrasts.arg = model$contrasts
  )
  intercept_pos <- match("(Intercept)", colnames(x_mat), nomatch = 0L)
  if (intercept_pos > 0L) {
    x_mat <- x_mat[, -intercept_pos, drop = FALSE]
  }

  beta_hat <- stats::coef(model)
  if (length(beta_hat) == 0L) {
    eta <- rep(0, nrow(x_mat))
  } else {
    if (!all(names(beta_hat) %in% colnames(x_mat))) {
      stop(paste(strwrap("Error: The prediction step for the mediator model could not align the new-data design matrix with the coefficients in the fitted model for M."), collapse = "\n"))
    }
    x_mat <- x_mat[, names(beta_hat), drop = FALSE]
    eta <- drop(x_mat %*% beta_hat)
  }

  q <- length(model$zeta)
  zeta_mat <- matrix(model$zeta, nrow = nrow(x_mat), ncol = q, byrow = TRUE)
  pfun <- switch(
    model$method,
    logistic = stats::plogis,
    probit   = stats::pnorm,
    loglog   = pgumbel,
    cloglog  = pGumbel,
    cauchit  = stats::pcauchy
  )
  cumpr <- matrix(pfun(zeta_mat - eta), nrow = nrow(x_mat), ncol = q)
  probs <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
  colnames(probs) <- model$lev

  if (!all(m_level_names %in% colnames(probs))) {
    stop(paste(strwrap("Error: The fitted ordered-logit model for the mediator did not return probabilities for all observed levels of M."), collapse = "\n"))
  }
  probs[, m_level_names, drop = FALSE]
}


ipwvent_inner <- function(
  data,
  D,
  M,
  L,
  Y,
  D_model,
  L_model,
  M_model,
  base_weights_name = NULL,
  censor = TRUE,
  censor_low = 0.01,
  censor_high = 0.99,
  minimal = FALSE
) {
  # load data
  df <- data

  # assign base weights
  if (is.null(base_weights_name)) {
    base_weights <- rep(1, nrow(df))
  }
  else {
    base_weights <- df[[base_weights_name]]
  }
  if (!is.null(base_weights_name) && any(is.na(base_weights))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the base weights variable (identified by the string argument base_weights_name in data). If that observation should not receive a positive weight, please replace the NA value with a zero before proceeding."), collapse = "\n"))
  }
  # BUG 0.12: reject negative base weights
  if (any(base_weights < 0, na.rm = TRUE)) {
    stop(paste(strwrap("Error: base_weights_name contains negative values. All weights must be non-negative."), collapse = "\n"))
  }

  # rescale base weights
  base_weights_rsc <- base_weights / mean(base_weights)

  # reset formula environments so glm/polr can find locally-defined objects
  # (e.g. base_weights_rsc). R copy-on-modify keeps the caller's objects unchanged.
  f_D <- formula(D_model)
  environment(f_D) <- environment()
  f_M <- formula(M_model)
  environment(f_M) <- environment()

  l_values <- .ipwvent_get_l_values(df[[L]])
  if (length(l_values) < 2L) {
    stop(paste(strwrap("Error: The exposure-induced confounder L must have at least two observed levels in the estimation sample."), collapse = "\n"))
  }

  # fit exposure model
  d_model <- glm(
    f_D,
    data = df,
    family = D_model$family,
    weights = base_weights_rsc
  )

  # fit model for the exposure-induced confounder
  if (inherits(L_model, "polr")) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop(paste(strwrap("Error: The ordered-logit model for the exposure-induced confounder requires the MASS package, but MASS is not installed."), collapse = "\n"))
    }
    df[[".ipwvent_L"]] <- .ipwvent_prepare_l_response(df[[L]], L_model = L_model)
    l_model <- .ipwvent_refit(L_model, df, response_var = ".ipwvent_L",
                              base_weights_rsc = base_weights_rsc)
  } else {
    l_model <- .ipwvent_refit(L_model, df, base_weights_rsc = base_weights_rsc)
  }

  # fit mediator model -- model type encoded in M_model class (glm or polr)
  # Determine M levels from data (re-computed inside inner for bootstrap safety)
  m_vals_inner <- sort(unique(df[[M]]))

  if (inherits(M_model, "polr")) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop(paste(strwrap("Error: The ordered-logit model for the mediator requires the MASS package, but MASS is not installed."), collapse = "\n"))
    }
    df_m <- df
    df_m[[".ipwvent_M"]] <- ordered(df_m[[M]], levels = m_vals_inner)
    m_model <- .ipwvent_refit(M_model, df_m, response_var = ".ipwvent_M",
                              base_weights_rsc = base_weights_rsc)
  } else {
    m_model <- glm(
      f_M,
      data = df,
      family = M_model$family,
      weights = base_weights_rsc
    )
  }

  # enforcing a no-missing-data rule
  if (nrow(d_model$model) != nrow(df) |
      nrow(l_model$model) != nrow(df) |
      nrow(m_model$model) != nrow(df)) {
    stop(paste(strwrap("Error: Please remove observations with missing/NA values for the exposure, mediator, exposure-induced confounder, outcome, or covariates."), collapse = "\n"))
  }

  # predict exposure probabilities
  ps_D1_C <- predict(d_model, newdata = df, type = "response")
  ps_D0_C <- 1 - ps_D1_C
  # Marginal P(D) -- empirical weighted mean
  marg_prob_D1 <- weighted.mean(df[[D]], base_weights_rsc)
  marg_prob_D0 <- 1 - marg_prob_D1

  # identify support groups (hardcoded: d=1, dstar=0)
  group_dstar <- df[[D]] == 0 & base_weights_rsc > 0
  group_d     <- df[[D]] == 1 & base_weights_rsc > 0

  # create counterfactual data sets used in the slide algorithm
  idataD0 <- df
  idataD0[[D]] <- 0
  idataD1 <- df
  idataD1[[D]] <- 1

  # predict the distribution of L under D = 1 and D = 0
  p_L_D1C <- .ipwvent_predict_l_probs(
    model = l_model,
    newdata = idataD1,
    l_values = l_values
  )
  p_L_D0C <- .ipwvent_predict_l_probs(
    model = l_model,
    newdata = idataD0,
    l_values = l_values
  )

  # predict mediator probability mass at the observed M value.
  # For binary M (glm):   P(M = m_obs | D, L, C) via predict(glm, type="response")
  # For ordinal M (polr): P(M = m_obs | D, L, C) via polr probs column extraction
  m_obs_char <- as.character(df[[M]])

  if (!inherits(m_model, "polr")) {
    # P(M=1|D=1,L_obs,C) and P(M=1|D=0,L_obs,C)
    p1_D1    <- predict(m_model, newdata = idataD1, type = "response")
    p_M_D1LC <- ifelse(df[[M]] == m_vals_inner[2L], p1_D1, 1 - p1_D1)
    p1_D0    <- predict(m_model, newdata = idataD0, type = "response")
    p_M_D0LC <- ifelse(df[[M]] == m_vals_inner[2L], p1_D0, 1 - p1_D0)
    col_idx  <- NULL  # not used for glm
  } else {
    # ordinal M: pre-compute column index for observed M (reused in loop)
    probs_D1 <- .ipwvent_predict_m_polr_probs(m_model, idataD1, m_vals_inner)
    col_idx  <- match(m_obs_char, colnames(probs_D1))
    p_M_D1LC <- probs_D1[cbind(seq_len(nrow(probs_D1)), col_idx)]

    probs_D0 <- .ipwvent_predict_m_polr_probs(m_model, idataD0, m_vals_inner)
    p_M_D0LC <- probs_D0[cbind(seq_len(nrow(probs_D0)), col_idx)]
  }

  # numerator: sum_l P(M=m_obs | D=d, L=l, C) * P(L=l | D=d, C)
  numer_D1 <- rep(0, nrow(df))
  numer_D0 <- rep(0, nrow(df))
  for (j in seq_along(l_values)) {
    l_value_j <- l_values[j]
    l_name_j  <- as.character(l_value_j)

    idataD0Lj <- idataD0
    idataD0Lj[[L]] <- .ipwvent_set_l_value(df[[L]], value = l_value_j, n = nrow(df))
    idataD1Lj <- idataD1
    idataD1Lj[[L]] <- .ipwvent_set_l_value(df[[L]], value = l_value_j, n = nrow(df))

    if (!inherits(m_model, "polr")) {
      p1_D0Lj   <- predict(m_model, newdata = idataD0Lj, type = "response")
      p_M_D0LjC <- ifelse(df[[M]] == m_vals_inner[2L], p1_D0Lj, 1 - p1_D0Lj)
      p1_D1Lj   <- predict(m_model, newdata = idataD1Lj, type = "response")
      p_M_D1LjC <- ifelse(df[[M]] == m_vals_inner[2L], p1_D1Lj, 1 - p1_D1Lj)
    } else {
      pr_D0Lj   <- .ipwvent_predict_m_polr_probs(m_model, idataD0Lj, m_vals_inner)
      p_M_D0LjC <- pr_D0Lj[cbind(seq_len(nrow(pr_D0Lj)), col_idx)]
      pr_D1Lj   <- .ipwvent_predict_m_polr_probs(m_model, idataD1Lj, m_vals_inner)
      p_M_D1LjC <- pr_D1Lj[cbind(seq_len(nrow(pr_D1Lj)), col_idx)]
    }

    numer_D0 <- numer_D0 + (p_M_D0LjC * p_L_D0C[, l_name_j])
    numer_D1 <- numer_D1 + (p_M_D1LjC * p_L_D1C[, l_name_j])
  }

  # create inverse probability weights exactly as on the lecture slides
  # sw1: targets E[Y(0, G_0)] -- applied to D=0 observations
  sw1 <- ifelse(
    group_dstar,
    (marg_prob_D0 / ps_D0_C) *
      (1 / p_M_D0LC) *
      numer_D0,
    0
  )
  # sw2: targets E[Y(1, G_1)] -- applied to D=1 observations
  sw2 <- ifelse(
    group_d,
    (marg_prob_D1 / ps_D1_C) *
      (1 / p_M_D1LC) *
      numer_D1,
    0
  )
  # sw3: targets E[Y(1, G_0)] -- applied to D=1 observations, dstar mediator distribution
  sw3 <- ifelse(
    group_d,
    (marg_prob_D1 / ps_D1_C) *
      (1 / p_M_D1LC) *
      numer_D0,
    0
  )

  # censor IPWs among observations with non-zero weights
  if (censor) {
    sw1[group_dstar] <- trimQ(sw1[group_dstar], low = censor_low, high = censor_high)
    sw2[group_d] <- trimQ(sw2[group_d], low = censor_low, high = censor_high)
    sw3[group_d] <- trimQ(sw3[group_d], low = censor_low, high = censor_high)
  }

  # multiply IPWs by rescaled base weights
  final_w1 <- sw1 * base_weights_rsc
  final_w2 <- sw2 * base_weights_rsc
  final_w3 <- sw3 * base_weights_rsc

  # estimate interventional means and effects
  Ehat_Ydstar_Gdstar <- weighted.mean(df[[Y]], final_w1)
  Ehat_Yd_Gd <- weighted.mean(df[[Y]], final_w2)
  Ehat_Yd_Gdstar <- weighted.mean(df[[Y]], final_w3)
  ATE <- Ehat_Yd_Gd - Ehat_Ydstar_Gdstar
  IDE <- Ehat_Yd_Gdstar - Ehat_Ydstar_Gdstar
  IIE <- Ehat_Yd_Gd - Ehat_Yd_Gdstar

  # compile and output
  if (minimal) {
    out <- list(
      ATE = ATE,
      IDE = IDE,
      IIE = IIE
    )
  }
  else {
    out <- list(
      ATE = ATE,
      IDE = IDE,
      IIE = IIE,
      weights1 = final_w1,
      weights2 = final_w2,
      weights3 = final_w3,
      model_d = d_model,
      model_l = l_model,
      model_m = m_model,
      Ehat_Ydstar_Gdstar = Ehat_Ydstar_Gdstar,
      Ehat_Yd_Gd = Ehat_Yd_Gd,
      Ehat_Yd_Gdstar = Ehat_Yd_Gdstar
    )
  }
  return(out)
}


#' Interventional-effects inverse probability weighting (IPW) estimator
#'
#' @description
#' `ipwvent()` uses inverse probability weighting (IPW) to estimate the total
#' effect (ATE), interventional direct effect (IDE), and interventional
#' indirect effect (IIE).
#'
#' @details
#' `ipwvent()` performs causal mediation analysis using inverse probability
#' weighting and computes inferential statistics using the nonparametric
#' bootstrap. The function is designed for settings with a binary exposure `D`,
#' a discrete mediator `M` (binary or ordered categorical), and a single
#' exposure-induced confounder `L` of the mediator-outcome relation.
#'
#' To construct the weights, `ipwvent()` fits three nuisance models:
#'
#' 1. A logit model for the exposure conditional on baseline covariates,
#'    \eqn{P(D = 1 \mid C)}.
#' 2. A model for the discrete exposure-induced confounder conditional on the
#'    exposure and baseline covariates, \eqn{P(L = l \mid D, C)} for each
#'    observed level \eqn{l} of `L`.
#' 3. A logistic model (binary M) or ordered logistic model (ordinal M) for
#'    \eqn{P(M = m \mid D, L, C)}.
#'
#' These fitted models are combined to estimate three counterfactual means:
#' \deqn{E[Y_{0, \tilde{G}_{0}}], \; E[Y_{1, \tilde{G}_{1}}], \; E[Y_{1, \tilde{G}_{0}}],}
#' where \eqn{\tilde{G}_{a}} denotes the mediator distribution under \eqn{D = a}
#' after integrating over the exposure-induced confounder distribution under
#' \eqn{D = a} given baseline covariates.
#'
#' The resulting decomposition is
#' \deqn{ATE(1,0) = IDE(1,0) + IIE(1,0)}
#' with
#' \deqn{IDE(1,0) = E[Y_{1, \tilde{G}_{0}}] - E[Y_{0, \tilde{G}_{0}}]}
#' and
#' \deqn{IIE(1,0) = E[Y_{1, \tilde{G}_{1}}] - E[Y_{1, \tilde{G}_{0}}].}
#'
#' @param data A data frame.
#' @param D A character scalar identifying the name of the exposure variable in
#'   `data`. The exposure must be numeric, binary, and coded 0/1.
#' @param M A character scalar identifying the name of the mediator variable in
#'   `data`. M must be discrete: either binary (coded 0/1) or ordered
#'   categorical (integer values, 3--20 levels). Continuous, non-integer, or
#'   standardised mediators are rejected.
#' @param L A character scalar identifying the name of the exposure-induced
#'   confounder variable in `data`.
#' @param Y A character scalar identifying the name of the outcome variable in
#'   `data`. The outcome variable must be numeric.
#' @param D_model A fitted \code{glm()} object for the exposure given baseline
#'   covariates (binomial or quasibinomial family). The response variable must
#'   match the \code{D} argument.
#'   E.g., \code{D_model = glm(att22 ~ female + black + paredu, data = nlsy2, family = binomial())}.
#' @param L_model A fitted \code{glm()} (binomial or quasibinomial family) or
#'   \code{MASS::polr()} object for the time-varying confounder given exposure
#'   and baseline covariates. The response variable must match the \code{L}
#'   argument. \code{lm()} is not accepted because \code{ipwvent()} sums over
#'   discrete levels of L.
#' @param M_model A fitted \code{glm()} (binomial or quasibinomial family) or
#'   \code{MASS::polr()} object for the mediator given exposure, time-varying
#'   confounder, and baseline covariates. The response variable must match the
#'   \code{M} argument. \code{lm()} is not accepted because \code{ipwvent()}
#'   sums over discrete levels of M.
#' @param base_weights_name A character scalar identifying the name of the base
#'   weights variable in `data`, if applicable.
#' @param censor A logical scalar indicating whether the IPW weights should be
#'   censored.
#' @param censor_low,censor_high Numeric scalars in \eqn{[0, 1]} specifying the
#'   quantile bounds for weight censoring.
#' @param boot A logical scalar indicating whether to perform the nonparametric
#'   bootstrap.
#' @param boot_reps An integer scalar for the number of bootstrap replications. In practice, we recommend a minimum of 1000 replications.
#' @param boot_conf_level A numeric scalar for the confidence level of the
#'   bootstrap interval.
#' @param boot_seed An integer scalar specifying the random-number seed.
#' @param boot_parallel A logical scalar for parallelised bootstrapping.
#' @param boot_cores An integer scalar for the number of CPU cores.
#'
#' @returns A list with elements ATE, IDE, IIE, weights1--3, model_d, model_l,
#'   model_m, and (if boot = TRUE) ci_ATE, ci_IDE, ci_IIE, pvalue_ATE,
#'   pvalue_IDE, pvalue_IIE, boot_ATE, boot_IDE, boot_IIE.
#'
#' @examples
#' \dontrun{
#' ## Prepare data
#' load(url(paste0(
#'   "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/",
#'   "refs/heads/main/data/Brader_et_al2008/Brader_et_al2008.RData"
#' )))
#'
#' D <- "tone_eth"
#' M <- "emo_ord"
#' L <- "p_harm_ord"
#' Y <- "std_immigr"
#' C <- c("ppage", "female", "hs", "sc", "ba", "ppincimp")
#'
#' # keep complete cases for the raw source columns; standardize outcome
#' key_vars <- c("immigr", "emo", "p_harm", D, C)
#' Brader1 <- Brader[complete.cases(Brader[, key_vars]), ]
#' Brader1$std_immigr <-
#'   (Brader1$immigr - mean(Brader1$immigr)) / sd(Brader1$immigr)
#'
#' # create ordered-factor versions of L and M
#' Brader1$p_harm_ord <- ordered(Brader1$p_harm)
#' Brader1$emo_ord    <- ordered(Brader1$emo)
#'
#' # helper to build the RHS of a formula from a character vector
#' rhs <- function(xs) paste(xs, collapse = " + ")
#'
#' ## Fit nuisance models
#' D_model <- glm(
#'   as.formula(paste(D, "~", rhs(C))),
#'   data   = Brader1,
#'   family = binomial()
#' )
#' L_model <- MASS::polr(
#'   as.formula(paste(L, "~", rhs(c(D, C)))),
#'   data = Brader1,
#'   Hess = TRUE
#' )
#' M_model <- MASS::polr(
#'   as.formula(paste(M, "~", rhs(c(D, L, C)))),
#'   data = Brader1,
#'   Hess = TRUE
#' )
#'
#' ## Estimate interventional effects (point estimates)
#' out1 <- ipwvent(
#'   data    = Brader1,
#'   D       = D,
#'   M       = M,
#'   L       = L,
#'   Y       = Y,
#'   D_model = D_model,
#'   L_model = L_model,
#'   M_model = M_model
#' )
#' out1[c("ATE", "IDE", "IIE")]
#'
#' ## Perform a nonparametric bootstrap
#' out2 <- ipwvent(
#'   data      = Brader1,
#'   D         = D,
#'   M         = M,
#'   L         = L,
#'   Y         = Y,
#'   D_model   = D_model,
#'   L_model   = L_model,
#'   M_model   = M_model,
#'   boot      = TRUE,
#'   boot_reps = 1000,
#'   boot_seed = 1234
#' )
#' out2[c(
#'   "ATE", "IDE", "IIE",
#'   "ci_ATE", "ci_IDE", "ci_IIE",
#'   "pvalue_ATE", "pvalue_IDE", "pvalue_IIE"
#' )]
#' }
#'
#' @export
ipwvent <- function(
  data,
  D,
  M,
  L,
  Y,
  D_model,
  L_model,
  M_model,
  base_weights_name = NULL,
  censor = TRUE,
  censor_low = 0.01,
  censor_high = 0.99,
  boot = FALSE,
  boot_reps = 200,
  boot_conf_level = 0.95,
  boot_seed = NULL,
  boot_parallel = FALSE,
  boot_cores = max(c(parallel::detectCores() - 2L, 1L), na.rm = TRUE)
) {
  # load data
  data_outer <- data

  # create adjusted boot_parallel logical
  boot_parallel_rev <- isTRUE(boot_parallel) && boot_cores > 1

  # preliminary error/warning checks for the bootstrap
  if (boot) {
    # BUG 0.6: validate boot_reps
    if (!is.numeric(boot_reps) || length(boot_reps) != 1 ||
        boot_reps < 2 || boot_reps != round(boot_reps))
      stop("Error: boot_reps must be a single integer >= 2.")
    # BUG 0.8: validate boot_conf_level
    if (!is.numeric(boot_conf_level) ||
        boot_conf_level <= 0 || boot_conf_level >= 1)
      stop("Error: boot_conf_level must be a numeric value in (0, 1).")
    if (boot_parallel && boot_cores == 1) {
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

  # BUG 0.7: validate censor bounds
  if (censor) {
    if (!is.numeric(censor_low)  || censor_low  <  0 ||
        !is.numeric(censor_high) || censor_high >  1 ||
        censor_low >= censor_high)
      stop(paste(strwrap("Error: censor_low and censor_high must satisfy 0 <= censor_low < censor_high <= 1. Note: unlike Stata (which takes integer percentiles such as 1 99), ipwvent() uses proportions such as 0.01 0.99."), collapse = "\n"))
  }

  # BUG 0.5: check that D, M, L, Y and (optionally) base_weights_name exist
  for (v in c(D, M, L, Y)) {
    if (!v %in% names(data_outer))
      stop(paste0("Error: '", v, "' was not found as a column in data."))
  }
  if (!is.null(base_weights_name) && !base_weights_name %in% names(data_outer))
    stop(paste0("Error: '", base_weights_name,
                "' (base_weights_name) was not found in data."))

  # other error/warning checks
  # --- model-object validation ---
  if (!inherits(D_model, "glm"))
    stop("D_model must be a glm() object fit with binomial or quasibinomial family.")
  if (!D_model$family$family %in% c("binomial", "quasibinomial"))
    stop("D_model must use binomial or quasibinomial family (logistic regression for binary exposure).")

  if (!inherits(L_model, "glm") && !inherits(L_model, "polr"))
    stop("L_model must be a glm() or MASS::polr() object. lm() is not supported: ipwvent() sums over discrete levels of L.")
  if (inherits(L_model, "glm") && !L_model$family$family %in% c("binomial", "quasibinomial"))
    stop("L_model glm() must use binomial or quasibinomial family.")

  if (!inherits(M_model, "glm") && !inherits(M_model, "polr"))
    stop("M_model must be a glm() or MASS::polr() object. lm() is not supported: ipwvent() sums over discrete levels of M.")
  if (inherits(M_model, "glm") && !M_model$family$family %in% c("binomial", "quasibinomial"))
    stop("M_model glm() must use binomial or quasibinomial family.")

  # check that each model targets the intended response
  if (all.vars(formula(D_model))[1] != D)
    stop(paste(strwrap("Error: The response variable in D_model does not match the D argument."), collapse = "\n"))
  if (all.vars(formula(L_model))[1] != L)
    stop(paste(strwrap("Error: The response variable in L_model does not match the L argument."), collapse = "\n"))
  if (all.vars(formula(M_model))[1] != M)
    stop(paste(strwrap("Error: The response variable in M_model does not match the M argument."), collapse = "\n"))

  # --- data variable checks ---
  if (!is.numeric(data_outer[[D]])) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be numeric."), collapse = "\n"))
  }
  if (!is.numeric(data_outer[[Y]])) {
    stop(paste(strwrap("Error: The outcome variable (identified by the string argument Y in data) must be numeric."), collapse = "\n"))
  }
  if (!(is.numeric(data_outer[[L]]) || is.factor(data_outer[[L]]) ||
        is.ordered(data_outer[[L]]) || is.character(data_outer[[L]]))) {
    stop(paste(strwrap("Error: The exposure-induced confounder variable L must be numeric, character, a factor, or an ordered factor."), collapse = "\n"))
  }
  if (any(is.na(data_outer[[D]]))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the exposure variable (identified by the string argument D in data)."), collapse = "\n"))
  }
  if (any(!data_outer[[D]] %in% c(0, 1))) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be a numeric variable consisting only of the values 0 or 1. There is at least one observation in the data that does not meet this criterion."), collapse = "\n"))
  }
  # BUG 0.2: reject continuous L (non-integer numeric values)
  if (is.numeric(data_outer[[L]])) {
    l_vals_check <- sort(unique(data_outer[[L]]))
    if (!all(l_vals_check == round(l_vals_check))) {
      stop(paste(strwrap("Error: The exposure-induced confounder L must be discrete (binary or ordered categorical). Continuous or non-integer numeric values are not supported."), collapse = "\n"))
    }
    if (inherits(L_model, "glm") && !(length(l_vals_check) == 2L && setequal(l_vals_check, c(0, 1)))) {
      stop(paste(strwrap("Error: When L_model is a glm(), L must be binary and coded 0/1."), collapse = "\n"))
    }
  }
  if (!is.numeric(data_outer[[M]]) && !is.ordered(data_outer[[M]])) {
    stop(paste(strwrap("Error: The mediator variable M must be numeric or an ordered factor."), collapse = "\n"))
  }
  # BUG 0.17: reject continuous or non-integer M; require discrete M with <= 20 levels
  if (is.numeric(data_outer[[M]])) {
    m_vals_check <- sort(unique(data_outer[[M]]))
    m_gaps <- if (length(m_vals_check) > 1L) diff(m_vals_check) else numeric(0)
    m_is_integer_valued <- all(m_vals_check == round(m_vals_check))
    m_is_equally_spaced <- length(m_gaps) == 0L || diff(range(m_gaps)) < 1e-8
    if (!m_is_integer_valued || !m_is_equally_spaced || length(m_vals_check) > 20L) {
      stop(paste(strwrap("Error: The mediator M must be discrete -- either binary (2 levels) or ordered categorical (3 or more ordered integer levels). Continuous or non-integer numeric mediators are not supported. ipwvent() requires discrete M because the weight numerator sums over levels of M analogously to L. If M is a Likert-scale or similar ordinal variable, ensure it has not been mean-centred or standardised before passing it to ipwvent()."), collapse = "\n"))
    }
    if (m_is_integer_valued && m_is_equally_spaced &&
        length(m_vals_check) > 10L && length(m_vals_check) <= 20L) {
      warning(paste(strwrap("Warning: The mediator M has more than 10 ordered levels. If M is a count variable rather than an ordered categorical variable, results may be unreliable. ipwvent() is designed for binary or ordered categorical M."), collapse = "\n"))
    }
    if (length(m_vals_check) == 2L && !setequal(m_vals_check, c(0, 1))) {
      stop(paste(strwrap("Error: Binary M (exactly 2 unique values) must be coded 0/1."), collapse = "\n"))
    }
  }

  # compute point estimates
  est <- ipwvent_inner(
    data = data_outer,
    D = D,
    M = M,
    L = L,
    Y = Y,
    D_model = D_model,
    L_model = L_model,
    M_model = M_model,
    base_weights_name = base_weights_name,
    censor = censor,
    censor_low = censor_low,
    censor_high = censor_high,
    minimal = FALSE
  )

  # bootstrap, if requested
  if (boot) {
    boot_fnc <- function() {
      boot_data <- data_outer[sample(nrow(data_outer), size = nrow(data_outer), replace = TRUE), ]

      ipwvent_inner(
        data = boot_data,
        D = D,
        M = M,
        L = L,
        Y = Y,
        D_model = D_model,
        L_model = L_model,
        M_model = M_model,
        base_weights_name = base_weights_name,
        censor = censor,
        censor_low = censor_low,
        censor_high = censor_high,
        minimal = TRUE
      )
    }

    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type = "PSOCK")
      doParallel::registerDoParallel(cl = x_cluster)
      parallel::clusterExport(
        cl = x_cluster,
        varlist = c(
          "ipwvent_inner",
          ".ipwvent_build_rhs_formula",
          ".ipwvent_build_polr_start",
          ".ipwvent_fit_polr",
          ".ipwvent_refit",
          ".ipwvent_get_l_values",
          ".ipwvent_prepare_l_response",
          ".ipwvent_set_l_value",
          ".ipwvent_as_prob_matrix",
          ".ipwvent_predict_l_probs",
          ".ipwvent_predict_m_polr_probs",
          "trimQ",
          "D_model",
          "L_model",
          "M_model"
        ),
        envir = environment()
      )
      `%dopar%` <- foreach::`%dopar%`
    }

    if (!is.null(boot_seed)) {
      set.seed(boot_seed)
      # BUG 0.10: only register doRNG when a real cluster was actually created
      if (boot_parallel_rev) {
        doRNG::registerDoRNG(boot_seed)
      }
    }

    if (boot_parallel_rev) {
      boot_res <- foreach::foreach(i = 1:boot_reps, .combine = comb_list_vec) %dopar% {
        boot_fnc()
      }
      boot_ATE <- boot_res$ATE
      boot_IDE <- boot_res$IDE
      boot_IIE <- boot_res$IIE
    }
    else {
      boot_ATE <- rep(NA_real_, boot_reps)
      boot_IDE <- rep(NA_real_, boot_reps)
      boot_IIE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_iter <- tryCatch(boot_fnc(), error = function(e) NULL)
        if (!is.null(boot_iter)) {
          boot_ATE[i] <- boot_iter$ATE
          boot_IDE[i] <- boot_iter$IDE
          boot_IIE[i] <- boot_iter$IIE
        }
      }
      n_fail <- sum(is.na(boot_ATE))
      if (n_fail > 0L) {
        warning(paste(strwrap(paste0(
          "Warning: ", n_fail, " of ", boot_reps,
          " bootstrap iterations failed to converge (model optimisation error) ",
          "and were excluded from CI and p-value computation. ",
          "Results are based on ", boot_reps - n_fail, " successful replications."
        )), collapse = "\n"))
      }
    }

    if (boot_parallel_rev) {
      parallel::stopCluster(x_cluster)
      rm(x_cluster)
    }

    boot_alpha <- 1 - boot_conf_level
    boot_ci_probs <- c(
      boot_alpha / 2,
      1 - boot_alpha / 2
    )
    boot_ci <- function(x) {
      quantile(x, probs = boot_ci_probs, na.rm = TRUE)
    }
    ci_ATE <- boot_ci(boot_ATE)
    ci_IDE <- boot_ci(boot_IDE)
    ci_IIE <- boot_ci(boot_IIE)

    boot_pval <- function(x) {
      2 * min(
        mean(x < 0, na.rm = TRUE),
        mean(x > 0, na.rm = TRUE)
      )
    }
    pvalue_ATE <- boot_pval(boot_ATE)
    pvalue_IDE <- boot_pval(boot_IDE)
    pvalue_IIE <- boot_pval(boot_IIE)

    boot_out <- list(
      ci_ATE = ci_ATE,
      ci_IDE = ci_IDE,
      ci_IIE = ci_IIE,
      pvalue_ATE = pvalue_ATE,
      pvalue_IDE = pvalue_IDE,
      pvalue_IIE = pvalue_IIE,
      boot_ATE = boot_ATE,
      boot_IDE = boot_IDE,
      boot_IIE = boot_IIE
    )
  }

  out <- est
  if (boot) {
    out <- append(out, boot_out)
  }
  return(out)
}
