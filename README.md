# causalMedR
This repository contains R functions for conducting causal mediation analysis using the methods described in Wodtke and Zhou "Causal Mediation Analysis."

## Table of Contents
- [linmed – mediation analysis using linear models](#linmed-mediation-analysis-using-linear-models)
- [medsim – mediation analysis using a simulation approach](#medsim-mediation-analysis-using-a-simulation-approach)
- [ipwmed – mediation analysis using inverse probability weights](#ipwmed-mediation-analysis-using-inverse-probability-weights)
- [impcde – a regression imputation estimator for controlled direct effects](#impcde-a-regression-imputation-estimator-for-controlled-direct-effects)
- [ipwcde – an inverse probability weighting estimator for controlled direct effects](#ipwcde-an-inverse-probability-weighting-estimator-for-controlled-direct-effects)
- [rwrlite – regression-with-residuals estimation for interventional effects](#rwrlite-regression-with-residuals-estimation-for-interventional-effects)
- [linpath – analysis of path-specific effects using linear models](#linpath-analysis-of-path-specific-effects-using-linear-models)
- [ipwpath – analysis of path-specific effects using inverse probability weights](#ipwpath-analysis-of-path-specific-effects-using-inverse-probability-weights)
- [utils – utility functions](#utils-utility-functions)

## `linmed`: mediation analysis using linear models

The `linmed` function estimates natural direct, natural indirect, controlled direct, and total effects using linear models. It supports mediation analysis with a single mediator and with multiple mediators, and includes built-in support for nonparametric bootstrapping with optional parallel computation.

### Function

```r
linmed(
  data,
  D,
  M,
  Y,
  C = NULL,
  d,
  dstar,
  m = NULL,
  interaction_DM = FALSE,
  interaction_DC = FALSE,
  interaction_MC = FALSE,
  weights_name = NULL,
  boot = FALSE,
  boot_reps = 1000,
  boot_conf_level = 0.95,
  boot_seed = NULL,
  boot_parallel = FALSE,
  boot_cores = NULL
)
```

### Arguments

| Argument | Description |
|----------|-------------|
| `data` | A data frame. |
| `D` | Name of the exposure variable (must be numeric). |
| `M` | Name(s) of mediator variable(s). Can be a character vector of one or more numeric variables. |
| `Y` | Name of the outcome variable (must be numeric). |
| `C` | (Optional) Covariates to include in both mediator and outcome models. |
| `d`, `dstar` | Numeric values specifying the exposure contrast of interest (`d - dstar`). |
| `m` | (Optional) Values for mediators when estimating the controlled direct effect (CDE). Must be the same length as `M`. |
| `interaction_DM` | Whether to include exposure × mediator interactions in the outcome model. |
| `interaction_DC` | Whether to include exposure × covariate interactions in **both** mediator and outcome models. |
| `interaction_MC` | Whether to include mediator × covariate interactions in the outcome model. |
| `weights_name` | (Optional) Name of the weights variable, if using sampling weights. |
| `boot` | Whether to perform nonparametric bootstrap. |
| `boot_reps` | Number of bootstrap replications (default: 1000). |
| `boot_conf_level` | Confidence level for bootstrap intervals (default: 0.95). |
| `boot_seed` | Random seed for reproducibility. |
| `boot_parallel` | Whether to parallelize bootstrap (requires `doParallel`, `doRNG`, `foreach`). |
| `boot_cores` | Number of CPU cores to use when parallelizing. |

### Returns

A list with the following elements:

- `ATE`: Estimated total average treatment effect.
- `NDE`: Estimated natural direct effect.
- `NIE`: Estimated natural indirect effect.
- `CDE`: Estimated controlled direct effect (for specified `m`).
- `model_m`: List of fitted mediator model objects.
- `model_y`: Fitted outcome model object.
- `miss_summary`: Summary of missing vs. non-missing observations.

If `boot = TRUE`, the list also includes:

- `ci_ATE`, `ci_NDE`, `ci_NIE`, `ci_CDE`: Bootstrap confidence intervals.
- `pvalue_ATE`, `pvalue_NDE`, `pvalue_NIE`, `pvalue_CDE`: Bootstrap-based p-values from tests of no effect.
- `boot_ATE`, `boot_NDE`, `boot_NIE`, `boot_CDE`: Bootstrap replicate estimates.

### Examples

#### Single Mediator, No Interactions

```r
linmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0
)
```

#### With Interaction Terms

```r
linmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  interaction_DM = TRUE,
  interaction_DC = TRUE,
  interaction_MC = TRUE,
  d = 1,
  dstar = 0
)
```

#### Controlled Direct Effect (CDE)

```r
linmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0,
  m = 1
)
```

#### Multiple Mediators

```r
linmed(
  data = nlsy,
  D = "att22",
  M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0
)
```

#### Bootstrapped Estimates

```r
linmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0,
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234
)
```

#### Parallelized Bootstrap

```r
linmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0,
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234,
  boot_parallel = TRUE
)
```


## `medsim`: mediation analysis using a simulation approach

The `medsim` function estimates natural direct and indirect effects, interventional direct and indirect effects, controlled direct effects, and path-specific effects using a simulation approach. It supports a wide variety of models (including the entire family of GLMs, multinomial logit, and ordered logistic), and includes optional support for bootstrapping with parallel processing.

### Function

```r
medsim(
  data,
  num_sim = 2000,
  cat_list,
  treatment,
  intv_med,
  model_spec,
  weights = NULL,
  boot = FALSE,
  reps = 100,
  resv_core = 1,
  seed = NULL
)
```

### Arguments

| Argument      | Description |
|---------------|-------------|
| `data`        | A data frame containing all variables referenced in the model specifications. |
| `num_sim`     | Integer. Number of Monte Carlo simulation draws (default: `2000`). |
| `cat_list`    | A vector of treatment levels to compare (default: `c(0, 1)`). |
| `treatment`   | Name of the treatment variable (character). |
| `intv_med`    | A list specifying the intervention(s) on the mediators. Set to `NULL` if interventional or controlled direct effects are not of interest. |
| `model_spec`  | A list of lists defining the models for the mediators and the outcome. Each model must include: <br> • `func`: model-fitting function (e.g., `"glm"`, `"polr"`) <br> • `formula`: model formula <br> • `args`: (optional) list of additional arguments to pass to the function. |
| `weights`     | (Optional) Name of the variable containing weights to use in model fitting. If `NULL`, no weights are applied. |
| `boot`        | Logical. If `TRUE`, performs nonparametric bootstrap to obtain confidence intervals and p-values (default: `FALSE`). Requires `doParallel`, `doRNG`, and `foreach`. |
| `reps`        | Integer. Number of bootstrap replications (default: `100`). |
| `resv_core`   | Integer. Number of CPU cores to reserve (i.e., not use) when bootstrapping via parallel processing (default: `1`). |
| `seed`        | Integer or `NULL`. Seed for reproducibility. |

### Returns

- If `boot = FALSE`:  
  A list of point estimates for the mediation effects of interest.

- If `boot = TRUE`:  
  A data frame containing:
  - Point estimates for the mediation effects of interest  
  - 95% bootstrap confidence intervals  
  - P-values from tests of no effect

### Examples

#### Estimate Natural Direct and Indirect Effects

```r
# Specify models for M (logit) and Y (normal linear)
formula_M <- paste(M, "~", paste(c(D, C), collapse = " + "))
formula_Y <- paste(Y, "~", paste(c(D, M, C), collapse = " + "))

specs <- list(
  list(func = "glm", formula = as.formula(formula_M), args = list(family = "binomial")),
  list(func = "lm", formula = as.formula(formula_Y))
)

# Compute estimates
sim_nat <- medsim(
  data = nlsy,
  num_sim = 1000,
  treatment = D,
  intv_med = NULL,
  model_spec = specs,
  seed = 60637
)
```

#### Estimate Interventional and Controlled Direct Effects, Include Bootstrap CIs and P-values

```r
# Specify models for the mediator, exposure-induced confounder, and the outcome
Lform <- ever_unemp_age3539 ~ att22 * (female + black + hispan +
  paredu + parprof + parinc_prank + famsize + afqt3)

Mform <- log_faminc_adj_age3539 ~ att22 * (female + black + hispan +
  paredu + parprof + parinc_prank + famsize + afqt3)

Yform <- std_cesd_age40 ~ (log_faminc_adj_age3539 * att22) *
  (female + black + hispan + paredu + parprof + parinc_prank +
   famsize + afqt3 + ever_unemp_age3539)

specs <- list(
  list(func = "glm", formula = as.formula(Lform), args = list(family = "binomial")),
  list(func = "lm", formula = as.formula(Mform)),
  list(func = "lm", formula = as.formula(Yform))
)

# Estimate interventional effects
sim_ie <- medsim(
  data = nlsy,
  num_sim = 1000,
  treatment = D,
  intv_med = M,
  model_spec = specs,
  boot = TRUE,
  reps = 2000,
  seed = 60637
)

# Estimate controlled direct effect setting M (income) to log(50000)
sim_cde <- medsim(
  data = nlsy,
  num_sim = 1000,
  treatment = D,
  intv_med = paste0(M, "=log(5e4)"),
  model_spec = specs,
  boot = TRUE,
  reps = 2000,
  seed = 60637
)
```


## `ipwmed`: mediation analysis using inverse probability weights

The `ipwmed` function implements the inverse probability weighting (IPW) estimator for the total effect, natural direct effect, and natural indirect effect. This function supports both univariate and multivariate mediators and includes built-in support for bootstrap inference and parallelized computation.

### Function

```r
ipwmed(
  data,
  D,
  M,
  Y,
  formula1_string,
  formula2_string,
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
  boot_cores = NULL
)
```

### Arguments

| Argument | Description |
|----------|-------------|
| `data` | A data frame. |
| `D` | Name of the exposure variable (character scalar; must identify a binary variable). |
| `M` | Name(s) of mediator variable(s); a character vector (length 1 or more) identifying numeric variables. |
| `Y` | Name of the outcome variable (character scalar; must identify a numeric variable). |
| `formula1_string` | A formula (as a string) for estimating the probability of exposure given baseline covariates using a logit model, i.e., _f(D \| C)_. |
| `formula2_string` | A formula (as a string) for estimating the probability of exposure given baseline covariates and mediators using a logit model, i.e., _s(D \| C, M)_. |
| `base_weights_name` | (Optional) Name of a variable containing sampling or base weights. |
| `stabilize` | Logical. If `TRUE`, uses stabilized weights (default: `TRUE`). |
| `censor` | Logical. If `TRUE`, applies weight censoring (default: `TRUE`). |
| `censor_low`, `censor_high` | Quantile cutoffs for censoring weights (default: 0.01 and 0.99, respectively). |
| `boot` | Logical. If `TRUE`, performs a bootstrap to return confidence intervals and p-values (default: `FALSE`). |
| `boot_reps` | Number of bootstrap replications (default: `1000`). |
| `boot_conf_level` | Confidence level for bootstrap intervals (default: `0.95`). |
| `boot_seed` | Integer seed for reproducibility. |
| `boot_parallel` | Logical. If `TRUE`, parallelizes the bootstrap (requires `doParallel`, `doRNG`, and `foreach`). |
| `boot_cores` | Number of CPU cores for parallel bootstrap. If `NULL`, defaults to available cores minus 2. |

### Returns

- If `boot = FALSE`:
  A **list** with the following elements:
  - `ATE`, `NDE`, `NIE`: Estimated effects
  - `weights1`, `weights2`, `weights3`: the inverse probability weights
  - `model_d1`, `model_d2`: Fitted logit models from `formula1_string` and `formula2_string`

- If `boot = TRUE`, the return includes the above, plus:
  - `ci_ATE`, `ci_NDE`, `ci_NIE`: Bootstrap confidence intervals
  - `pvalue_ATE`, `pvalue_NDE`, `pvalue_NIE`: Two-sided p-values
  - `boot_ATE`, `boot_NDE`, `boot_NIE`: Vectors of replicate bootstrap estimates

### Examples

#### Natural Effects Through a Single Mediator

```r
ipw_nat <- ipwmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
  formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3+ever_unemp_age3539"
)
```

#### Natural Effects Through Multiple Mediators

```r
ipw_mnat <- ipwmed(
  data = nlsy,
  D = "att22",
  M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
  formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3+ever_unemp_age3539+log_faminc_adj_age3539"
)
```

#### Bootstrap

```r
ipw_nat_boot <- ipwmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
  formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3+ever_unemp_age3539",
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234
)
```

#### Parallel Bootstrap

```r
ipw_nat_bootpar <- ipwmed(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  formula1_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3",
  formula2_string = "att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3+ever_unemp_age3539",
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234,
  boot_parallel = TRUE
)
```


## `impcde`: a regression imputation estimator for controlled direct effects

The `impcde` function estimates controlled direct effects using a regression imputation approach. It requires a fitted outcome model and computes the CDE by comparing predicted outcomes under different exposure levels, while holding the mediator fixed at a specified value. This function supports optional sampling weights and the nonparametric bootstrap for computing confidence intervals and p-values. Parallelized bootstrap computation is also supported.

### Function

```r
impcde(
  data,
  model_y,
  D,
  M,
  d,
  dstar,
  m,
  weights_name = NULL,
  boot = FALSE,
  boot_reps = 1000,
  boot_conf_level = 0.95,
  boot_seed = NULL,
  boot_parallel = FALSE,
  boot_cores = NULL
)
```

### Arguments

| Argument | Description |
|----------|-------------|
| `data` | A data frame containing the variables used in `model_y`, `D`, and `M`. |
| `model_y` | A fitted outcome model object (e.g., from `lm()` or `glm()`). |
| `D` | Name of the exposure variable (character scalar). |
| `M` | Name of the mediator variable (character scalar). |
| `d`, `dstar` | Numeric values representing two levels of the exposure. The exposure contrast of interest is `d - dstar`. |
| `m` | Numeric value at which to fix the mediator for CDE estimation. |
| `weights_name` | (Optional) Name of a variable containing sampling weights. |
| `boot` | Logical. If `TRUE`, use the nonparametric bootstrap to obtain confidence intervals and p-values. |
| `boot_reps` | Number of bootstrap replications (default: `1000`). |
| `boot_conf_level` | Confidence level for bootstrap interval (default: `0.95`). |
| `boot_seed` | Integer seed for reproducibility. |
| `boot_parallel` | Logical. If `TRUE`, parallelizes the bootstrap (requires `doParallel`, `doRNG`, and `foreach`). |
| `boot_cores` | Number of CPU cores for parallel bootstrap. If `NULL`, defaults to available cores minus 2. |

### Returns

- If `boot = FALSE`:  
  A numeric scalar representing the estimated CDE for the contrast `d - dstar` at mediator value `m`.

- If `boot = TRUE`:  
  A list with the following elements:
  - `CDE`: Point estimate
  - `ci_CDE`: Bootstrap confidence interval
  - `pvalue_CDE`: Two-sided bootstrap p-value
  - `boot_CDE`: Vector of bootstrap replicate estimates

### Examples

#### Basic Usage

```r
mod1 <- lm(
  std_cesd_age40 ~ att22 * (ever_unemp_age3539 +  + female + black + hispan + paredu +
    parprof + parinc_prank + famsize + afqt3),
  data = nlsy
)

impcde(
  data = nlsy,
  model_y = mod1,
  D = "att22",
  M = "ever_unemp_age3539",
  d = 1,
  dstar = 0,
  m = 1
)
```

#### Bootstrap 

```r
impcde(
  data = nlsy,
  model_y = mod1,
  D = "att22",
  M = "ever_unemp_age3539",
  d = 1,
  dstar = 0,
  m = 1,
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234
)
```

#### Parallel Bootstrap

```r
impcde(
  data = nlsy,
  model_y = mod1,
  D = "att22",
  M = "ever_unemp_age3539",
  d = 1,
  dstar = 0,
  m = 1,
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234,
  boot_parallel = TRUE
)
```


## `ipwcde`: an inverse probability weighting estimator for controlled direct effects

The `ipwcde` function estimates controlled direct effects using inverse probability weighting. This function supports only a single mediator, which must be numeric and binary. It supports optional weight stabilization, weight censoring, bootstrap confidence intervals, and parallelized bootstrap computation.

### Function

```r
ipwcde(
  data,
  D,
  M,
  Y,
  m,
  formula_D_string,
  formula_M_string,
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
  boot_cores = NULL
)
```

### Arguments

| Argument | Description |
|----------|-------------|
| `data` | A data frame. |
| `D` | Name of the exposure variable (character). Must be numeric and binary (0/1). |
| `M` | Name of the single mediator variable (character). Must be numeric and binary (0/1). |
| `Y` | Name of the numeric outcome variable (character). |
| `m` | Numeric value at which to fix the mediator for CDE estimation. |
| `formula_D_string` | A string representing the formula for a logit model for the exposure, e.g., `"att22 ~ female + black + paredu"` (used to estimate _f(D \| C)_). |
| `formula_M_string` | A string representing the formula for a logit model for the mediator, e.g., `"M ~ D + C"` (used to estimate _g(M \| C, D)_). |
| `base_weights_name` | (Optional) Name of the base weights variable. |
| `stabilize` | Logical. If `TRUE`, stabilizes the IPW weights (default: `TRUE`). |
| `censor` | Logical. If `TRUE`, applies weight censoring  (default: `TRUE`). |
| `censor_low`, `censor_high` | Quantiles for censoring the weights (default: 0.01 and 0.99). |
| `boot` | Logical. If `TRUE`, performs a bootstrap to compute CIs and p-values. |
| `boot_reps` | Number of bootstrap replications (default: `1000`). |
| `boot_conf_level` | Confidence level for bootstrap interval (default: `0.95`). |
| `boot_seed` | Integer seed for reproducibility. |
| `boot_parallel` | Logical. If `TRUE`, runs bootstrap in parallel (requires `doParallel`, `doRNG`, `foreach`). |
| `boot_cores` | Number of CPU cores for parallel bootstrap. Defaults to available cores minus two. |

### Returns

- If `boot = FALSE`:
  A list with:
  - `CDE`: Estimated controlled direct effect
  - `weights`: Final IPWs
  - `model_d`: Fitted logit model for the exposure (`formula_D_string`)
  - `model_m`: Fitted logit model for the mediator (`formula_M_string`)

- If `boot = TRUE`, the list also includes:
  - `ci_CDE`: Bootstrap confidence interval
  - `pvalue_CDE`: Two-sided p-value
  - `boot_CDE`: Vector of bootstrap replicate estimates

### Examples

#### Estimate the CDE

```r
cde_est <- ipwcde(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  m = 1,
  formula_D_string = "att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3",
  formula_M_string = "ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3"
)
```

#### Bootstrap

```r
cde_est_boot <- ipwcde(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  m = 1,
  formula_D_string = "att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3",
  formula_M_string = "ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3",
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234
)
```

#### Parallelized Bootstrap

```r
cde_est_bootpar <- ipwcde(
  data = nlsy,
  D = "att22",
  M = "ever_unemp_age3539",
  Y = "std_cesd_age40",
  m = 1,
  formula_D_string = "att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3",
  formula_M_string = "ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3",
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234,
  boot_parallel = TRUE
)
```


## `rwrlite`: regression-with-residuals estimation for interventional effects

The `rwrlite` function is a wrapper for two core functions in the [`rwrmed`](https://github.com/xiangzhou09/rwrmed) package. It implements regression-with-residuals (RWR) estimation to compute overall effects, interventional direct and indirect effects, and controlled direct effects.

To use this function, first install the `rwrmed` package using:

```r
devtools::install_github("xiangzhou09/rwrmed")
```

### Function

```r
rwrlite(
  data,
  D,
  C = NULL,
  d,
  dstar,
  m,
  Y_formula,
  M_formula,
  M_family,
  L_formula_list,
  weights = NULL,
  boot = FALSE,
  boot_reps = 1000,
  boot_conf_level = 0.95,
  boot_seed = NULL,
  boot_parallel = FALSE,
  boot_cores = NULL
)
```

### Arguments

| Argument | Description |
|----------|-------------|
| `data` | A data frame. |
| `D` | Name of the exposure variable (character scalar). |
| `C` | Covariates to include in all models (character vector). |
| `d`, `dstar` | Two numeric values of the exposure; the contrast of interest is `d - dstar`. |
| `m` | Numeric value at which to fix the mediator when estimating the CDE. |
| `Y_formula` | Formula object for the outcome model. |
| `M_formula` | Formula object for the mediator model. |
| `M_family` | Family for the mediator model passed to `glm` (e.g., `"gaussian"` or `binomial()`). |
| `L_formula_list` | A list of formulas for the exposure-induced confounder models. |
| `weights` | Optional numeric vector of weights used across models. |
| `boot` | Logical. If `TRUE`, performs a nonparametric bootstrap. |
| `boot_reps` | Number of bootstrap replications. |
| `boot_conf_level` | Confidence level for bootstrap intervals (default: 0.95). |
| `boot_seed` | Random seed for reproducibility. |
| `boot_parallel` | Logical. If `TRUE`, parallelizes the bootstrap. |
| `boot_cores` | Number of CPU cores to use. (default: all but two available cores) |

### Returns

- If `boot = FALSE`:  
  A list with:
  - `OE`, `IDE`, `IIE`, `CDE`: Estimated interventional and controlled effects
  - `models_L`: Fitted models for the exposure-induced confounders
  - `model_M`: Fitted mediator model
  - `model_Y`: Fitted outcome model
  - `data_ed`: Processed dataset with mean-centered covariates and residualized confounders

- If `boot = TRUE`, also includes:
  - `ci_*`: Bootstrap confidence intervals
  - `pvalue_*`: Bootstrap p-values from tests of no effect
  - `boot_*`: Vectors of bootstrap replicates for each effect

### Examples

#### Estimate Interventional Effects

```r
formula_L <- ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3
formula_M <- log_faminc_adj_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3
formula_Y <- std_cesd_age40 ~ att22 * log_faminc_adj_age3539 + ever_unemp_age3539 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

rwr_est <- rwrlite(
  data = nlsy,
  D = "att22",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0,
  m = log(5e4),
  Y_formula = formula_Y,
  M_formula = formula_M,
  L_formula_list = list(formula_L)
)
```

#### Bootstrap

```r
rwr_est_boot <- rwrlite(
  data = nlsy,
  D = "att22",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0,
  m = log(5e4),
  Y_formula = formula_Y,
  M_formula = formula_M,
  L_formula_list = list(formula_L),
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 60637
)
```

#### Parallelized Bootstrap

```r
rwr_est_bootpar <- rwrlite(
  data = nlsy,
  D = "att22",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  d = 1,
  dstar = 0,
  m = log(5e4),
  Y_formula = formula_Y,
  M_formula = formula_M,
  L_formula_list = list(formula_L),
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234,
  boot_parallel = TRUE
)
```


## `linpath`: analysis of path-specific effects using linear models

The `linpath` function estimates path-specific effects using a linear models for each mediator and the outcome.

### Function

```r
linpath(
  data,
  D,
  M,
  Y,
  C = NULL,
  d,
  dstar,
  interaction_DM = FALSE,
  interaction_DC = FALSE,
  interaction_MC = FALSE,
  weights_name = NULL,
  boot = FALSE,
  boot_reps = 1000,
  boot_conf_level = 0.95,
  boot_seed = NULL,
  boot_parallel = FALSE,
  boot_cores = NULL
)
```

### Arguments

| Argument         | Description |
|------------------|-------------|
| `data`           | A data frame. |
| `D`              | Name of the exposure variable (character). |
| `M`              | A character vector or list of vectors identifying mediators, in hypothesized causal order. See below for examples. |
| `Y`              | Name of the outcome variable (character). |
| `C`              | Optional character vector of baseline covariates to include in both the mediator and outcome models. |
| `d`, `dstar`     | Numeric values defining the exposure contrast of interest (`d - dstar`). |
| `interaction_DM` | Logical. If `TRUE`, includes exposure × mediator interactions in the outcome model. |
| `interaction_DC` | Logical. If `TRUE`, includes exposure × covariate interactions in both the mediator and outcome models. |
| `interaction_MC` | Logical. If `TRUE`, includes mediator × covariate interactions in the outcome model. |
| `weights_name`   | Optional name of sampling weights variable in `data`. |
| `boot` | Logical. If `TRUE`, use the nonparametric bootstrap to obtain confidence intervals and p-values. |
| `boot_reps` | Number of bootstrap replications (default: `1000`). |
| `boot_conf_level` | Confidence level for bootstrap interval (default: `0.95`). |
| `boot_seed` | Integer seed for reproducibility. |
| `boot_parallel` | Logical. If `TRUE`, parallelizes the bootstrap (requires `doParallel`, `doRNG`, and `foreach`). |
| `boot_cores` | Number of CPU cores for parallel bootstrap. If `NULL`, defaults to available cores minus 2. |


### Specifying Mediators with `M`

The `M` argument supports both univariate and multivariate mediators and encodes their assumed causal order.

- **Two ordered mediators:**

  ```r
  M = list("m1", "m2")
  ```

- **Partial ordering (treat two mediators as a block):**

  ```r
  M = list("m1", c("m2", "m3"))
  ```

- **Add a dummy-encoded factor mediator:**

  ```r
  M = list("m1", c("m2", "m3"), c("level2", "level3", "level4"))
  ```

The order of items within a group (e.g., `c("m2", "m3")`) doesn’t matter, but order of groups in the list does.

### Returns

By default, `linpath()` returns a list containing:

- `ATE`: Estimated average total effect for the contrast `d - dstar`.
- `PSE`: Named vector of path-specific effects.
- `miss_summary`: Counts of missing and non-missing values for `D`, `M`, `Y`, and `C`.
- `models_m`: Fitted mediator models.
- `models_y`: Outcome models, each with increasing numbers of mediators.

If `boot = TRUE`, it also returns:

- `ci_ATE`: Bootstrap confidence interval for ATE.
- `ci_PSE`: Bootstrap confidence intervals for PSEs.
- `pvalue_ATE`: P-value for test that the ATE is zero.
- `pvalue_PSE`: P-values for tests that the PSEs are zero.
- `boot_ATE`: ATE estimates from each bootstrap sample.
- `boot_PSE`: Matrix of PSE estimates from each bootstrap samples.

### Examples

#### Estimate Path-specific Effects through Two Mediators

```r
pse_est <- linpath(
  data = nlsy,
  D = "att22",
  M = list("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")
)
```

#### Add interactions

```r
pse_est_inter <- linpath(
  data = nlsy,
  D = "att22",
  M = list("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  interaction_DM = TRUE,
  interaction_DC = TRUE,
  interaction_MC = TRUE
)
```

#### Parallelized Bootstrap

```r
pse_est_boot <- linpath(
  data = nlsy,
  D = "att22",
  M = list("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 60637,
  boot_parallel = TRUE
)
```


## `ipwpath`: analysis of path-specific effects using inverse probability weights

The `ipwpath` function implements an inverse probability weighting approach to estimate path-specific effects.

### Function

```r
ipwpath(
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
  boot_cores = NULL
)
```

### Arguments

| Argument             | Description |
|----------------------|-------------|
| `data`               | A data frame. |
| `D`                  | Name of the exposure variable (character). |
| `M`                  | A character vector or list of vectors identifying mediators, in hypothesized causal order. |
| `Y`                  | Name of the outcome variable (character). |
| `C`                  | Optional vector of baseline covariates for the exposure models. |
| `base_weights_name`  | Optional name of a sampling weights variable. |
| `stabilize`          | Logical. If `TRUE`, stabilizes the weights using marginal probabilities of exposure (default: `TRUE`). |
| `censor` | Logical. If `TRUE`, applies weight censoring (default: `TRUE`). |
| `censor_low`, `censor_high` | Quantiles for censoring the weights (default: 0.01 and 0.99). |
| `boot` | Logical. If `TRUE`, use the nonparametric bootstrap to obtain confidence intervals and p-values. |
| `boot_reps` | Number of bootstrap replications (default: `1000`). |
| `boot_conf_level` | Confidence level for bootstrap interval (default: `0.95`). |
| `boot_seed` | Integer seed for reproducibility. |
| `boot_parallel` | Logical. If `TRUE`, parallelizes the bootstrap (requires `doParallel`, `doRNG`, and `foreach`). |
| `boot_cores` | Number of CPU cores for parallel bootstrap. If `NULL`, defaults to available cores minus 2. |

### Specifying Mediators with `M`

The `M` argument supports both univariate and multivariate mediators and encodes their assumed causal order.

- **Two ordered mediators:**

  ```r
  M = list("m1", "m2")
  ```

- **Partial ordering (treat two mediators as a block):**

  ```r
  M = list("m1", c("m2", "m3"))
  ```

- **Add a dummy-encoded factor mediator:**

  ```r
  M = list("m1", c("m2", "m3"), c("level2", "level3", "level4"))
  ```

The order of items within a group (e.g., `c("m2", "m3")`) doesn’t matter, but order of groups in the list does.

### Returns

By default, `ipwpath` returns:

- `ATE`: Estimated average total effect.
- `PSE`: Named vector of estimated path-specific effects.
- `weights1`, `weights2`: IP weights used for multivariate natural effect estimation.
- `weights3`: A matrix of IP weights for each set of mediators.
- `model_d1`: Fitted exposure model without mediators (`f(D|C)`).
- `models_d2`: List of exposure models with mediators (`s(D|C,M)`), one for each set of mediators.

If `boot = TRUE`, it also includes:

- `ci_ATE`: Bootstrap confidence interval for ATE.
- `ci_PSE`: Confidence intervals for PSEs.
- `pvalue_ATE`: P-value for test that ATE is zero.
- `pvalue_PSE`: P-values for tests that PSEs are zero.
- `boot_ATE`: ATE estimates from bootstrap samples.
- `boot_PSE`: PSE estimates from bootstrap samples.

### Examples

#### Estimate PSEs with two ordered mediators

```r
pse_est <- ipwpath(
  data = nlsy,
  D = "att22",
  M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")
)
```

#### Bootstrap

```r
pse_boot <- ipwpath(
  data = nlsy,
  D = "att22",
  M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234
)
```

#### Parallelized bootstrap

```r
pse_bootpar <- ipwpath(
  data = nlsy,
  D = "att22",
  M = c("ever_unemp_age3539", "log_faminc_adj_age3539"),
  Y = "std_cesd_age40",
  C = c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"),
  boot = TRUE,
  boot_reps = 2000,
  boot_seed = 1234,
  boot_parallel = TRUE
)
```


## `utils`: utility functions

This script defines helper functions used internally by many of the other other functions in this repository.

- **`demean(x, w)`**  
  Returns a weighted, mean-centered version of `x`. If no weights are provided, it defaults to equal weights.

- **`trim(x, min, max)`**  
  Top- and bottom-codes the vector `x` at fixed limits `min` and `max`. Useful for censoring inverse probability weights.

- **`trimQ(x, low, high)`**  
  Censors `x` at empirical quantiles defined by `low` and `high`. For example, `trimQ(x, 0.01, 0.99)` trims at the 1st and 99th percentiles. Useful for censoring inverse probability weights.

- **`comb_list_vec(...)`**  
  Combines multiple lists of vectors into a single list using `mapply()`.

These functions are are intended for internal use only.
