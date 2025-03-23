# causalMedR
This repository contains R functions for conducting causal mediation analysis using the methods described in Wodtke and Zhou "Causal Mediation Analysis."


## Table of Contents
- [linmed – mediation analysis using linear models](#linmed-mediation-analysis-using-linear-models)
- [medsim – mediation analysis using a simulation approach](#medsim-causal-mediation-analysis-using-a-simulation-approach)


## `linmed`: mediation analysis using linear models

The `linmed` function estimates **natural direct**, **natural indirect**, **controlled direct**, and **total effects** using linear models. It supports mediation analysis with a single mediator and with multiple mediators, and includes built-in support for **nonparametric bootstrapping** with optional **parallel computation**.

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

---

### Examples

#### Single Mediator, No Interactions

```r
linmed(
  data = nlsy1,
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
  data = nlsy1,
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
  data = nlsy1,
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
  data = nlsy2,
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
  data = nlsy1,
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
  data = nlsy1,
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


## `medsim`: causal mediation analysis using a simulation approach

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
| `cat_list`    | Character vector of treatment levels to compare (default: `c(0, 1)`). |
| `treatment`   | Name of the treatment variable (string). |
| `intv_med`    | A list specifying the intervention on the mediators. Set to `NULL` if interventional or controlled direct effects are not of interest. |
| `model_spec`  | A list of lists defining the models for the mediators and outcome. Each model must include: <br> • `func`: model-fitting function (e.g., `"glm"`, `"polr"`) <br> • `formula`: model formula <br> • `args`: (optional) list of additional arguments to pass to the function. |
| `weights`     | (Optional) Name of the variable containing weights to use in model fitting. If `NULL`, no weights are applied. |
| `boot`        | Logical. If `TRUE`, performs nonparametric bootstrap to obtain confidence intervals and p-values (default: `FALSE`). Requires `doParallel`, `doRNG`, and `foreach`. |
| `reps`        | Integer. Number of bootstrap replications (default: `100`). |
| `resv_core`   | Integer. Number of CPU cores to reserve (i.e., not use) when bootstrapping via parallel processing (default: `1`). |
| `seed`        | Integer or `NULL`. Seed for reproducibility. |

### Returns

- If `boot = FALSE`:  
  A **list** of point estimates for the mediation effects of interest.

- If `boot = TRUE`:  
  A **data frame** containing:
  - Point estimates  
  - 95% bootstrap confidence intervals  
  - P-values from test of no effect
