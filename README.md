# causalMedR
This repository contains R functions for conducting causal mediation analysis using the methods outlined in Wodtke and Zhou "Causal Mediation Analysis."


## Table of Contents
- [linmed – mediation analysis using linear models](#linmed--mediation-analysis-using-linear-models)


## `linmed()`: mediation analysis using linear models

The `linmed()` function estimates **natural direct**, **natural indirect**, **controlled direct**, and **total effects** using linear models. It supports mediation analysis with a single mediator and with multiple mediators, and includes built-in support for **nonparametric bootstrapping** with optional **parallel computation**.

### Function Signature

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
