
# immunoPlex

**Version:** 0.0.0.9000 (Development)

## Overview

`immunoPlex` is an R package for statistical analysis of multiplex cytokine assays with limit-of-detection (LOD) censoring. The package implements censoring-aware regression models, automated model comparison with estimand tracking, and multivariate pattern discovery methods suitable for immunological and biomarker research.

## Features

### Statistical Methods
- Gaussian generalized linear mixed models (GLMM) for uncensored continuous outcomes
- Gamma generalized linear mixed models (GLMM) with log link
- Tobit regression for left- and right-censored data
- Accelerated failure time (AFT) models with log-normal distribution
- Proper interval coding for censored observations
- Homogeneous vs. heterogeneous variance model comparison (LRT)
- Estimand specification (ratio of means, ratio of medians, mean difference)
- Model selection via Akaike Information Criterion (AIC)

### Multivariate Analysis
- Partial Least Squares Discriminant Analysis (PLS-DA)
- Orthogonal PLS-DA (OPLS-DA) for enhanced interpretation
- Permutation testing for model validation
- Variable importance in projection (VIP) scores

### Detection Analysis
- McNemar's exact test for paired detection frequency changes
- Visualization of detection patterns across conditions

### Preprocessing
- Conditional covariate selection based on prevalence thresholds
- IQR-based winsorization within groups for outlier management

### Diagnostics
- Censoring-aware residual plots
- DHARMa simulation-based residual diagnostics (uniformity, dispersion, outlier tests)
- Convergence assessment (Hessian positive-definiteness, gradient magnitude, NA SEs)
- Random effects diagnostics (ICC, variance components, Shapiro-Wilk normality)
- Leave-one-out sensitivity analysis for influential subject detection
- Quantile-quantile plots with randomized quantile residuals

## Installation

### Development Version

Install from source:

``` r
# Install devtools if not already available
install.packages("devtools")

# Install immunoPlex from local source or GitLab
devtools::install()
```

### Dependencies

**Required packages:**
- `glmmTMB` (≥ 1.0.0) - Generalized linear mixed models
- `survival` - Survival analysis and censored regression
- `ggplot2` - Graphics
- `patchwork` - Plot composition
- `dplyr`, `tidyr`, `tibble` - Data manipulation

**Optional packages:**
- `DHARMa` - Residual diagnostics for hierarchical models
- `broom.mixed` - Tidy extraction of mixed model coefficients
- `statmod` - Randomized quantile residuals
- `ropls` (Bioconductor) - Required for PLS-DA functionality
- `censReg` - Alternative Tobit implementation

## Usage Example

### Basic Workflow

``` r
library(immunoPlex)

# Load example dataset
data("immunoplex_example", package = "immunoPlex")

# Preprocess cytokine expression data
dat <- immuno_preprocess(
  expression = immunoplex_example$expression,
  metadata   = immunoplex_example$metadata,
  lod_lookup = immunoplex_example$lod_lookup
)

# Fit single model
fit <- fit_one(dat, family = "gamma", random = "")
summary(fit)
plot(fit, plot_type = "censor_aware", save_pdf = FALSE)

# Compare multiple model families
models <- fit_models(dat, families = c("gamma", "tobit", "aft"), 
                     ulod = TRUE, random = "")
print(models)
plot(models, plot_type = "comparison", save_pdf = FALSE)

# Examine model comparison results
models$comparison[, c("family", "estimand", "aic", "delta_aic")]
#>   family         estimand    aic delta_aic
#>    gamma   ratio_of_means  245.67      0.00
#>    tobit   ratio_of_means  248.32      2.65  
#>      aft ratio_of_medians  251.45      5.78
```

## Data Structure Requirements

### Using `immuno_preprocess()` (Recommended)

The `immuno_preprocess()` function formats data from wide-format expression matrices into the required structure for modeling.

### Manual Data Preparation

If preparing data manually, the data frame must contain:

**Required columns:**
- `value` - Numeric cytokine concentration
- `cens_lod` - Logical indicator for left-censoring at LOD
- `lod` - Numeric LOD threshold (when `cens_lod` is TRUE)

**Optional columns:**
- `cens_ulod` - Logical indicator for right-censoring at upper LOD
- `ulod` - Numeric upper LOD threshold (when `cens_ulod` is TRUE)
- `subject_id` - Subject identifier for random effects models
- Model covariates (default: `timepoint`, `disease`, `age`)

**Note:** The function `prepare_cytokine_data()` provides a convenience wrapper for the included example dataset.

## Statistical Models

### Available Model Families

All models report the specific estimand being estimated:

| Family | Implementation | Link Function | Estimand | Censoring Handling |
|--------|---------------|---------------|----------|-------------------|
| `gaussian` | `glmmTMB` | identity | Mean difference | No censoring (continuous) |
| `gamma` | `glmmTMB` | log | Ratio of means | Imputation at LOD |
| `tobit` | `survreg` | identity | Ratio of means | Interval coding `[-∞, LOD]` |
| `aft` | `survreg` | log | Ratio of medians | Interval coding `[ε, LOD]` |
| `tobit_censreg` | `censReg` | identity | Ratio of means | Single LOD only |
| `auto` | Heuristic | varies | varies | Data-driven selection |

### Model Specification

**Default fixed effects:** `timepoint*disease + age`

**Random effects:** `(1|subject_id)` for GLMM families (gaussian, gamma) when:
- Subject-level grouping variable exists
- Sample size exceeds minimum threshold (configurable)
- Adequate within-subject replication

**Dispersion modeling:** Gaussian models support a `dispformula` argument for modeling heterogeneous variance (e.g., `dispformula = ~timepoint`).

**Optimizer cascade:** Gaussian models use a multi-optimizer cascade (nlminb → BFGS → L-BFGS-B) with automatic fallback to fixed-effects-only if random effects cause convergence failures.

**Model comparison:** Akaike Information Criterion (AIC) with delta-AIC reported.

**Visualization:** The `plot()` method for `immuno_model_set` objects displays model comparison tables and diagnostic plots for the best-fitting model.

## Preprocessing Helpers

``` r
# Conditional covariate selection based on prevalence
covars <- choose_covariates(
  data = dat,
  base_covars = c("timepoint", "disease", "age"),
  candidates = c("HIV", "malaria"),
  rare_threshold = 0.03  # exclude if < 3% prevalence
)

# IQR-based winsorization within groups
dat_clean <- winsorize_by_group(
  data = dat,
  value_col = "value",
  group_cols = c("cytokine", "timepoint"),
  iqr_mult = 3
)
```

## Model Diagnostics

### Available Diagnostic Plots

**Residual plots:**
- Residuals vs. fitted values
- Quantile-quantile plots
- Censoring-aware visualizations (censored observations highlighted)

**Residual types:**
- Randomized quantile residuals (Dunn & Smyth, 1996) for censored models
- Pearson residuals for Gamma GLMM
- DHARMa scaled residuals when `DHARMa` package available

**Usage:**

``` r
# Diagnostic plots for best model
plot(models, plot_type = "best_residuals", save_pdf = FALSE)

# Diagnostic plots for all fitted models
plot(models, plot_type = "all_residuals", save_pdf = FALSE)

# DHARMa diagnostics
plot(fit, dharma = TRUE)
```

### glmmTMB Model Diagnostics

For Gaussian and Gamma GLMM models fitted via `glmmTMB`, dedicated diagnostic functions provide comprehensive model checking:

``` r
fit <- fit_one(dat, family = "gaussian",
               fixed = "timepoint + disease + age",
               random = "(1|subject_id)")

# DHARMa simulation-based diagnostics
diag <- dharma_diagnostics(fit, n_sim = 500)
diag$passed_all
diag$tests  # uniformity, dispersion, outlier, quantile deviation

# Convergence assessment
conv <- check_convergence(fit)
conv$converged      # overall verdict
conv$hessian_pd     # Hessian positive-definite?
conv$max_gradient   # largest absolute gradient component
conv$na_se_count    # number of NA standard errors

# Random effects diagnostics
re <- check_random_effects(fit)
re$icc              # intra-class correlation
re$re_variance      # random effect variance
re$re_normality_p   # Shapiro-Wilk p-value for RE normality
```

### Variance Model Comparison

Compare homogeneous vs. heterogeneous variance structures using a likelihood ratio test:

``` r
cmp <- compare_variance_models(
  data = dat,
  fixed = "timepoint + disease + age",
  random = "(1|subject_id)",
  dispformula_hetero = ~timepoint
)

cmp$preferred_model    # "homogeneous" or "heterogeneous"
cmp$preference_reason  # explanation (LRT p-value, AIC difference)
cmp$comparison         # tibble with AIC, BIC, LRT statistic and p-value
```

### Leave-One-Out Sensitivity

Identify influential subjects whose removal substantially changes coefficient estimates:

``` r
loo <- loo_sensitivity(
  data = dat,
  fixed = "timepoint + disease + age",
  target_terms = c("timepointT2", "timepointT3"),
  random = "(1|subject_id)",
  influence_threshold = 20  # flag if |% change| > 20%
)

# Results: tibble with excluded_subject, term, full_estimate,
#          loo_estimate, difference, pct_change, influential
loo[loo$influential, ]
```

## Configuration Options

Global configuration parameters can be set using standard R `options()`:

``` r
options(
  fit_one.min_subjects = 30L,  # Minimum subjects to retain random effects
  fit_one.min_reps     = 3L,   # Minimum replicates per subject for random effects

  # glmmTMB convergence thresholds
  immunoPlex.glmmTMB_gradient_warn = 0.01,   # Warning threshold for max gradient
  immunoPlex.glmmTMB_gradient_fail = 0.1,    # Failure threshold for max gradient

  # LOO sensitivity defaults
  immunoPlex.glmmTMB_loo_influence_threshold = 20  # % change to flag as influential
)
```

Use `get_immunoplex_config()` and `set_immunoplex_config()` to inspect and modify these settings programmatically.

## Comparing LOD Handling Approaches

The function `compare_lod_models()` compares preprocessing-based imputation methods (e.g., LOD/2, LOD/√2) with censoring-aware regression models:

``` r
cmp <- compare_lod_models(dat)
print(cmp)
plot(cmp, plot_best = TRUE)
```

This comparison helps evaluate the impact of different censoring handling strategies on inference.

## Multivariate Pattern Discovery

### PLS-DA and OPLS-DA

The package implements supervised multivariate methods for identifying cytokine expression patterns that discriminate between experimental groups.

**Methods available:**
- Partial Least Squares Discriminant Analysis (PLS-DA)
- Orthogonal PLS-DA (OPLS-DA) for enhanced interpretability

**Model validation:**
- Permutation testing (e.g., 1000 permutations)
- Cross-validation metrics (R²X, R²Y, Q²)
- Statistical significance testing for model performance

**Example workflow:**

``` r
# Data preprocessing
preprocessed <- plsda_preprocess(
  data = your_data,
  cytokine_cols = cytokines,
  metadata_cols = metadata,
  lod_lookup = lod_values,
  lod_method = "half"
)

# Model fitting with permutation testing
model <- plsda_fit(
  preprocessed_data = preprocessed,
  response_var = "group",
  n_components = 2,
  method = "PLS-DA",        # or "OPLS-DA"
  permutations = 1000
)

# Visualization
plot(model, plot_type = "scores", color_by = "group")
plot(model, plot_type = "vip", top_n = 20)
plot(model, plot_type = "biplot", color_by = "group")
plot(model, plot_type = "loadings")
```

**Variable importance:** VIP (Variable Importance in Projection) scores identify cytokines contributing most to group discrimination.

See `vignette("plsda-analysis")` for detailed methodology and interpretation.

## Example Dataset

The package includes a synthetic dataset for demonstration:

``` r
data("immunoplex_example", package = "immunoPlex")
```

**Components:**
- `expression` - Cytokine concentration matrix (samples × cytokines)
- `metadata` - Sample annotations (subject ID, timepoint, disease status, age)
- `lod_lookup` - Detection limits for each cytokine (LOD and ULOD)

## Documentation

**Function documentation:**
``` r
# Core modeling
?fit_one
?fit_models
?immuno_preprocess
?compare_lod_models

# Preprocessing helpers
?choose_covariates
?winsorize_by_group

# glmmTMB diagnostics
?dharma_diagnostics
?check_convergence
?check_random_effects
?compare_variance_models
?loo_sensitivity

# Multivariate
?plsda_preprocess
?plsda_fit
```

**Vignettes:**
``` r
vignette("mcnemar-detection-analysis")
vignette("plsda-analysis")
```

## Citation

To cite `immunoPlex` in publications, use:

``` r
citation("immunoPlex")
```

## License

This package is licensed under the MIT License. See `LICENSE` file for details.

## References

- Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., Skaug, H. J., Mächler, M., & Bolker, B. M. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. *The R Journal*, 9(2), 378-400. doi:10.32614/RJ-2017-066

- Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals. *Journal of Computational and Graphical Statistics*, 5(3), 236-244. doi:10.1080/10618600.1996.10474708

- Fay, M. P., & Lumbard, K. (2021). Confidence intervals for difference in proportions for matched pairs compatible with exact McNemar's or sign tests. *Statistics in Medicine*, 40(5), 1147-1159. doi:10.1002/sim.8829

- Hartig, F. (2022). DHARMa: Residual diagnostics for hierarchical (multi-level/mixed) regression models. R package version 0.4.6. https://CRAN.R-project.org/package=DHARMa

- Pedersen, T. L. (2024). *patchwork: The Composer of Plots*. R package version 1.2.0. https://CRAN.R-project.org/package=patchwork

- Therneau, T. M. (2024). *A Package for Survival Analysis in R*. R package version 3.5-8. https://CRAN.R-project.org/package=survival

- Thévenot, E. A., Roux, A., Xu, Y., Ezan, E., & Junot, C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index, and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, 14(8), 3322-3335. doi:10.1021/acs.jproteome.5b00354

- Tobin, J. (1958). Estimation of relationships for limited dependent variables. *Econometrica*, 26(1), 24-36. doi:10.2307/1907382

- Trygg, J., & Wold, S. (2002). Orthogonal projections to latent structures (O-PLS). *Journal of Chemometrics*, 16(3), 119-128. doi:10.1002/cem.695

- Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. https://ggplot2.tidyverse.org

- Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics. *Chemometrics and Intelligent Laboratory Systems*, 58(2), 109-130. doi:10.1016/S0169-7439(01)00155-1
