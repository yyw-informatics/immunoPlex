# Regenerate data/immunoplex_example.rda
#
# Builds the package's worked-example dataset on top of
# simulate_immunoassay(). Parameters are tuned against the v4 benchmark
# reports (B0, B5): a target ~25% censoring rate where Tobit/AFT clearly
# beat LOD/2 substitution but the analytes are not all wiped out.
#
# Run:
#   Rscript data-raw/immunoplex_example.R
#
# Produces data/immunoplex_example.rda via usethis::use_data().

suppressPackageStartupMessages({
  pkgload::load_all()           # picks up simulate_immunoassay() + analyte_library
  library(tidyr)
  library(dplyr)
  library(usethis)
})

set.seed(20260427)

# ---- Simulate ----------------------------------------------------------------

sim <- simulate_immunoassay(
  # Design ----
  n_subjects   = 120,                          # 60/group, well above 25/cond
  n_timepoints = 2,                            # T1 / T2
  n_analytes   = 20,                           # full built-in analyte library
  design       = "pre_post",
  group_levels = c("No", "Yes"),               # disease = No (reference) vs Yes
  group_allocation = c(1, 1),

  # Biology ----
  signal_analytes     = c("IL-6", "IL-8", "TNF-a", "IFN-g",
                          "MCP-1", "IP-10"),  # mix of pro-inflam / Th1 / chemokine
  group_effects       = 1.2,                   # log-scale effect, visible separation
  effect_direction    = "up",                  # disease=Yes elevated
  time_effects        = 0.3,                   # mild drift T1 -> T2
  interaction_effects = 0.4,                   # disease effect grows at T2
  re_intercept_sd     = 0.5,
  residual_sd         = 0.8,
  analyte_correlation = "block",               # Th1 / Th2 / chemokine clusters
  block_rho           = 0.5,

  # Subject covariates (effects on log-conc) ----
  covariates = list(
    age      = list(type = "continuous", effect = 0.15),
    preexist = list(type = "binary",     effect = 0.30)
  ),

  # Assay layer ----
  lod_quantile = 0.25,                         # ~25% per-analyte censoring
  ulod_values  = c("IL-8" = 8000,              # high-abundance analytes saturate
                   "RANTES" = 12000,
                   "IP-10" = 8000,
                   "MIP-1b" = 6000),
  n_replicates = 2,
  replicate_sd = 0.10,
  n_plates     = 1,
  n_lots       = 1,                            # single LOD per analyte

  seed = 20260427
)

cat("Realized per-analyte censoring (LOD):\n")
print(round(sort(sim$meta$realized_censoring, decreasing = TRUE), 3))
cat(sprintf("\nOverall LOD censoring: %.1f%%\n",
            100 * mean(sim$data$cens_lod, na.rm = TRUE)))
cat(sprintf("Overall ULOD censoring: %.1f%%\n",
            100 * mean(sim$data$cens_ulod, na.rm = TRUE)))


# ---- Reshape into expression / metadata / lod_lookup -------------------------

d <- sim$data

# Realistic ages: rescale the standardized covariate to integer years.
# (The latent values were already generated with the z-scaled covariate, so
# this rescale only affects how `age` is reported in metadata.)
subj_age <- d %>%
  distinct(subject_id, age) %>%
  mutate(age = pmin(80L, pmax(20L, as.integer(round(50 + 15 * age)))))

subj_pre <- d %>%
  distinct(subject_id, preexist) %>%
  mutate(preexist = ifelse(preexist == 1L, "Yes", "No"))

samples <- d %>%
  distinct(subject_id, timepoint, replicate_id, group) %>%
  arrange(subject_id, timepoint, replicate_id) %>%
  mutate(
    sample_id = sprintf("%s_%s_R%d",
                        as.character(subject_id),
                        as.character(timepoint),
                        replicate_id)
  ) %>%
  left_join(subj_age, by = "subject_id") %>%
  left_join(subj_pre, by = "subject_id") %>%
  rename(disease = group, replicate = replicate_id) %>%
  mutate(
    subject_id = as.character(subject_id),
    timepoint  = as.character(timepoint),
    disease    = as.character(disease)
  )

metadata <- samples[, c("sample_id", "subject_id", "timepoint",
                        "disease", "age", "preexist", "replicate")]

# Wide expression matrix on the raw (linear) scale -- what immuno_preprocess
# expects as input.
#
# The simulator clamps censored raw values *at* the LOD, but downstream
# `prepare_cytokine_data()` flags censoring with strict `value < lod`.
# For censored rows, swap the clamped value for a sub-LOD reading drawn
# from the latent truth (capped at lod * 0.999) -- a plausible "instrument
# reads just below detection" stand-in that ensures the strict-less-than
# flag fires.
d_export <- d
cens_idx <- which(d_export$cens_lod)
if (length(cens_idx) > 0) {
  truth_below <- sim$truth$true_raw[cens_idx]
  d_export$value_raw[cens_idx] <- pmin(truth_below,
                                       d_export$lod[cens_idx] * 0.999)
}

expr_long <- d_export %>%
  mutate(sample_id = sprintf("%s_%s_R%d",
                             as.character(subject_id),
                             as.character(timepoint),
                             replicate_id)) %>%
  select(sample_id, cytokine, value_raw)

expression <- pivot_wider(expr_long,
                          names_from = cytokine,
                          values_from = value_raw)
expression <- as.data.frame(expression)
rownames(expression) <- expression$sample_id
expression$sample_id <- NULL
expression <- expression[metadata$sample_id, , drop = FALSE]

# LOD lookup: single LOD per analyte (n_lots = 1, so all rows share it).
analyte_names <- colnames(expression)
lod_vec  <- sim$meta$dgp_params$assay$lod_values[analyte_names]
ulod_vec <- sim$meta$dgp_params$assay$ulod_values[analyte_names]
ulod_vec[!is.finite(ulod_vec)] <- NA_real_

lod_lookup <- data.frame(
  cytokine = analyte_names,
  lod      = unname(lod_vec),
  ulod     = unname(ulod_vec),
  stringsAsFactors = FALSE
)


# ---- Sanity checks -----------------------------------------------------------

stopifnot(
  identical(rownames(expression), metadata$sample_id),
  all(metadata$timepoint %in% c("T1", "T2")),
  all(metadata$disease   %in% c("Yes", "No")),
  all(metadata$preexist  %in% c("Yes", "No")),
  nrow(metadata) == nrow(expression),
  ncol(expression) == 20,
  nrow(lod_lookup) == 20
)

# Confirm per-condition n
n_by_cond <- metadata %>%
  distinct(subject_id, disease) %>%
  count(disease)
cat("\nSubjects per disease group:\n")
print(n_by_cond)
stopifnot(all(n_by_cond$n >= 25))


# ---- Assemble and save -------------------------------------------------------

immunoplex_example <- list(
  expression = expression,
  metadata   = metadata,
  lod_lookup = lod_lookup
)

use_data(immunoplex_example, overwrite = TRUE, compress = "xz")

cat("\nSaved data/immunoplex_example.rda\n")
cat(sprintf("  expression: %d samples x %d cytokines\n",
            nrow(expression), ncol(expression)))
cat(sprintf("  metadata:   %d rows, cols = %s\n",
            nrow(metadata), paste(colnames(metadata), collapse = ", ")))
cat("  cytokines:  ", paste(analyte_names, collapse = ", "), "\n")
