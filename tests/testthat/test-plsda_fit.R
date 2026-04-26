# Tests for plsda_fit function

# Helper function to create test data
create_test_plsda_data <- function(n_samples = 30, n_cytokines = 10) {
  test_data <- data.frame(
    SubjectID = paste0("S", 1:n_samples),
    group = rep(c("Control", "Treatment"), each = n_samples/2),
    site = rep(c("A", "B", "C"), length.out = n_samples),
    stringsAsFactors = FALSE
  )
  
  # Add cytokines
  for (i in 1:n_cytokines) {
    # Create some group separation in the data
    control_vals <- rnorm(n_samples/2, mean = 100, sd = 20)
    treatment_vals <- rnorm(n_samples/2, mean = 120, sd = 20)
    test_data[[paste0("cyt", i)]] <- c(control_vals, treatment_vals)
  }
  
  lod_df <- data.frame(
    cytokine = paste0("cyt", 1:n_cytokines),
    lod = rep(10, n_cytokines)
  )
  
  list(data = test_data, lod = lod_df)
}


test_that("plsda_fit requires ropls package", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  # Should work if ropls is available
  expect_no_error(
    plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "group",
      n_components = 2,
      permutations = 0,  # Skip for speed
      verbose = FALSE
    )
  )
})


test_that("plsda_fit validates inputs correctly", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  # Should error with wrong input type
  expect_error(
    plsda_fit(
      preprocessed_data = list(a = 1, b = 2),
      response_var = "group",
      verbose = FALSE
    ),
    "must be a plsda_preprocessed object"
  )
  
  # Should error with missing response variable
  expect_error(
    plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "nonexistent_var",
      verbose = FALSE
    ),
    "not found in metadata"
  )
})


test_that("plsda_fit returns correct structure", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    permutations = 0,
    verbose = FALSE
  )
  
  # Check class
  expect_s3_class(model, "plsda_model")
  expect_type(model, "list")
  
  # Check components
  expect_true("model" %in% names(model))
  expect_true("scores" %in% names(model))
  expect_true("loadings" %in% names(model))
  expect_true("vip_scores" %in% names(model))
  expect_true("model_stats" %in% names(model))
  expect_true("response_var" %in% names(model))
  
  # Check scores structure
  expect_true(is.data.frame(model$scores))
  expect_true("p1" %in% names(model$scores))
  expect_true("p2" %in% names(model$scores))
  expect_true("group" %in% names(model$scores))
  
  # Check loadings structure
  expect_true(is.data.frame(model$loadings))
  expect_true("cytokine" %in% names(model$loadings))
  
  # Check VIP scores structure
  expect_true(is.data.frame(model$vip_scores))
  expect_true("cytokine" %in% names(model$vip_scores))
  expect_true("vip_score" %in% names(model$vip_scores))
  expect_true("importance" %in% names(model$vip_scores))
  
  # Check model stats
  expect_true(is.data.frame(model$model_stats))
  expect_true("R2X" %in% names(model$model_stats))
  expect_true("R2Y" %in% names(model$model_stats))
  expect_true("Q2" %in% names(model$model_stats))
})


test_that("plsda_fit handles OPLS-DA correctly", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  expect_equal(model$method, "OPLS-DA")
  
  # OPLS-DA should have orthogonal components
  expect_true("o1" %in% names(model$scores))
})


test_that("plsda_fit print and summary methods work", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    permutations = 0,
    verbose = FALSE
  )
  
  # Print should not error
  expect_output(print(model), "PLS-DA Model")
  expect_output(print(model), "Method: PLS-DA")
  
  # Summary should not error
  expect_output(summary(model), "PLS-DA Model Summary")
  expect_output(summary(model), "Top 10 Cytokines")
})


# ============================================================================
# Tests for Cross-Validation Functionality
# ============================================================================

test_that("recommend_cv_method returns correct structure", {
  rec <- recommend_cv_method(n_total = 50, n_min_group = 15, verbose = FALSE)
  
  # Check return structure
  expect_type(rec, "list")
  expect_true("method" %in% names(rec))
  expect_true("description" %in% names(rec))
  expect_true("rationale" %in% names(rec))
  expect_true("considerations" %in% names(rec))
  
  # Check field types
  expect_true(is.character(rec$method) || is.numeric(rec$method))
  expect_type(rec$description, "character")
  expect_type(rec$rationale, "character")
  expect_type(rec$considerations, "character")
})


test_that("recommend_cv_method recommends LOOCV for small samples", {
  # Small total sample size
  rec1 <- recommend_cv_method(n_total = 40, n_min_group = 15, verbose = FALSE)
  expect_equal(rec1$method, "loocv")
  expect_match(rec1$description, "LOOCV")
  
  # Small group size (overrides large total)
  rec2 <- recommend_cv_method(n_total = 80, n_min_group = 12, verbose = FALSE)
  expect_equal(rec2$method, "loocv")
  
  # Edge case: exactly 50 samples but small group
  rec3 <- recommend_cv_method(n_total = 50, n_min_group = 14, verbose = FALSE)
  expect_equal(rec3$method, "loocv")
})


test_that("recommend_cv_method recommends k-fold for larger samples", {
  # Moderate sample size
  rec1 <- recommend_cv_method(n_total = 75, n_min_group = 20, verbose = FALSE)
  expect_equal(rec1$method, 7)
  expect_match(rec1$description, "7-fold")
  
  # Large sample size
  rec2 <- recommend_cv_method(n_total = 150, n_min_group = 40, verbose = FALSE)
  expect_equal(rec2$method, 10)
  expect_match(rec2$description, "10-fold")
  
  # Edge case: exactly 100 samples with adequate groups
  rec3 <- recommend_cv_method(n_total = 100, n_min_group = 20, verbose = FALSE)
  expect_equal(rec3$method, 10)
})


test_that("recommend_cv_method verbose mode works", {
  # Should print detailed output
  expect_output(
    recommend_cv_method(n_total = 47, n_min_group = 15, verbose = TRUE),
    "Cross-Validation Method Recommendation"
  )
  expect_output(
    recommend_cv_method(n_total = 47, n_min_group = 15, verbose = TRUE),
    "Recommended:"
  )
  expect_output(
    recommend_cv_method(n_total = 47, n_min_group = 15, verbose = TRUE),
    "Considerations:"
  )
})


test_that("plsda_fit accepts different cross_validation parameters", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data(n_samples = 30)
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  # Test NULL (auto-select)
  expect_no_error(
    model1 <- plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "group",
      n_components = 2,
      cross_validation = NULL,
      permutations = 0,
      verbose = FALSE
    )
  )
  
  # Test LOOCV string
  expect_no_error(
    model2 <- plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "group",
      n_components = 2,
      cross_validation = "loocv",
      permutations = 0,
      verbose = FALSE
    )
  )
  
  # Test "loo" alias
  expect_no_error(
    model3 <- plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "group",
      n_components = 2,
      cross_validation = "loo",
      permutations = 0,
      verbose = FALSE
    )
  )
  
  # Test integer k-fold
  expect_no_error(
    model4 <- plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "group",
      n_components = 2,
      cross_validation = 7,
      permutations = 0,
      verbose = FALSE
    )
  )
})


test_that("plsda_fit validates cross_validation parameter", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  # Invalid string
  expect_error(
    plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "group",
      cross_validation = "invalid",
      verbose = FALSE
    ),
    "must be an integer.*or 'loo'/'loocv'"
  )
  
  # Too small integer
  expect_error(
    plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "group",
      cross_validation = 2,
      verbose = FALSE
    ),
    "must be >= 3"
  )
})


test_that("plsda_fit auto-selects CV correctly", {
  skip_if_not_installed("ropls")
  
  # Small sample: should select LOOCV
  test_obj_small <- create_test_plsda_data(n_samples = 40)
  preprocessed_small <- plsda_preprocess(
    data = test_obj_small$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj_small$lod,
    verbose = FALSE
  )
  
  model_small <- plsda_fit(
    preprocessed_data = preprocessed_small,
    response_var = "group",
    n_components = 2,
    cross_validation = NULL,  # Auto-select
    permutations = 0,
    verbose = FALSE
  )
  
  # Should have selected LOOCV (n = total samples for LOOCV)
  expect_equal(model_small$cv_folds, 40)
  expect_true(model_small$cv_method == "loocv" || model_small$cv_folds == 40)
})


test_that("plsda_fit stores CV information in result", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    cross_validation = "loocv",
    permutations = 0,
    verbose = FALSE
  )
  
  # Check CV info is stored
  expect_true("cv_method" %in% names(model))
  expect_true("cv_folds" %in% names(model))
  expect_equal(model$cv_method, "loocv")
  expect_equal(model$cv_folds, 30)  # Total samples
})


test_that("plsda_fit print shows CV method", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_test_plsda_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  # Test LOOCV display
  model_loocv <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    cross_validation = "loocv",
    permutations = 0,
    verbose = FALSE
  )
  
  expect_output(print(model_loocv), "Cross-validation.*LOOCV")
  
  # Test k-fold display
  model_kfold <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    cross_validation = 7,
    permutations = 0,
    verbose = FALSE
  )
  
  expect_output(print(model_kfold), "Cross-validation.*7 -fold")
})


# ---- rep_col / replicate_agg (pre-aggregation path) -------------------------

test_that("plsda_fit with rep_col matches pre-aggregated data with rep_col = NULL", {
  skip_if_not_installed("ropls")

  set.seed(20260424)
  n_subj <- 30L
  n_cyt  <- 8L
  n_rep  <- 2L

  subj_ids <- paste0("S", seq_len(n_subj))
  grp      <- rep(c("Control", "Treatment"), length.out = n_subj)

  # Per-subject latent means shared across replicates.
  latent <- matrix(rnorm(n_subj * n_cyt, mean = 100, sd = 20),
                   nrow = n_subj, ncol = n_cyt)
  latent[grp == "Treatment", 1:(n_cyt %/% 2)] <-
    latent[grp == "Treatment", 1:(n_cyt %/% 2)] + 25

  rep_df <- data.frame(
    SubjectID = rep(subj_ids, each = n_rep),
    group     = rep(grp,      each = n_rep),
    rep_id    = rep(seq_len(n_rep), times = n_subj),
    stringsAsFactors = FALSE
  )
  # Replicate-level noise on top of the latent matrix.
  cyto_mat <- matrix(NA_real_, nrow = n_subj * n_rep, ncol = n_cyt,
                     dimnames = list(NULL, paste0("cyt", seq_len(n_cyt))))
  for (i in seq_len(n_subj)) {
    for (r in seq_len(n_rep)) {
      row <- (i - 1L) * n_rep + r
      cyto_mat[row, ] <- latent[i, ] + rnorm(n_cyt, 0, 3)
    }
  }
  rep_df <- cbind(rep_df, as.data.frame(cyto_mat))

  lod_df <- data.frame(cytokine = paste0("cyt", seq_len(n_cyt)),
                       lod      = rep(5, n_cyt))

  # Preprocess with replicates kept (replicate_col = NULL in plsda_preprocess).
  pre_reps <- plsda_preprocess(
    data          = rep_df,
    cytokine_cols = paste0("cyt", seq_len(n_cyt)),
    metadata_cols = c("SubjectID", "group", "rep_id"),
    lod_lookup    = lod_df,
    scale_data    = TRUE,
    verbose       = FALSE
  )
  expect_equal(nrow(pre_reps$expression), 2L * n_subj)

  # --- Path A: plsda_fit aggregates internally via rep_col ---
  mod_A <- plsda_fit(
    preprocessed_data = pre_reps,
    response_var      = "group",
    n_components      = 2,
    permutations      = 0,
    rep_col           = "rep_id",
    replicate_agg     = "mean",
    verbose           = FALSE
  )

  # --- Path B: aggregate_replicates on preprocessed expression, then fit ---
  expr <- pre_reps$expression
  meta <- pre_reps$metadata
  long_df <- data.frame(
    .grp_key = rep(paste(meta$SubjectID, meta$group, sep = "|"),
                   times = ncol(expr)),
    .cyto    = rep(colnames(expr), each = nrow(expr)),
    .rep     = rep(meta$rep_id, times = ncol(expr)),
    .value   = as.vector(expr),
    stringsAsFactors = FALSE
  )
  agg_long <- aggregate_replicates(
    data = long_df, rep_col = ".rep",
    subject_col = ".grp_key", cytokine_col = ".cyto",
    timepoint_col = NULL, value_col = ".value", rule = "mean"
  )
  unique_keys <- unique(long_df$.grp_key)
  expr_agg <- matrix(NA_real_,
                     nrow = length(unique_keys), ncol = ncol(expr),
                     dimnames = list(NULL, colnames(expr)))
  row_idx <- match(agg_long$.grp_key, unique_keys)
  col_idx <- match(agg_long$.cyto,    colnames(expr))
  expr_agg[cbind(row_idx, col_idx)] <- agg_long$.value

  first_idx <- match(unique_keys, paste(meta$SubjectID, meta$group, sep = "|"))
  meta_agg <- meta[first_idx, c("SubjectID", "group"), drop = FALSE]
  rownames(meta_agg) <- NULL

  pre_agg <- pre_reps
  pre_agg$expression <- expr_agg
  pre_agg$metadata   <- meta_agg
  pre_agg$sample_ids <- pre_reps$sample_ids[first_idx]

  mod_B <- plsda_fit(
    preprocessed_data = pre_agg,
    response_var      = "group",
    n_components      = 2,
    permutations      = 0,
    rep_col           = NULL,
    verbose           = FALSE
  )

  # Row-order in the two wide matrices is identical by construction, so the
  # ropls fits should land on the same scores / loadings / VIP vectors.
  expect_equal(nrow(mod_A$scores),   nrow(mod_B$scores))
  expect_equal(mod_A$scores$p1,      mod_B$scores$p1,      tolerance = 1e-8)
  expect_equal(mod_A$scores$p2,      mod_B$scores$p2,      tolerance = 1e-8)
  expect_equal(mod_A$loadings[, 1],  mod_B$loadings[, 1],  tolerance = 1e-8)
  expect_equal(mod_A$vip_scores$VIP, mod_B$vip_scores$VIP, tolerance = 1e-8)
})

test_that("plsda_fit: response_var divergence creates distinct grp_key groups (no silent average)", {
  skip_if_not_installed("ropls")

  # grp_key includes every non-rep_col column, so a replicate whose response
  # label flips relative to its peers becomes its own biological-sample group
  # rather than being silently first-row-averaged into the peer label. Verify
  # this so a future refactor that changes grp_key composition can't
  # reintroduce the silent-drop failure mode.
  set.seed(42)
  n_cyt <- 6L
  d <- data.frame(
    SubjectID = rep(c("S1", "S2", "S3", "S4"), each = 2L),
    group     = c("Control", "Treatment",          # S1 diverges
                  "Control",   "Control",
                  "Treatment", "Treatment",
                  "Treatment", "Treatment"),
    rep_id    = rep(1:2, times = 4L),
    stringsAsFactors = FALSE
  )
  cyto_mat <- matrix(rnorm(nrow(d) * n_cyt, mean = 100, sd = 10),
                     nrow = nrow(d),
                     dimnames = list(NULL, paste0("cyt", seq_len(n_cyt))))
  d <- cbind(d, as.data.frame(cyto_mat))
  lod <- data.frame(cytokine = paste0("cyt", seq_len(n_cyt)), lod = rep(1, n_cyt))
  pre <- plsda_preprocess(
    data = d, cytokine_cols = paste0("cyt", seq_len(n_cyt)),
    metadata_cols = c("SubjectID", "group", "rep_id"),
    lod_lookup = lod, scale_data = FALSE, log_transform = FALSE,
    verbose = FALSE
  )
  mod <- plsda_fit(pre, response_var = "group", permutations = 0,
                   rep_col = "rep_id", verbose = FALSE)
  # 3 subjects collapse to 1 group each (2 reps -> 1); S1 becomes 2 groups
  # because its 2 reps have different 'group' labels.
  expect_equal(nrow(mod$scores), 5L)
  # S1 now appears with BOTH labels.
  s1_rows <- mod$scores[mod$scores$SubjectID == "S1", , drop = FALSE]
  expect_equal(sort(as.character(s1_rows$group)), c("Control", "Treatment"))
})

test_that("plsda_fit updates preprocessing_info to reflect post-aggregation n_samples", {
  skip_if_not_installed("ropls")

  set.seed(1)
  n_subj <- 6L; n_cyt <- 4L; n_rep <- 2L
  subj <- paste0("S", seq_len(n_subj))
  grp  <- rep(c("A", "B"), length.out = n_subj)
  d <- data.frame(
    SubjectID = rep(subj, each = n_rep),
    group     = rep(grp,  each = n_rep),
    rep_id    = rep(seq_len(n_rep), times = n_subj),
    stringsAsFactors = FALSE
  )
  cyto_mat <- matrix(rnorm(nrow(d) * n_cyt, 100, 5),
                     nrow = nrow(d),
                     dimnames = list(NULL, paste0("cyt", seq_len(n_cyt))))
  d <- cbind(d, as.data.frame(cyto_mat))
  lod <- data.frame(cytokine = paste0("cyt", seq_len(n_cyt)), lod = rep(1, n_cyt))
  pre <- plsda_preprocess(
    data = d, cytokine_cols = paste0("cyt", seq_len(n_cyt)),
    metadata_cols = c("SubjectID", "group", "rep_id"),
    lod_lookup = lod, scale_data = FALSE, log_transform = FALSE,
    verbose = FALSE
  )
  expect_equal(pre$preprocessing_info$n_samples, n_subj * n_rep)

  mod <- plsda_fit(pre, response_var = "group", n_components = 2,
                   permutations = 0, rep_col = "rep_id",
                   replicate_agg = "mean", verbose = FALSE)
  expect_equal(mod$preprocessing_info$n_samples, n_subj)
  expect_true(isTRUE(mod$preprocessing_info$replicate_averaged))
  expect_equal(mod$preprocessing_info$replicate_col, "rep_id")
  expect_equal(mod$preprocessing_info$replicate_agg, "mean")
  # n_cytokines unchanged
  expect_equal(mod$preprocessing_info$n_cytokines, n_cyt)
})

test_that("plsda_fit rejects rep_col not in metadata and rep_col == response_var", {
  skip_if_not_installed("ropls")
  test_obj <- create_test_plsda_data()
  pre <- plsda_preprocess(
    data          = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup    = test_obj$lod,
    verbose       = FALSE
  )
  expect_error(
    plsda_fit(pre, response_var = "group", permutations = 0,
              rep_col = "nonexistent", verbose = FALSE),
    "rep_col 'nonexistent' not found"
  )
  expect_error(
    plsda_fit(pre, response_var = "group", permutations = 0,
              rep_col = "group", verbose = FALSE),
    "rep_col and response_var cannot be the same"
  )
})

test_that("plsda_fit replicate_agg = 'median' / 'geomean' dispatch without error", {
  skip_if_not_installed("ropls")
  # Tiny 3-rep dataset; we only need a no-error check, the parity test above
  # exercises the numerics for 'mean'.
  set.seed(1)
  n_subj <- 12L; n_cyt <- 6L; n_rep <- 3L
  subj <- paste0("S", seq_len(n_subj))
  grp  <- rep(c("A", "B"), length.out = n_subj)

  d <- data.frame(
    SubjectID = rep(subj, each = n_rep),
    group     = rep(grp,  each = n_rep),
    rep_id    = rep(seq_len(n_rep), times = n_subj),
    stringsAsFactors = FALSE
  )
  cyto_mat <- matrix(NA_real_, nrow = n_subj * n_rep, ncol = n_cyt,
                     dimnames = list(NULL, paste0("cyt", seq_len(n_cyt))))
  for (i in seq_len(n_subj)) {
    for (r in seq_len(n_rep)) {
      row <- (i - 1L) * n_rep + r
      cyto_mat[row, ] <- rnorm(n_cyt, 100 + 10 * (grp[i] == "B"), 10)
    }
  }
  d <- cbind(d, as.data.frame(cyto_mat))

  lod <- data.frame(cytokine = paste0("cyt", seq_len(n_cyt)), lod = rep(1, n_cyt))
  pre <- plsda_preprocess(
    data = d, cytokine_cols = paste0("cyt", seq_len(n_cyt)),
    metadata_cols = c("SubjectID", "group", "rep_id"),
    lod_lookup = lod, scale_data = FALSE, log_transform = FALSE,
    verbose = FALSE
  )
  # geomean needs strictly positive values; preprocessed expression can be
  # centered or contain negatives depending on transform. With
  # log_transform=FALSE + scale_data=FALSE the matrix is strictly positive.
  expect_no_error(
    plsda_fit(pre, response_var = "group", permutations = 0,
              rep_col = "rep_id", replicate_agg = "median", verbose = FALSE)
  )
  expect_no_error(
    plsda_fit(pre, response_var = "group", permutations = 0,
              rep_col = "rep_id", replicate_agg = "geomean", verbose = FALSE)
  )
})
