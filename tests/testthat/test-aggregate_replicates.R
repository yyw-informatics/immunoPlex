# Tests for aggregate_replicates()
#
# Covers:
#   - numeric path: mean / median / geomean across single-rep and multi-rep
#   - boolean detection path: majority_vote / any_detected / all_detected
#   - boolean censoring path: semantics match detection path (inverted)
#   - mutually-exclusive target columns validated
#   - NA handling (groups with partial / full missingness)
#   - degenerate single-replicate-per-group case = identity
#   - timepoint_col NULL (cross-sectional) supported
#   - keep_sd emits companion SD column for value_col only
#   - factor levels on group keys are preserved

library(testthat)

# ---- fixture helpers --------------------------------------------------------

make_numeric_long <- function() {
  d <- expand.grid(
    rep_id     = c(1L, 2L),
    timepoint  = c("T1", "T2", "T3"),
    cytokine   = c("A", "B"),
    subject_id = c("S1", "S2"),
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )
  # Pin S1/A/T1 to (10, 20) so mean=15, median=15, geomean=sqrt(200), sd=sqrt(50).
  is_s1_a_t1 <- d$subject_id == "S1" & d$cytokine == "A" & d$timepoint == "T1"
  d$intensity <- ifelse(is_s1_a_t1,
                        ifelse(d$rep_id == 1L, 10, 20),
                        100 + seq_len(nrow(d)))
  d[, c("subject_id", "cytokine", "timepoint", "rep_id", "intensity")]
}

make_bool_long <- function() {
  data.frame(
    subject_id = rep(c("S1", "S2", "S3"), each = 4),
    cytokine   = rep("A", times = 12),
    timepoint  = rep(c("T1", "T1", "T2", "T2"), times = 3),
    rep_id     = rep(1:2, times = 6),
    # S1: T1 = (T,F) disagreement, T2 = (T,T) both detected
    # S2: T1 = (F,F), T2 = (T,F)
    # S3: T1 = (T,T), T2 = (F,F)
    detected   = c(TRUE, FALSE, TRUE, TRUE,
                   FALSE, FALSE, TRUE, FALSE,
                   TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

# ---- numeric path -----------------------------------------------------------

test_that("numeric path: mean/median/geomean each aggregate correctly", {
  d <- make_numeric_long()

  r_mean <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "mean"
  )
  expect_equal(nrow(r_mean), 2 * 2 * 3)  # 2 subjects * 2 cytokines * 3 timepoints
  expect_named(r_mean, c("subject_id", "cytokine", "timepoint", "intensity"))

  s1_a_t1 <- r_mean[r_mean$subject_id == "S1" &
                    r_mean$cytokine   == "A" &
                    r_mean$timepoint  == "T1", ]
  expect_equal(s1_a_t1$intensity, 15)

  r_med <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "median"
  )
  # median of (10, 20) == 15
  expect_equal(r_med$intensity[r_med$subject_id == "S1" &
                               r_med$cytokine   == "A" &
                               r_med$timepoint  == "T1"], 15)

  r_geo <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "geomean"
  )
  expect_equal(r_geo$intensity[r_geo$subject_id == "S1" &
                               r_geo$cytokine   == "A" &
                               r_geo$timepoint  == "T1"],
               sqrt(10 * 20))
})

test_that("numeric path: default rule is 'mean'", {
  d <- make_numeric_long()
  r_default <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity"
  )
  r_mean <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "mean"
  )
  expect_equal(r_default, r_mean)
})

test_that("numeric path: geomean warns and NA-collapses on non-positive input", {
  d <- data.frame(
    subject_id = rep("S1", 4),
    cytokine   = rep(c("A", "B"), each = 2),
    rep_id     = rep(1:2, times = 2),
    intensity  = c(10, 20, -1, 5),
    stringsAsFactors = FALSE
  )
  expect_warning(
    r <- aggregate_replicates(
      d, rep_col = "rep_id",
      subject_col = "subject_id", cytokine_col = "cytokine",
      value_col = "intensity", rule = "geomean"
    ),
    "geomean"
  )
  # B has a non-positive replicate -> NA; A is fine.
  expect_equal(r$intensity[r$cytokine == "A"], sqrt(10 * 20))
  expect_true(is.na(r$intensity[r$cytokine == "B"]))
})

test_that("numeric path: geomean warning is emitted once, not once per offending group", {
  # Pre-fix, each offending group fired its own warning(); a realistic
  # panel of many cytokines / subjects with LOD-adjacent signals could
  # produce hundreds of duplicate warnings. Now the warning is
  # accumulated and emitted once at the end with the offending group
  # count.
  set.seed(1)
  n_grp <- 50L
  d <- data.frame(
    subject_id = paste0("S", seq_len(n_grp)),
    cytokine   = "A",
    rep_id     = 1L,
    intensity  = -seq_len(n_grp),  # every group non-positive
    stringsAsFactors = FALSE
  )
  warns <- testthat::capture_warnings(
    r <- aggregate_replicates(
      d, rep_col = "rep_id",
      subject_col = "subject_id", cytokine_col = "cytokine",
      value_col = "intensity", rule = "geomean"
    )
  )
  expect_length(warns, 1L)
  expect_match(warns, "50 of 50 group")
  expect_true(all(is.na(r$intensity)))
})

test_that("numeric path: keep_sd appends <value_col>_sd column", {
  d <- make_numeric_long()
  r <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "mean", keep_sd = TRUE
  )
  expect_true("intensity_sd" %in% names(r))
  # sd of c(10, 20) == 7.0710678...
  expect_equal(r$intensity_sd[r$subject_id == "S1" &
                              r$cytokine   == "A" &
                              r$timepoint  == "T1"],
               stats::sd(c(10, 20)))
})

# ---- boolean detection path -------------------------------------------------

test_that("detection path: majority_vote / any_detected / all_detected match mcnemar_detection semantics", {
  d <- make_bool_long()

  r_maj <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected",
    rule = "majority_vote"
  )
  r_any <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected",
    rule = "any_detected"
  )
  r_all <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected",
    rule = "all_detected"
  )

  pick <- function(x, sub, tp) {
    x$detected[x$subject_id == sub & x$timepoint == tp]
  }

  # S1,T1 = (T,F). majority_vote: round(0.5)=0 -> FALSE. any=TRUE. all=FALSE.
  expect_false(pick(r_maj, "S1", "T1"))
  expect_true (pick(r_any, "S1", "T1"))
  expect_false(pick(r_all, "S1", "T1"))

  # S1,T2 = (T,T). all rules TRUE.
  expect_true(pick(r_maj, "S1", "T2"))
  expect_true(pick(r_any, "S1", "T2"))
  expect_true(pick(r_all, "S1", "T2"))

  # S2,T1 = (F,F). all rules FALSE.
  expect_false(pick(r_maj, "S2", "T1"))
  expect_false(pick(r_any, "S2", "T1"))
  expect_false(pick(r_all, "S2", "T1"))
})

test_that("detection path: default rule is 'majority_vote'", {
  d <- make_bool_long()
  r_def <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected"
  )
  r_maj <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected",
    rule = "majority_vote"
  )
  expect_equal(r_def, r_maj)
})

test_that("censoring path inverts detection: same group resolves as !detected", {
  d <- make_bool_long()
  d_cens <- d
  d_cens$cens_lod <- !d$detected
  d_cens$detected <- NULL

  r_det_any <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected",
    rule = "any_detected"
  )
  r_cens_any <- aggregate_replicates(
    d_cens, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", censoring_col = "cens_lod",
    rule = "any_detected"
  )

  # Align sort order, then test inversion.
  key <- function(x) paste(x$subject_id, x$timepoint, sep = "|")
  r_det_any  <- r_det_any[order(key(r_det_any)), ]
  r_cens_any <- r_cens_any[order(key(r_cens_any)), ]
  expect_equal(r_cens_any$cens_lod, !r_det_any$detected)
})

# ---- NA handling ------------------------------------------------------------

test_that("fully-missing groups: all three boolean rules return NA", {
  # Prior behavior: majority_vote -> NA, all_detected -> NA, any_detected
  # -> FALSE (via any(..., na.rm=TRUE) on an empty set). The asymmetry
  # made "any_detected" silently code fully-missing subjects as below-LOD
  # rather than dropping them. Now all three align on NA.
  d <- data.frame(
    subject_id = c("S1", "S1"),
    cytokine   = c("A", "A"),
    timepoint  = c("T1", "T1"),
    rep_id     = 1:2,
    detected   = c(NA, NA),
    stringsAsFactors = FALSE
  )
  for (rl in c("majority_vote", "any_detected", "all_detected")) {
    r <- aggregate_replicates(
      d, rep_col = "rep_id",
      subject_col = "subject_id", cytokine_col = "cytokine",
      timepoint_col = "timepoint", detection_col = "detected",
      rule = rl
    )
    expect_true(is.na(r$detected),
                info = paste0("rule = ", rl))
  }
})

test_that("mean ignores NA replicates via na.rm=TRUE", {
  d <- data.frame(
    subject_id = c("S1", "S1", "S1"),
    cytokine   = c("A", "A", "A"),
    timepoint  = c("T1", "T1", "T1"),
    rep_id     = 1:3,
    intensity  = c(10, NA, 20),
    stringsAsFactors = FALSE
  )
  r <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "mean"
  )
  expect_equal(r$intensity, 15)
})

# ---- degenerate / structural cases ------------------------------------------

test_that("single-replicate-per-group input is identity for every rule", {
  d_num <- data.frame(
    subject_id = c("S1", "S1", "S2"),
    cytokine   = c("A", "B", "A"),
    timepoint  = c("T1", "T1", "T1"),
    rep_id     = c(1L, 1L, 1L),
    intensity  = c(10, 20, 30),
    stringsAsFactors = FALSE
  )
  for (rl in c("mean", "median", "geomean")) {
    r <- aggregate_replicates(
      d_num, rep_col = "rep_id",
      subject_col = "subject_id", cytokine_col = "cytokine",
      timepoint_col = "timepoint", value_col = "intensity",
      rule = rl
    )
    expect_equal(nrow(r), 3)
    expect_equal(sort(r$intensity), c(10, 20, 30))
  }

  d_bool <- data.frame(
    subject_id = c("S1", "S1", "S2"),
    cytokine   = c("A", "B", "A"),
    timepoint  = c("T1", "T1", "T1"),
    rep_id     = c(1L, 1L, 1L),
    detected   = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  for (rl in c("majority_vote", "any_detected", "all_detected")) {
    r <- aggregate_replicates(
      d_bool, rep_col = "rep_id",
      subject_col = "subject_id", cytokine_col = "cytokine",
      timepoint_col = "timepoint", detection_col = "detected",
      rule = rl
    )
    expect_equal(nrow(r), 3)
    # Order-agnostic: every row must match its input.
    key_in  <- paste(d_bool$subject_id, d_bool$cytokine, d_bool$timepoint)
    key_out <- paste(r$subject_id,      r$cytokine,      r$timepoint)
    r_ord  <- r[match(key_in, key_out), ]
    expect_equal(r_ord$detected, d_bool$detected)
  }
})

test_that("timepoint_col = NULL groups by (subject, cytokine) only", {
  d <- data.frame(
    subject_id = c("S1", "S1", "S1", "S1"),
    cytokine   = c("A", "A", "B", "B"),
    rep_id     = c(1L, 2L, 1L, 2L),
    intensity  = c(10, 20, 30, 40),
    stringsAsFactors = FALSE
  )
  r <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = NULL, value_col = "intensity", rule = "mean"
  )
  expect_equal(nrow(r), 2)
  expect_false("timepoint" %in% names(r))
  expect_equal(sort(r$intensity), c(15, 35))
})

test_that("mixed-type data frames: non-target columns are dropped from output", {
  d <- data.frame(
    subject_id = c("S1", "S1"),
    cytokine   = c("A", "A"),
    timepoint  = c("T1", "T1"),
    rep_id     = 1:2,
    intensity  = c(10, 20),
    detected   = c(TRUE, FALSE),
    batch      = c("b1", "b1"),
    stringsAsFactors = FALSE
  )
  r <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity", rule = "mean"
  )
  expect_named(r, c("subject_id", "cytokine", "timepoint", "intensity"))
  expect_false("detected" %in% names(r))
  expect_false("batch" %in% names(r))
  expect_false("rep_id" %in% names(r))
})

test_that("factor levels AND values on group keys are preserved", {
  d <- data.frame(
    subject_id = factor(c("S2", "S1", "S2", "S1"), levels = c("S2", "S1")),
    cytokine   = factor(c("A", "A", "A", "A")),
    timepoint  = factor(c("T2", "T2", "T1", "T1"), levels = c("T2", "T1")),
    rep_id     = c(1L, 1L, 1L, 1L),
    intensity  = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  r <- aggregate_replicates(
    d, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity", rule = "mean"
  )
  expect_true(is.factor(r$subject_id))
  expect_equal(levels(r$subject_id), c("S2", "S1"))
  expect_true(is.factor(r$timepoint))
  expect_equal(levels(r$timepoint), c("T2", "T1"))
  # Values (not just levels) survive the round trip. Pre-fix, typeof(factor)
  # returned "integer" and `vector(mode = "integer")` + assignment stripped
  # the labels, so every key came back as NA after re-factoring against
  # character levels.
  expect_equal(sort(as.character(r$subject_id)), c("S1", "S1", "S2", "S2"))
  expect_equal(sort(as.character(r$timepoint)),  c("T1", "T1", "T2", "T2"))
  # And the aggregated values end up in the correct group rows.
  picks <- with(r, paste(subject_id, timepoint))
  expect_equal(r$intensity[match("S1 T2", picks)], 2)
  expect_equal(r$intensity[match("S2 T1", picks)], 3)
})

# ---- input validation -------------------------------------------------------

test_that("exactly one of value_col / detection_col / censoring_col is required", {
  d <- make_numeric_long()
  expect_error(
    aggregate_replicates(d, rep_col = "rep_id",
                         subject_col = "subject_id",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint"),
    "exactly one"
  )
  expect_error(
    aggregate_replicates(d, rep_col = "rep_id",
                         subject_col = "subject_id",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint",
                         value_col = "intensity",
                         detection_col = "intensity"),
    "exactly one"
  )
})

test_that("missing column names are caught up front", {
  d <- make_numeric_long()
  expect_error(
    aggregate_replicates(d, rep_col = "rep_id",
                         subject_col = "nope",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint",
                         value_col = "intensity"),
    "Missing required columns"
  )
})

test_that("value_col must be numeric; detection/censoring must be logical", {
  d <- make_numeric_long()
  d$label <- as.character(d$intensity)
  expect_error(
    aggregate_replicates(d, rep_col = "rep_id",
                         subject_col = "subject_id",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint",
                         value_col = "label"),
    "must be numeric"
  )

  d2 <- make_numeric_long()
  expect_error(
    aggregate_replicates(d2, rep_col = "rep_id",
                         subject_col = "subject_id",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint",
                         detection_col = "intensity"),
    "must be logical"
  )
})

test_that("rule accepts unambiguous partial matches (match.arg semantics)", {
  # API consistency with mcnemar_detection(replicate_agg = ...), which uses
  # match.arg and therefore accepts partials like "maj" -> "majority_vote".
  d_num <- make_numeric_long()
  r_med <- aggregate_replicates(
    d_num, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "med"
  )
  r_full <- aggregate_replicates(
    d_num, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", value_col = "intensity",
    rule = "median"
  )
  expect_equal(r_med, r_full)

  d_b <- make_bool_long()
  r_any <- aggregate_replicates(
    d_b, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected",
    rule = "any"
  )
  r_full_any <- aggregate_replicates(
    d_b, rep_col = "rep_id",
    subject_col = "subject_id", cytokine_col = "cytokine",
    timepoint_col = "timepoint", detection_col = "detected",
    rule = "any_detected"
  )
  expect_equal(r_any, r_full_any)
})

test_that("rule is validated against type-appropriate allow-list", {
  d <- make_numeric_long()
  expect_error(
    aggregate_replicates(d, rep_col = "rep_id",
                         subject_col = "subject_id",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint",
                         value_col = "intensity",
                         rule = "any_detected"),
    "rule.*must be one of"
  )
  d_b <- make_bool_long()
  expect_error(
    aggregate_replicates(d_b, rep_col = "rep_id",
                         subject_col = "subject_id",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint",
                         detection_col = "detected",
                         rule = "geomean"),
    "rule.*must be one of"
  )
})

test_that("keep_sd with boolean path is rejected", {
  d <- make_bool_long()
  expect_error(
    aggregate_replicates(d, rep_col = "rep_id",
                         subject_col = "subject_id",
                         cytokine_col = "cytokine",
                         timepoint_col = "timepoint",
                         detection_col = "detected",
                         keep_sd = TRUE),
    "keep_sd"
  )
})
