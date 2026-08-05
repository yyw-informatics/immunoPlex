# Data documentation for immunoPlex package

#' Example immunoPlex cytokine dataset
#'
#' A synthetic multiplex immunoassay dataset built with
#' \code{\link{simulate_immunoassay}} to demonstrate and test the
#' immunoPlex workflow. Parameters were chosen against the v4 benchmark
#' suite (B0, B5) so the dataset highlights where censoring-aware models
#' (Tobit, AFT) outperform LOD-substitution alternatives. The simulation
#' targets ~25% per-analyte left censoring -- the regime where
#' substitution methods incur visible bias and coverage loss while Tobit
#' tracks the truth.
#'
#' @format A list with three components:
#' \describe{
#'   \item{expression}{Data frame, 480 samples x 20 cytokines, raw
#'     concentrations on the linear scale. Sample identifiers are
#'     stored in \code{rownames}. Cytokines: IL-1b, IL-2, IL-4, IL-6,
#'     IL-8, IL-10, IL-12p70, IL-13, IL-17A, TNF-a, IFN-g, MCP-1,
#'     MIP-1a, MIP-1b, RANTES, IP-10, Eotaxin, G-CSF, VEGF, EGF.}
#'   \item{metadata}{Data frame with 480 rows and 7 columns:
#'     \itemize{
#'       \item \code{sample_id}: unique sample identifier (rowname of
#'         \code{expression})
#'       \item \code{subject_id}: subject identifier shared across
#'         timepoints and replicates
#'       \item \code{timepoint}: \code{"T1"} or \code{"T2"}
#'       \item \code{disease}: \code{"Yes"} or \code{"No"} (60 subjects each)
#'       \item \code{age}: integer age in years (covariate)
#'       \item \code{preexist}: \code{"Yes"} or \code{"No"} pre-existing
#'         condition flag (binary covariate)
#'       \item \code{replicate}: technical replicate index (1 or 2)
#'     }}
#'   \item{lod_lookup}{Data frame with 20 rows (one per cytokine) and
#'     columns \code{cytokine}, \code{lod}, \code{ulod}. Upper LOD is
#'     set for high-abundance saturating analytes (IL-8, MIP-1b,
#'     RANTES, IP-10) and \code{NA} otherwise.}
#' }
#'
#' @details
#' The simulation design (see \code{data-raw/immunoplex_example.R}):
#' \itemize{
#'   \item 120 subjects (60 per disease group), 2 timepoints, 2
#'     technical replicates -- 480 samples total.
#'   \item Block-correlated analytes following the immunoassay
#'     categories (Th1 / Th2 / chemokine / pro-inflammatory /
#'     anti-inflammatory / growth factor) at \code{block_rho = 0.5}.
#'   \item Six signal analytes (IL-6, IL-8, TNF-a, IFN-g, MCP-1, IP-10)
#'     with disease-up effect (\code{group_effects = 1.2} on log scale),
#'     a mild T1 -> T2 trend, and a positive disease x timepoint
#'     interaction. The remaining 14 analytes are null.
#'   \item Per-analyte LODs set at the 25% quantile of the marginal
#'     distribution, yielding ~25% left censoring per cytokine.
#'   \item Censored expression values fall just below the reported LOD
#'     so the strict left-censoring flag in
#'     \code{\link{prepare_cytokine_data}} fires correctly.
#' }
#'
#' Use this dataset to exercise the package end-to-end:
#' \itemize{
#'   \item Preprocessing with \code{\link{immuno_preprocess}}
#'   \item Censoring-aware fits with \code{\link{fit_one}} /
#'     \code{\link{fit_models}}
#'   \item LOD-handling comparisons with
#'     \code{\link{compare_lod_models}}
#'   \item Detection-frequency tests with
#'     \code{\link{mcnemar_detection}}
#'   \item Rank-based ANCOVA with \code{\link{ancova_one}} /
#'     \code{\link{ancova_fit}}
#'   \item Multivariate discrimination with \code{\link{plsda_fit}}
#' }
#'
#' @seealso \code{\link{simulate_immunoassay}},
#'   \code{\link{prepare_cytokine_data}},
#'   \code{\link{list_cytokines}}, \code{\link{immuno_preprocess}}
#'
#' @examples
#' data("immunoplex_example", package = "immunoPlex")
#'
#' # Available cytokines and their realised censoring
#' list_cytokines()
#'
#' # Single-cytokine modeling frame
#' dat <- prepare_cytokine_data("IL-6")
#' attr(dat, "pct_censored")
#'
#' \dontrun{
#' processed <- immuno_preprocess(
#'   expr        = immunoplex_example$expression,
#'   meta        = immunoplex_example$metadata,
#'   lod_lookup  = immunoplex_example$lod_lookup,
#'   lod_methods = "half"
#' )
#' }
"immunoplex_example"
