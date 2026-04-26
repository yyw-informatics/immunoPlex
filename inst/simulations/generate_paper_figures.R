# =============================================================================
# generate_paper_figures.R — Main-paper simulation figures (Figs 2–4)
# =============================================================================
#
# Produces publication-quality PDFs:
#   - paper_fig2.pdf: Censoring & Preprocessing (4 panels, A–D)
#   - paper_fig3.pdf: Per-Analyte Method Calibration (4 panels, A–D)
#   - paper_fig4.pdf: Panel-Level FDR Control (3 panels, A–C)
#
# Input:   v3 benchmark summary .rds files in cache/
# Output:  cache/figures/paper_fig{2,3,4}.pdf
#
# Usage:
#   Rscript inst/simulations/generate_paper_figures.R [--cache_dir PATH]
#
# =============================================================================

# ---- Dependencies -----------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# Resolve script directory
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
  }
)
source(file.path(SCRIPT_DIR, "helpers_benchmark.R"))


# ---- CLI Arguments ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
ci <- match("--cache_dir", args)
CACHE_DIR <- if (!is.na(ci) && ci < length(args)) args[ci + 1] else
  file.path(SCRIPT_DIR, "cache")

FIG_DIR <- file.path(CACHE_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- Constants ---------------------------------------------------------------

FIG_W  <- 180 / 25.4          # 180 mm -> 7.09 in
FIG_H  <- FIG_W * 3 / 4       # 4:3 aspect -> 5.31 in
BASE   <- 8                   # font size for 180 mm figure

# Publication theme (extends theme_benchmark from helpers)
theme_pub <- function(...) {
  theme_benchmark(base_size = BASE) +
    theme(plot.tag = element_text(face = "bold", size = BASE + 2), ...)
}

# Acceptance band for coverage plots
coverage_band <- function() {
  list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.93, ymax = 0.97,
             fill = "grey90", alpha = 0.5),
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50")
  )
}

# Calibration band for type I error plots
type1_band <- function() {
  list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.025, ymax = 0.075,
             fill = "grey90", alpha = 0.5),
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey50")
  )
}


# ---- Load Data ---------------------------------------------------------------

rds <- function(f) {
  p <- file.path(CACHE_DIR, f)
  if (!file.exists(p)) {
    warning("Missing: ", p, " — panels using this data will be placeholders.")
    return(NULL)
  }
  readRDS(p)
}

message("Loading summary data from: ", CACHE_DIR)

b0         <- rds("preprocessing_v3_summary.rds")
b1         <- rds("censoring_v3_summary.rds")
b1_panel   <- rds("censoring_v3_panel_summary.rds")
b2         <- rds("mcnemar_v2_summary.rds")
b3         <- rds("ancova_v3_summary.rds")
b4         <- rds("plsda_v3_summary.rds")
b5         <- rds("substitution_v3_summary.rds")
b5_panel   <- rds("substitution_v3_panel_summary.rds")
b5_analyte <- rds("substitution_v3_panel_analyte_summary.rds")

# Track figure completeness
fig_status <- list(fig2 = "complete", fig3 = "complete", fig4 = "complete")

message("Data loaded. Missing: ",
        paste(c(if (is.null(b0)) "b0", if (is.null(b1)) "b1",
                if (is.null(b1_panel)) "b1_panel", if (is.null(b2)) "b2",
                if (is.null(b3)) "b3", if (is.null(b4)) "b4",
                if (is.null(b5)) "b5", if (is.null(b5_panel)) "b5_panel",
                if (is.null(b5_analyte)) "b5_analyte"),
              collapse = ", ") %||% "none",
        "\n")


# =============================================================================
# FIGURE 2: Censoring & Preprocessing
# =============================================================================

message("--- Figure 2: Censoring & Preprocessing ---")

## Panels A & B share data: B1, effect = 0.5, backbone n, four key methods
fig2_methods <- c("Oracle", "Tobit", "Gaussian+LOD/2", "Raw")

d2ab <- b1 %>%
  filter(effect_size == 0.5,
         method %in% fig2_methods,
         n_subjects %in% c(40, 100, 200)) %>%
  mutate(cens_pct = censoring_target * 100,
         n_label  = factor(paste0("n = ", n_subjects),
                           levels = paste0("n = ", c(40, 100, 200))))

## A -- Bias vs. censoring
fig2a <- ggplot(d2ab, aes(cens_pct, bias, colour = method, shape = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.2) +
  geom_errorbar(aes(ymin = bias - 1.96 * bias_mcse,
                     ymax = bias + 1.96 * bias_mcse),
                width = 3, linewidth = 0.25) +
  facet_wrap(~n_label, nrow = 1) +
  scale_colour_method() +
  labs(x = "Censoring (%)", y = "Bias", tag = "A") +
  theme_pub(legend.position = "none")

## B -- Coverage vs. censoring
fig2b <- ggplot(d2ab, aes(cens_pct, coverage, colour = method, shape = method)) +
  coverage_band() +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.2) +
  geom_errorbar(aes(ymin = coverage - 1.96 * coverage_mcse,
                     ymax = coverage + 1.96 * coverage_mcse),
                width = 3, linewidth = 0.25) +
  facet_wrap(~n_label, nrow = 1) +
  scale_colour_method() +
  coord_cartesian(ylim = c(NA, 1)) +
  labs(x = "Censoring (%)", y = "95% CI Coverage", tag = "B") +
  theme_pub(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key.size = unit(3, "mm"))

## C -- Method ranking (B5, n = 100, 30% censoring, effect = 0.5, sd = 1.0)
d2c <- b5 %>%
  filter(censoring_target == 0.30,
         n_subjects == 100,
         effect_size == 0.5,
         residual_sd == 1.0) %>%
  arrange(rmse) %>%
  mutate(method = factor(method, levels = method))

fig2c <- ggplot(d2c, aes(rmse, method, colour = method)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = rmse - 1.96 * rmse_mcse,
                      xmax = rmse + 1.96 * rmse_mcse),
                 height = 0.3, linewidth = 0.3) +
  scale_colour_method() +
  labs(x = "RMSE", y = NULL, tag = "C") +
  theme_pub(legend.position = "none")

## D -- Coverage by artifact severity (B0, n = 100, cens = 25%)
d2d <- b0 %>%
  filter(n_subjects == 100,
         censoring_target == 0.25,
         effect_size == 0.5) %>%
  mutate(artifact_severity = factor(artifact_severity,
                                     levels = c("none", "moderate", "heavy")))

fig2d <- ggplot(d2d, aes(artifact_severity, coverage,
                          colour = method, group = method)) +
  coverage_band() +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = coverage - 1.96 * coverage_mcse,
                     ymax = coverage + 1.96 * coverage_mcse),
                width = 0.15, linewidth = 0.25) +
  scale_colour_method() +
  coord_cartesian(ylim = c(NA, 1)) +
  labs(x = "Artifact Severity", y = "95% CI Coverage", tag = "D") +
  theme_pub(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key.size = unit(3, "mm"))

## Assemble Figure 2
fig2 <- (fig2a + fig2b) / (fig2c + fig2d)

ggsave(file.path(FIG_DIR, "paper_fig2.pdf"), fig2,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
message("  Saved: paper_fig2.pdf\n")


# =============================================================================
# FIGURE 3: Per-Analyte Method Calibration
# =============================================================================

message("--- Figure 3: Per-Analyte Method Calibration ---")

## A -- ANCOVA FPR: confounding vs. none (B3, null, no outliers)
d3a <- b3 %>%
  filter(effect_size == 0,
         outlier_rate == 0,
         heteroscedastic == FALSE,
         nonlinear_cov == FALSE,
         outlier_direction == "symmetric",
         confounder_type %in% c("none", "confounding")) %>%
  mutate(confounder_label = ifelse(confounder_type == "none",
                                    "No confounding", "Confounding"))

fig3a <- ggplot(d3a, aes(factor(n_per_group), rejection_rate,
                          colour = method, shape = method)) +
  type1_band() +
  geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = rejection_rate - 1.96 * rejection_mcse,
                     ymax = rejection_rate + 1.96 * rejection_mcse),
                width = 0.3, linewidth = 0.25,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~confounder_label) +
  scale_colour_method() +
  labs(x = "n per group", y = "Type I Error Rate", tag = "A") +
  theme_pub(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key.size = unit(3, "mm")) +
  guides(colour = guide_legend(nrow = 2))

## B -- ANCOVA coverage under outliers (B3, effect = 0.5, confounding, n = 100)
d3b <- b3 %>%
  filter(effect_size == 0.5,
         confounder_type == "confounding",
         heteroscedastic == FALSE,
         nonlinear_cov == FALSE,
         outlier_direction == "symmetric",
         n_per_group == 100)

fig3b <- ggplot(d3b, aes(factor(outlier_rate), coverage,
                          colour = method, shape = method)) +
  coverage_band() +
  geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = coverage - 1.96 * coverage_mcse,
                     ymax = coverage + 1.96 * coverage_mcse),
                width = 0.3, linewidth = 0.25,
                position = position_dodge(width = 0.5)) +
  scale_colour_method() +
  coord_cartesian(ylim = c(NA, 1)) +
  labs(x = "Outlier Rate", y = "95% CI Coverage", tag = "B") +
  theme_pub(legend.position = "none")

## C -- McNemar Type I error calibration (B2, null, averaged over conditions)
if (!is.null(b2)) {
  d3c <- b2 %>%
    filter(delta_rate == 0) %>%
    group_by(method, n_subjects) %>%
    summarise(type1_mean = mean(mean_fdr, na.rm = TRUE),
              type1_se   = sd(mean_fdr, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  fig3c <- ggplot(d3c, aes(n_subjects, type1_mean,
                            colour = method, shape = method)) +
    type1_band() +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = type1_mean - 1.96 * type1_se,
                       ymax = type1_mean + 1.96 * type1_se),
                  width = 5, linewidth = 0.25) +
    scale_colour_method() +
    scale_x_continuous(breaks = c(15, 30, 50, 100, 200)) +
    labs(x = "Sample Size (n)", y = "Type I Error Rate", tag = "C") +
    theme_pub(legend.position = "bottom",
              legend.title = element_blank(),
              legend.key.size = unit(3, "mm"))
} else {
  fig3c <- plot_spacer() +
    plot_annotation(subtitle = "Panel C: McNemar data pending (B2)") &
    theme_pub()
  fig_status$fig3 <- "partial (missing panel C — B2 data pending)"
  message("  [SKIP] Panel C — mcnemar_v2_summary.rds not available")
}

## D -- PLS-DA VIP AUROC (B4, nested CV, representative conditions)
d3d <- b4 %>%
  filter(method == "PLS-DA nested CV",
         n_discriminatory == 5,
         correlation == "none",
         censoring == 0.01,
         imbalance == 1,
         effect_size > 0) %>%
  mutate(effect_label = paste0("d = ", effect_size))

fig3d <- ggplot(d3d, aes(n_per_group, vip_auroc,
                          colour = effect_label, shape = effect_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = vip_auroc - 1.96 * vip_auroc_mcse,
                     ymax = vip_auroc + 1.96 * vip_auroc_mcse),
                width = 3, linewidth = 0.25) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(x = "n per group", y = "VIP AUROC",
       colour = "Effect", shape = "Effect", tag = "D") +
  theme_pub(legend.position = "bottom",
            legend.key.size = unit(3, "mm"))

## Assemble Figure 3
fig3 <- (fig3a + fig3b) / (fig3c + fig3d)

ggsave(file.path(FIG_DIR, "paper_fig3.pdf"), fig3,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
message("  Saved: paper_fig3.pdf\n")


# =============================================================================
# FIGURE 4: Panel-Level FDR Control
# =============================================================================

message("--- Figure 4: Panel-Level FDR Control ---")

## A -- Panel FDR by censoring and method (B1*, effect = 0.5, n = 100)
d4a <- b1_panel %>%
  filter(effect_size == 0.5, n_subjects == 100)

fig4a <- ggplot(d4a, aes(factor(censoring_target), panel_fdr,
                          colour = method, shape = method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = panel_fdr - 1.96 * panel_fdr_mcse,
                     ymax = panel_fdr + 1.96 * panel_fdr_mcse),
                width = 0.25, linewidth = 0.25,
                position = position_dodge(width = 0.4)) +
  scale_colour_method() +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Censoring Rate", y = "Panel FDR", tag = "A") +
  theme_pub(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key.size = unit(3, "mm"))

## B -- Panel sensitivity by method (B1*, effect > 0, n = 100)
## Note: uses b1_panel (always available when Fig 4A builds)
if (!is.null(b1_panel)) {
  d4b <- b1_panel %>%
    filter(effect_size > 0, n_subjects == 100) %>%
    mutate(effect_label = paste0("d = ", effect_size))

  fig4b <- ggplot(d4b, aes(factor(censoring_target), panel_sensitivity,
                            colour = method, shape = method)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = panel_sensitivity - 1.96 * panel_sensitivity_mcse,
                       ymax = panel_sensitivity + 1.96 * panel_sensitivity_mcse),
                  width = 0.25, linewidth = 0.25,
                  position = position_dodge(width = 0.4)) +
    facet_wrap(~effect_label) +
    scale_colour_method() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Censoring Rate", y = "Panel Sensitivity", tag = "B") +
    theme_pub(legend.position = "none")
} else {
  fig4b <- plot_spacer()
  message("  [SKIP] Panel B — b1_panel data not available")
}

## C -- Per-analyte bias vs. realized censoring (B5*, n = 100, effect = 0.5)
if (!is.null(b5_analyte)) {
  fig4c_methods <- c("Tobit", "LOD/2 + Gaussian",
                      "LOD/sqrt2 + Gaussian", "zero + Gaussian")

  d4c <- b5_analyte %>%
    filter(n_subjects == 100,
           effect_size == 0.5,
           method %in% fig4c_methods) %>%
    mutate(analyte_type = ifelse(is_signal, "Signal", "Null"),
           realized_pct = mean_realized_cens * 100)

  fig4c <- ggplot(d4c, aes(realized_pct, mean_bias,
                            colour = analyte_type, shape = analyte_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1.3, alpha = 0.8) +
    facet_wrap(~method, nrow = 1) +
    scale_colour_manual(values = c("Signal" = "#E69F00", "Null" = "#56B4E9")) +
    labs(x = "Realized Censoring (%)", y = "Mean Bias",
         colour = "Analyte", shape = "Analyte", tag = "C") +
    theme_pub(legend.position = "bottom",
              legend.key.size = unit(3, "mm"))
} else {
  fig4c <- plot_spacer()
  fig_status$fig4 <- "partial (missing panel C — B5 panel analyte data pending)"
  message("  [SKIP] Panel C — substitution_v3_panel_analyte_summary.rds not available")
}

## Assemble Figure 4 (2-over-1 layout)
fig4 <- (fig4a | fig4b) / fig4c + plot_layout(heights = c(1, 0.8))

ggsave(file.path(FIG_DIR, "paper_fig4.pdf"), fig4,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
message("  Saved: paper_fig4.pdf\n")


# ---- Summary ----------------------------------------------------------------

message("=== Figure Generation Summary ===")
for (fig in names(fig_status)) {
  message("  ", fig, ": ", fig_status[[fig]])
}
message("\nAll figures saved to: ", FIG_DIR)
