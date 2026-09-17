# ============================================================================
# manuscript_figures.R — every manuscript figure, in citation order.
# One function per figure; the plotting helpers live in Plots.R (and
# timeshift_plots.R for Fig 5). Run from the repository root:
#
#   Rscript manuscript_figures.R                 # all 16 figures
#   Rscript manuscript_figures.R fig3 figS8      # selected figures
#
# Running each figure in its own Rscript call avoids an intermittent PDF
# "write failed" when many large figures share one R session:
#   for f in fig1 fig2 fig3; do Rscript manuscript_figures.R $f; done
#
# Data (see README for the scripts that produce each tree):
#   results/minimal, results/acute   main runs, 4 reps       scripts/run_main.sh
#   results/tracking                 tracking sweep, 8 reps  scripts/run_tracking_sweep.sh
#   results/timeshift                time-shift summaries    scripts/run_timeshift.sh
#   results_gamma/                   tempo sweep, N_H = N_P  scripts/run_tempo_sweep.sh
#   results_zoom/                    write_every = 1 runs    scripts/run_zoom.sh
# Output: figures/<name>.pdf and .png
# ============================================================================

source("Plots.R")

options(ggt.results_root = "results")
OUT_DIR <- "figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SIGMA <- 0.01

ZOOM_ROOT     <- "results_zoom"      # write_every = 1, tagged "zoom"
TEMPO_ROOT    <- "results_gamma"
TIMESHIFT_DIR <- "results/timeshift"

# Stop early, naming the script to run, when a figure's input is missing
require_data <- function(path, script) {
  if (!file.exists(path)) stop(path, " not found; run ", script, " first", call. = FALSE)
}

save_fig <- function(p, name, width, height) {
  safe_ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
  safe_ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  cat("Saved:", name, "\n")
  invisible(p)
}

# Bold A, B, C ... tags in reading order
tagged <- function(p) {
  p + plot_annotation(tag_levels = "A") &
    .tag_theme &
    theme(plot.margin = margin(14, 14, 6, 18))
}

inset_legend <- theme(
  legend.position        = "inside",
  legend.position.inside = c(0.98, 0.98),
  legend.justification   = c(1, 1),
  legend.background      = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
  legend.key             = element_blank(),
  legend.margin          = margin(2, 4, 2, 4)
)


# ---------------------------------------------------------------------------
# Main text
# ---------------------------------------------------------------------------

# Fig 1: model and geometry (analytical, no simulation data)
fig1 <- function() {
  p <- fig_grouped_1_setup("minimal", filename = NULL)
  save_fig(p, "Fig1_model", 13, 9)
}

# Fig 2: minimal-model time series, v, c, W_P, W_H x four scenarios.
# The x axis counts substitutions (gen), not evolutionary time: over the same
# 100K substitutions ET/ET covers ~1e4 times more evolutionary time than ER/ER.
fig2 <- function() {
  require_data("results/minimal", "scripts/run_main.sh")
  p <- fig_timeseries("minimal", sigma = SIGMA, diploid = TRUE,
                      replicates = TRUE, max_pts = 400, highlight_rep = 1,
                      show_omega = FALSE, show_title = FALSE)
  save_fig(p, "Fig2_timeseries", 11, 2.2 * 4)
}

# Fig 3: mechanism of destabilisation, drift -> amplification -> threshold -> boundary
fig3 <- function(model = "minimal") {
  require_data(file.path("results", model), "scripts/run_main.sh")
  ev <- load_neutral_events(model, SIGMA)
  pA <- fig_neutral_fraction(ev)
  pB <- fig_decoupling(ev)
  pC <- fig_step_sizes(model, sigma = SIGMA, diploid = TRUE, filename = NULL) +
    labs(title = NULL)
  pD <- fig_slope_distribution(model_name = model, sigma = SIGMA, diploid = TRUE,
                               title_label = NULL, pct_inline = TRUE, lim_q = 0.05) +
    inset_legend + theme(aspect.ratio = 1)
  pE <- fig_nash_violation_map(model_name = model, sigma = SIGMA, diploid = TRUE,
                               resolution = 100) +
    labs(title = NULL)
  pF <- fig_trait_occupancy(model_name = model, sigma = SIGMA, diploid = TRUE,
                            weight = "time")
  p <- tagged(((pA | pB) / pC / (pD | pE) / pF) +
                plot_layout(heights = c(1, 1, 1.1, 0.9)))
  save_fig(p, "Fig3_mechanism", 13, 21)
}

# Fig 4: selection and tempo.
#   A omega/2N distributions, B per-replicate medians, C aligned zoom,
#   D-E tempo sweep
fig4 <- function(model = "minimal") {
  require_data(file.path("results", model), "scripts/run_main.sh")
  require_data(ZOOM_ROOT, "scripts/run_zoom.sh")
  require_data(TEMPO_ROOT, "scripts/run_tempo_sweep.sh")
  pA <- fig_omega(model_name = model, sigma = SIGMA, diploid = TRUE,
                  filename = NULL) +
    labs(title = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))
  pB <- fig_discriminator(model_name = model, sigma = SIGMA, diploid = TRUE,
                          tag_prefix = NULL, filename = NULL)
  pC <- wrap_elements(full = fig_aligned_zoom(model, SIGMA, results_root = ZOOM_ROOT,
                                              tag_prefix = "zoom"))
  tp <- tempo_panels(TEMPO_ROOT)
  # One legend for both tempo panels (same shapes and linetypes)
  p <- tagged(((pA | pB) / pC /
                 (tp$fitness | (tp$boundary + theme(legend.position = "none")))) +
                plot_layout(heights = c(1, 1.3, 1)))
  save_fig(p, "Fig4_selection_tempo", 13, 22)
}

# Fig 5: time-shift assays (zoom summaries: panels A-C; main: panels D-E)
fig5 <- function(model = "minimal") {
  for (f in c("timeshift_zoom.csv", "timeshift_zoom_pairs.csv",
              "timeshift_main.csv", "timeshift_main_pairs.csv"))
    require_data(file.path(TIMESHIFT_DIR, f), "scripts/run_timeshift.sh main|zoom")
  source("timeshift_plots.R")
  TS_DIR <<- TIMESHIFT_DIR

  # This layout was designed around timeshift_plots.R's 12-pt theme; the larger
  # manuscript fonts make its row titles and +-200 tick labels collide
  p <- fig_timeshift(model = model, filename = NULL) &
    theme(axis.text   = element_text(size = 13, colour = "black"),
          axis.title  = element_text(size = 14),
          strip.text  = element_text(size = 14),
          plot.title  = element_text(size = 15),
          legend.text = element_text(size = 13))
  save_fig(p, "Fig5_timeshift", 14, 15)
}


# ---------------------------------------------------------------------------
# Supplementary
# ---------------------------------------------------------------------------

# S1: realised-trait (top) and fitness (bottom) density hexbins, minimal model
figS1 <- function() {
  # Count legend hidden: its labels are wrong on the log fill scale; the caption
  # states that colour is log-scaled density
  p <- fig_grouped_2_dynamics("minimal", sigma = SIGMA, diploid = TRUE,
                              tag_prefix = NULL, filename = NULL) &
    theme(legend.position = "none")
  save_fig(p, "FigS1_density", 13, 8)
}

# S2: zoomed time series, first 1,000 substitutions, from the zoom runs
figS2 <- function() {
  require_data(ZOOM_ROOT, "scripts/run_zoom.sh")
  p <- fig_timeseries("minimal", sigma = SIGMA, diploid = TRUE, tag_prefix = "zoom",
                      results_root = ZOOM_ROOT, window = c(0, 1000),
                      highlight_rep = 1, show_omega = FALSE, show_title = FALSE)
  save_fig(p, "FigS2_timeseries_zoom", 11, 2.2 * 4)
}

# S3: time-series statistics, trimmed (mean, SD, spectral slope; corr. length for c)
figS3 <- function() {
  save_fig(fig_ts_stats_si("minimal", SIGMA), "FigS3_ts_stats", 11, 13)
}

# S4: response-rule snapshots, ER/ER
figS4 <- function() {
  p <- fig_snapshots(model_name = "minimal", sigma = SIGMA, diploid = TRUE,
                     filename = NULL, show_title = FALSE)
  save_fig(p, "FigS4_rule_snapshots", 13, 9)
}

# S5: response-rule parameter time series, ER/ER
figS5 <- function() {
  # Pseudo-log y: rare slope and intercept spikes of several hundred would
  # otherwise flatten the everyday drift near 0 (replaces the linear limits)
  p <- fig_strategy_evolution(model_name = "minimal", sigma = SIGMA, diploid = TRUE,
                              filename = NULL) &
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(-1000, -100, -10, 0, 10, 100, 1000))
  save_fig(p, "FigS5_rule_parameters", 12, 8)
}

# S6: full 100K omega / 2N time series, four scenarios
figS6 <- function(model = "minimal") {
  p <- tagged(fig_omega_timeseries(model, SIGMA, normalize = TRUE))
  save_fig(p, "FigS6_omega_timeseries", 13, 7)
}

# S7: tempo symmetry check at R = 1 (doubles as a code check)
figS7 <- function() {
  require_data(TEMPO_ROOT, "scripts/run_tempo_sweep.sh")
  save_fig(fig_tempo_symmetry(TEMPO_ROOT, SIGMA), "FigS7_tempo_symmetry", 11, 9)
}

# S8: tracking sweep, A virulence SD, B corr(W_H, W_P), C ER-host fitness vs ET/ET
figS8 <- function() {
  require_data("results/tracking", "scripts/run_tracking_sweep.sh")
  save_fig(tagged(fig_tracking_sweep_si(TRACKING_KS, SIGMA)), "FigS8_tracking_sweep", 16, 5.5)
}

# S9: tracking-model ER/ER time series, one column per k
figS9 <- function() {
  require_data("results/tracking", "scripts/run_tracking_sweep.sh")
  p <- fig_tracking_timeseries_ksweep(condition = "ERhost_ERpath", ks = TRACKING_KS,
                                      sigma = SIGMA, replicates = TRUE, max_pts = 400,
                                      highlight_rep = 1, show_omega = FALSE,
                                      show_title = FALSE, filename = NULL) &
    theme(plot.tag = element_blank())   # 28 panels outrun A-Z; columns carry k
  save_fig(p, "FigS9_tracking_timeseries", 2.5 * length(TRACKING_KS) + 1, 2.2 * 4)
}

# S10: acute-model time series, same layout as Fig 2
figS10 <- function() {
  require_data("results/acute", "scripts/run_main.sh")
  p <- fig_timeseries("acute", sigma = SIGMA, diploid = TRUE,
                      replicates = TRUE, max_pts = 400, highlight_rep = 1,
                      show_omega = FALSE, show_title = FALSE)
  save_fig(p, "FigS10_acute_timeseries", 11, 2.2 * 4)
}

# S11: acute-model diagnostics, analogues of Fig 3C, Fig 3F and Fig 4B
figS11 <- function(model = "acute") {
  require_data(file.path("results", model), "scripts/run_main.sh")
  pA <- fig_step_sizes(model, sigma = SIGMA, diploid = TRUE, filename = NULL) +
    labs(title = NULL)
  pB <- fig_trait_occupancy(model_name = model, sigma = SIGMA, diploid = TRUE,
                            weight = "time")
  pC <- fig_discriminator(model_name = model, sigma = SIGMA, diploid = TRUE,
                          tag_prefix = NULL, filename = NULL)
  p <- tagged((pA / (pB | pC)) + plot_layout(heights = c(1, 1)))
  save_fig(p, "FigS11_acute_diagnostics", 14, 13)
}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
FIGURES <- list(fig1 = fig1, fig2 = fig2, fig3 = fig3, fig4 = fig4, fig5 = fig5,
                figS1 = figS1, figS2 = figS2, figS3 = figS3, figS4 = figS4,
                figS5 = figS5, figS6 = figS6, figS7 = figS7, figS8 = figS8,
                figS9 = figS9, figS10 = figS10, figS11 = figS11)

if (sys.nframe() == 0L) {           # run as a script, not when source()d
  todo <- commandArgs(trailingOnly = TRUE)
  if (length(todo) == 0) todo <- names(FIGURES)
  unknown <- setdiff(todo, names(FIGURES))
  if (length(unknown)) stop("Unknown figure(s): ", paste(unknown, collapse = ", "))
  for (f in todo) FIGURES[[f]]()
}
