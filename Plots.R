# ============================================================================
# Goldstein et al. Game Theory — Plotting & Analysis
# Canan Karakoc
# Refactored: February 2026
# ============================================================================
#
# STRUCTURE:
#   §0  Setup (libraries, paths, theme)
#   §1  Fitness functions (all 4 models, defined ONCE)
#   §2  Utility functions (helpers used across figures)
#   §3  Data loading (unified loader for all experiments)
#   §4  FIGURE 1  — Fitness landscapes (any model)
#   §5  FIGURE 2  — Strategy lines & Nash equilibria
#   §6  FIGURE 3  — Time series (v, c, W_H, W_P, omega)
#   §7  FIGURE 4  — Phase-space density (hex plots on landscapes)
#   §8  FIGURE 5  — Snapshots of strategy lines (fig_snapshots)
#   §9  FIGURE 6  — Nash region & slope distribution (fig_nash_combined)
#   §10 FIGURE 7  — Strategy parameter evolution & stability
#   §11 Time series statistics (CV, spectral slope, ACF, memory)
#   §12 Step size analysis
#   §13 Neutral drift analysis
#   §14 Boundary analysis
#   §15 Cross-condition comparison figures (ts stats, steps, drift, dwell)
#
# ADDING A NEW FITNESS MODEL:
#   1. Add functions in §1 (mortality, host_fitness, path_fitness)
#   2. Register in the FITNESS_MODELS list at end of §1
#   3. Data loading in §3 picks it up automatically if directory matches
#
# ============================================================================


# ============================================================================
# §0  SETUP
# ============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(pracma)
library(scales)
library(zoo)

# --- Paths ---
# Change this one line to point at your repo root:
REPO_ROOT <- "~/Documents/GitHub/GoldsteinGameTheory"
setwd(REPO_ROOT)
dir.create("figures", showWarnings = FALSE)

# --- Global theme ---
mytheme <- theme_bw() +
  theme(
    axis.ticks.length   = unit(0.2, "cm"),
    legend.text         = element_text(size = 12),
    axis.text           = element_text(size = 14, color = "black"),
    axis.title          = element_text(size = 14),
    plot.title          = element_text(size = 14),
    panel.border        = element_rect(fill = NA, colour = "black", linewidth = 1),
    strip.text.x        = element_text(size = 14),
    strip.background    = element_blank(),
    legend.title        = element_blank(),
    axis.text.x.top     = element_blank(),
    axis.title.x.top    = element_blank(),
    axis.text.y.right   = element_blank(),
    axis.title.y.right  = element_blank(),
    axis.title.x        = element_text(margin = margin(16, 0, 0)),
    axis.title.y        = element_text(margin = margin(0, 16, 0, 0)),
    axis.text.x         = element_text(margin = margin(16, 0, 0, 0)),
    axis.text.y         = element_text(margin = margin(0, 16, 0, 0))
  )

set_theme(mytheme)


# ============================================================================
# §1  FITNESS FUNCTIONS — all models defined ONCE
# ============================================================================
#
# Notation (paper -> code):
#   c (clearance)  -> s      v (virulence) -> v
#   nS, nV         -> cost parameters
#   d0             -> baseline mortality
#   eps            -> boundary softness
#   beta           -> transmission exponent
#
# Convention: all functions take (v, s) with v first, s second.
# This matches the simulation output column order.

# ----- Shared parameter sets -----

# Acute model parameters (from paper)
PAR_ACUTE <- list(
  d0      = 0.1,
  nS      = 0.1,      # clearance cost (immunopathology)
  nV      = 1.0,      # virulence cost (10x clearance -- biologically motivated)
  eps     = 1e-3,     # boundary softness
  one_eps = 1 + 1e-3, # 1 + eps precomputed
  beta    = 1.0       # transmission exponent (linear)
)

# Chronic model parameters (from paper)
PAR_CHRONIC <- list(
  d0      = 0.1,
  nS      = 0.1,
  nV      = 1.0,
  eps     = 1e-3,
  one_eps = 1 + 1e-3,
  beta    = 1.0
)

# Minimal model (no extra parameters)
PAR_MINIMAL <- list()

# Taylor model
PAR_TAYLOR <- list(
  b  = 1.0,
  m0 = 1.0,
  n  = 0.75
)


# ----- Mortality functions -----

mortality_acute <- function(v, s, p = PAR_ACUTE) {
  p$d0 +
    p$nS * p$one_eps * s / (p$one_eps - s) +
    p$nV * p$one_eps * v / (p$one_eps - v)
}

mortality_chronic <- function(v, s, p = PAR_CHRONIC) {
  # Immunity modulates virulence damage via (1-s)
  p$d0 +
    p$nS * p$one_eps * s / (p$one_eps - s) +
    (1 - s) * p$nV * p$one_eps * v / (p$one_eps - v)
}


# ----- Host fitness -----

fH_acute <- function(v, s, p = PAR_ACUTE) {
  m <- mortality_acute(v, s, p)
  s / (s + m)
}

fH_chronic <- function(v, s, p = PAR_CHRONIC) {
  m <- mortality_chronic(v, s, p)
  1.0 / m
}

fH_minimal <- function(v, s, ...) {
  s * (1 - s) * (1 - v)
}

fH_taylor <- function(v, s, p = PAR_TAYLOR) {
  (s / (v + s)) * (p$b / (p$m0 + s))
}


# ----- Pathogen fitness -----

fP_acute <- function(v, s, p = PAR_ACUTE) {
  m <- mortality_acute(v, s, p)
  v^p$beta / (s + m)
}

fP_chronic <- function(v, s, p = PAR_CHRONIC) {
  m <- mortality_chronic(v, s, p)
  (1 - s) * v^p$beta / m
}

fP_minimal <- function(v, s, ...) {
  v * (1 - v) * (1 - s)
}

fP_taylor <- function(v, s, p = PAR_TAYLOR) {
  v^p$n / (v + s)
}


# ----- Registry: look up functions by model name -----
# Each entry: list(fH, fP, params, label)

FITNESS_MODELS <- list(
  acute = list(
    fH     = fH_acute,
    fP     = fP_acute,
    params = PAR_ACUTE,
    label  = "Acute"
  ),
  chronic = list(
    fH     = fH_chronic,
    fP     = fP_chronic,
    params = PAR_CHRONIC,
    label  = "Chronic"
  ),
  minimal = list(
    fH     = fH_minimal,
    fP     = fP_minimal,
    params = PAR_MINIMAL,
    label  = "Minimal"
  ),
  taylor = list(
    fH     = fH_taylor,
    fP     = fP_taylor,
    params = PAR_TAYLOR,
    label  = "Taylor"
  )
)

# Trait domain per model.  Taylor uses rates (unbounded); others use [0,1].
TRAIT_DOMAIN <- list(
  acute   = c(0.001, 0.999),
  chronic = c(0.001, 0.999),
  minimal = c(0.001, 0.999),
  taylor  = c(0.01,  30.0)     # Nash ≈ (v*=9, c*=3)
)

# Clean axis limits for plotting (not the simulation clamp bounds)
TRAIT_DISPLAY <- list(
  acute   = c(0, 1),
  chronic = c(0, 1),
  minimal = c(0, 1),
  taylor  = c(0, 30)
)

# ============================================================================
# §2  UTILITY FUNCTIONS
# ============================================================================

# Model-aware clamping (defaults to [0,1])
clamp_trait <- function(x, model_name = NULL) {
  if (is.null(model_name)) {
    return(pmin(1, pmax(0, x)))
  }
  dom <- TRAIT_DOMAIN[[model_name]]
  pmin(dom[2], pmax(dom[1], x))
}

# Keep old name as alias for backward compat
clamp01 <- function(x) pmin(1, pmax(0, x))

# --- Axis stripping helpers (for multi-panel layouts) ---
strip_y <- function(p) {
  p + labs(y = NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

strip_x <- function(p) {
  p + labs(x = NULL) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
}

# --- Fitness landscape grid ---
make_fitness_grid <- function(model_name, resolution = 300) {
  mod <- FITNESS_MODELS[[model_name]]
  dom <- TRAIT_DOMAIN[[model_name]]
  vs <- seq(dom[1], dom[2], length.out = resolution)
  ss <- seq(dom[1], dom[2], length.out = resolution)
  grid <- expand.grid(v = vs, s = ss)
  grid$fH <- mapply(mod$fH, grid$v, grid$s)
  grid$fP <- mapply(mod$fP, grid$v, grid$s)
  grid$joint <- grid$fH * grid$fP
  grid
}

# --- Best response curves ---
calc_best_responses <- function(model_name, n = 300) {
  mod <- FITNESS_MODELS[[model_name]]
  dom <- TRAIT_DOMAIN[[model_name]]
  
  v_vals <- seq(dom[1], dom[2], length.out = n)
  s_vals <- seq(dom[1], dom[2], length.out = n)
  
  host_br <- data.frame(
    v = v_vals,
    s = sapply(v_vals, function(vv) {
      optimize(function(ss) -mod$fH(vv, ss), dom)$minimum
    })
  )
  
  path_br <- data.frame(
    s = s_vals,
    v = sapply(s_vals, function(ss) {
      optimize(function(vv) -mod$fP(vv, ss), dom)$minimum
    })
  )
  
  list(host = host_br, path = path_br)
}


# --- Nash equilibrium (brute-force intersection of best responses) ---
find_nash <- function(model_name, n = 300) {
  br <- calc_best_responses(model_name, n)
  host_df <- br$host
  path_df <- br$path
  
  dom <- TRAIT_DOMAIN[[model_name]]
  best_dist <- Inf
  nash <- c(v = mean(dom), s = mean(dom))
  
  for (i in seq_len(nrow(host_df))) {
    vi <- host_df$v[i]
    si <- host_df$s[i]
    dists <- (path_df$v - vi)^2 + (path_df$s - si)^2
    j <- which.min(dists)
    d <- sqrt(dists[j])
    if (d < best_dist) {
      best_dist <- d
      nash <- c(v = vi, s = si)
    }
  }
  
  data.frame(v = nash["v"], s = nash["s"], row.names = NULL)
}

# --- Thin simulation data for faster plotting ---
# Adaptive: targets ~max_pts points. Short runs keep all data.
thin_for_plot <- function(df, every = NULL, max_pts = 2000) {
  out <- df %>% filter(event == "post")
  n <- nrow(out)
  if (is.null(every)) {
    every <- max(1, floor(n / max_pts))
  }
  if (every <= 1) return(out)
  out %>%
    slice(seq(every, n(), by = every))
}

# --- Numerical gradient ---
num_grad <- function(f, v, s, h = 1e-5) {
  v1 <- clamp01(v - h); v2 <- clamp01(v + h)
  s1 <- clamp01(s - h); s2 <- clamp01(s + h)
  dv <- (f(v2, s) - f(v1, s)) / max(v2 - v1, 1e-12)
  ds <- (f(v, s2) - f(v, s1)) / max(s2 - s1, 1e-12)
  list(dv = dv, ds = ds)
}


# ============================================================================
# §3  DATA LOADING (auto-discovery from config.json)
# ============================================================================
#
# Scans results/ for any directory containing config.json + simulation.csv.
# Extracts condition, fitness model, step size, diploid flag, pinned traits
# directly from the JSON — no hardcoded paths to maintain.
#
# Usage:
#   catalog <- discover_experiments()          # scan everything
#   catalog <- discover_experiments("results") # explicit root
#   View(catalog)                              # see what's available
#
#   # Load one experiment:
#   df <- load_sim_by_row(catalog, 1)
#
#   # Load a filtered set:
#   acute_et <- catalog %>% filter(fitness == "acute", condition == "EThost_ETpath")
#   dfs <- load_sim_set(acute_et)
#
#   # Load everything (careful with memory):
#   all_data <- load_all_discovered()

library(jsonlite)

# Global: minimum generation to include when loading data.
# Quick runs (10K gens) need 0; full runs (1M gens) use 10000.
# run_all_figures() auto-sets this; or set manually before loading.
MIN_GEN_CUTOFF <- 0

# Cached catalog — avoid re-scanning filesystem on every load
.catalog_cache <- new.env(parent = emptyenv())
.catalog_cache$data <- NULL
.catalog_cache$root <- NULL

#' Scan results/ tree and build a catalog of all experiments.
#' Results are cached; call discover_experiments(refresh = TRUE) to re-scan.
discover_experiments <- function(results_root = "results", refresh = FALSE) {
  
  # Return cache if valid
  if (!refresh && !is.null(.catalog_cache$data) &&
      identical(.catalog_cache$root, results_root)) {
    return(.catalog_cache$data)
  }
  
  # Find all config.json files
  configs <- list.files(results_root, pattern = "config\\.json$",
                        recursive = TRUE, full.names = TRUE)
  
  if (length(configs) == 0) {
    warning("No config.json files found in ", results_root,
            "\n  Working directory: ", getwd())
    return(tibble())
  }
  
  rows <- lapply(configs, function(cf) {
    tryCatch({
      cfg <- fromJSON(cf)
      dir_path <- dirname(cf)
      csv_path <- file.path(dir_path, "simulation.csv")
      
      tibble(
        dir           = dir_path,
        csv           = csv_path,
        csv_exists    = file.exists(csv_path),
        fitness       = cfg$fitness_model %||% NA_character_,
        condition     = cfg$condition %||% NA_character_,
        host_reactive = cfg$host_reactive %||% NA,
        path_reactive = cfg$path_reactive %||% NA,
        std_dev_move  = cfg$std_dev_move %||% NA_real_,
        gamma         = cfg$prob_host_mutate %||% NA_real_,
        diploid       = isTRUE(cfg$DIPLOID_KIMURA),
        rep           = if (is.null(cfg$rep)) NA_integer_ else as.integer(cfg$rep),
        tag           = cfg$tag %||% NA_character_,
        effective_seed = if (is.null(cfg$effective_seed)) NA_integer_ else as.integer(cfg$effective_seed),
        fix_host      = if (is.null(cfg$FIX_HOST_TRAIT) || 
                            isFALSE(cfg$FIX_HOST_TRAIT)) NA_real_ 
        else as.numeric(cfg$FIX_HOST_TRAIT),
        fix_path      = if (is.null(cfg$FIX_PATH_TRAIT) || 
                            isFALSE(cfg$FIX_PATH_TRAIT)) NA_real_ 
        else as.numeric(cfg$FIX_PATH_TRAIT),
        host_pop      = cfg$HOST_POP_N %||% NA_integer_,
        path_pop      = cfg$PATH_POP_N %||% NA_integer_,
        max_gens      = cfg$parameters$max_gens %||% NA_integer_,
        timestamp     = cfg$timestamp %||% NA_character_
      )
    }, error = function(e) {
      warning("Failed to parse: ", cf, " — ", e$message)
      NULL
    })
  })
  
  result <- bind_rows(rows) %>%
    filter(csv_exists) %>%
    arrange(fitness, condition, std_dev_move)
  
  # Cache
  .catalog_cache$data <- result
  .catalog_cache$root <- results_root
  
  cat("Discovered", nrow(result), "experiments in", results_root, "\n")
  result
}


#' Force re-scan of experiments (call after running new simulations)
refresh_catalog <- function(results_root = "results") {
  discover_experiments(results_root, refresh = TRUE)
}


#' Print a readable summary of all discovered experiments
list_experiments <- function(model = NULL, results_root = "results") {
  cat <- discover_experiments(results_root)
  if (!is.null(model)) cat <- cat %>% filter(fitness == model)
  
  if (nrow(cat) == 0) {
    message("No experiments found. Check working directory: ", getwd())
    return(invisible(cat))
  }
  
  cat("\n", strrep("─", 70), "\n")
  for (mod in unique(cat$fitness)) {
    sub <- cat %>% filter(fitness == mod)
    cat(sprintf("\n  %s  (%d experiments)\n", toupper(mod), nrow(sub)))
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      flags <- c()
      if (isTRUE(row$diploid)) flags <- c(flags, "diploid")
      if (!is.na(row$std_dev_move) && abs(row$std_dev_move - 0.1) > 1e-6)
        flags <- c(flags, sprintf("σ=%g", row$std_dev_move))
      if (!is.na(row$fix_host)) flags <- c(flags, sprintf("fixH=%.2f", row$fix_host))
      if (!is.na(row$fix_path)) flags <- c(flags, sprintf("fixP=%.2f", row$fix_path))
      if (!is.na(row$rep)) flags <- c(flags, sprintf("rep%d", row$rep))
      if (!is.na(row$tag)) flags <- c(flags, sprintf("tag=%s", row$tag))
      if (!is.na(row$max_gens) && row$max_gens != 1000000)
        flags <- c(flags, sprintf("%dK gens", row$max_gens / 1000))
      ftag <- if (length(flags) > 0) paste0("  [", paste(flags, collapse=", "), "]") else ""
      cat(sprintf("    %-25s %s\n", row$condition, ftag))
    }
  }
  cat("\n", strrep("─", 70), "\n")
  invisible(cat)
}


#' Human-readable label for an experiment row
experiment_label <- function(row) {
  parts <- c(row$condition)
  if (!is.na(row$std_dev_move) && row$std_dev_move != 0.1)
    parts <- c(parts, sprintf("σ=%.3g", row$std_dev_move))
  if (isTRUE(row$diploid))
    parts <- c(parts, "diploid")
  if (!is.na(row$fix_host))
    parts <- c(parts, sprintf("fixH=%.2g", row$fix_host))
  if (!is.na(row$fix_path))
    parts <- c(parts, sprintf("fixP=%.2g", row$fix_path))
  paste(parts, collapse = " | ")
}


#' Load simulation CSV for one catalog row
load_sim_by_row <- function(catalog, row_idx, min_gen = MIN_GEN_CUTOFF) {
  row <- catalog[row_idx, ]
  
  df <- read.csv(row$csv) %>%
    filter(event == "post", gen > min_gen) %>%
    mutate(
      fitness   = row$fitness,
      condition = row$condition,
      sigma     = row$std_dev_move,
      diploid   = row$diploid,
      label     = experiment_label(row),
      omegaPath = suppressWarnings(as.numeric(omegaPath)),
      omegaHost = suppressWarnings(as.numeric(omegaHost))
    )
  
  df
}


#' Load a filtered catalog subset into one data frame
load_sim_set <- function(catalog_subset, min_gen = MIN_GEN_CUTOFF) {
  bind_rows(
    lapply(seq_len(nrow(catalog_subset)), function(i) {
      load_sim_by_row(catalog_subset, i, min_gen)
    })
  )
}


#' Convenience: load everything in results/
load_all_discovered <- function(results_root = "results", min_gen = MIN_GEN_CUTOFF) {
  cat <- discover_experiments(results_root)
  cat("Found", nrow(cat), "experiments\n")
  load_sim_set(cat, min_gen)
}


# --- Backward-compatible scenario mapping ---
# Maps old scenario names to condition names for existing figure code
SCENARIO_TO_CONDITION <- c(
  "ET-ET"         = "EThost_ETpath",
  "ER-ER"         = "ERhost_ERpath",
  "ERpath-EThost" = "EThost_ERpath",
  "ERhost-ETpath" = "ERhost_ETpath"
)

#' Load like the old load_sim() but using auto-discovery
#' Works as drop-in replacement for existing figure functions
load_sim <- function(model, scenario, min_gen = MIN_GEN_CUTOFF,
                     sigma = NULL, diploid_filter = NULL) {
  
  cat <- discover_experiments()
  condition <- SCENARIO_TO_CONDITION[scenario]
  if (is.na(condition)) condition <- scenario
  
  subset <- cat %>% filter(fitness == model, condition == !!condition)
  
  if (nrow(subset) == 0) {
    warning(paste("No data for", model, "/", condition,
                  "\n  Available:", paste(unique(cat$condition[cat$fitness == model]),
                                          collapse = ", ")))
    return(NULL)
  }
  
  # Exclude pinned-trait runs (those belong in fig_pinned_comparison)
  subset <- subset %>% filter(is.na(fix_host), is.na(fix_path))
  
  # Apply explicit filters
  if (!is.null(sigma))
    subset <- subset %>% filter(abs(std_dev_move - sigma) < 1e-6)
  if (!is.null(diploid_filter))
    subset <- subset %>% filter(diploid == diploid_filter)
  
  # If still multiple matches, pick best: default sigma, default gamma, longest run
  if (nrow(subset) > 1) {
    subset <- subset %>%
      arrange(
        abs(std_dev_move - 0.1),                # prefer σ ≈ 0.1
        abs(ifelse(is.na(gamma), 0.01, gamma) - 0.01),  # prefer default γ ≈ 0.01
        desc(max_gens)                          # prefer longer
      ) %>%
      slice(1)
  }
  
  if (nrow(subset) == 0) {
    # Helpful message about what IS available
    available <- cat %>% filter(fitness == model, condition == !!condition)
    msg <- paste("No match for", model, "/", condition)
    if (!is.null(sigma)) msg <- paste0(msg, ", σ=", sigma)
    if (!is.null(diploid_filter)) msg <- paste0(msg, ", diploid=", diploid_filter)
    msg <- paste0(msg, "\n  Available variants:\n")
    for (i in seq_len(nrow(available))) {
      msg <- paste0(msg, "    ", experiment_label(available[i,]), "\n")
    }
    warning(msg)
    return(NULL)
  }
  
  cat("  Loading:", experiment_label(subset[1,]),
      "  (", basename(subset$dir[1]), ")\n")
  
  load_sim_by_row(subset, 1, min_gen) %>%
    mutate(scenario = scenario)  # keep old column name for figure code
}


#' Load all — backward compatible with old load_all_sims()
load_all_sims <- function(models = c("acute", "minimal"),
                          sigma = NULL, diploid_filter = NULL) {
  scenarios <- c("ET-ET", "ERpath-EThost", "ERhost-ETpath", "ER-ER")
  bind_rows(
    lapply(models, function(mod) {
      bind_rows(lapply(scenarios, function(sc) {
        load_sim(mod, sc, sigma = sigma, diploid_filter = diploid_filter)
      }))
    })
  ) %>%
    filter(nrow(.) > 0) %>%
    mutate(
      scenario = factor(scenario,
                        levels = c("ET-ET", "ERpath-EThost", "ERhost-ETpath", "ER-ER"))
    )
}


# ============================================================================
# Step-size comparison figure
# ============================================================================

#' Compare time series across step sizes for one condition
fig_step_size_comparison <- function(model_name = "acute",
                                     condition = "EThost_ETpath",
                                     diploid = NULL,
                                     width = NULL,
                                     height = NULL,
                                     filename = NULL) {
  cat <- discover_experiments() %>%
    filter(fitness == model_name, condition == !!condition,
           is.na(fix_host), is.na(fix_path)) %>%
    arrange(std_dev_move)
  if (!is.null(diploid)) {
    cat <- cat %>% filter(diploid == !!diploid)
  } else if (n_distinct(cat$diploid) > 1) {
    message("  Both diploid & haploid found for ", model_name, "/", condition,
            " — defaulting to haploid. Set diploid=TRUE to override.")
    cat <- cat %>% filter(!diploid)
  }
  
  if (nrow(cat) == 0) {
    warning("No experiments found for ", model_name, " / ", condition)
    return(NULL)
  }
  
  all_df <- load_sim_set(cat) %>%
    mutate(sigma_label = sprintf("σ = %g", sigma))
  
  # Auto-detect trait axis from model
  yax <- auto_trait_axis(model_name)
  tax <- auto_time_axis(all_df)
  
  compact_theme <- mytheme +
    theme(strip.text = element_text(size = 9),
          axis.text  = element_text(size = 9),
          axis.title = element_text(size = 11),
          panel.spacing.y = unit(2, "pt"))
  
  # Faceted time series: v and s
  p_v <- ggplot(thin_for_plot(all_df), aes(gen, v)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_wrap(~sigma_label, ncol = 1) +
    scale_x_continuous(breaks = c(1, 5e5, 1e6), labels = c("1", "500K", "1M")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(y = expression(italic(v)), x = NULL) +
    mytheme
  
  p_s <- ggplot(thin_for_plot(all_df), aes(gen, s)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_wrap(~sigma_label, ncol = 1) +
    scale_x_continuous(breaks = c(1, 5e5, 1e6), labels = c("1", "500K", "1M")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(y = expression(italic(c)), x = "Evolutionary time") +
    mytheme
  
  combined <- (p_v | p_s) +
    plot_annotation(
      title = paste0(model_name, " / ", condition, " — step size sweep")
    )
  
  if (!is.null(filename)) {
    n_sigmas <- length(unique(cat$std_dev_move))
    w <- if (!is.null(width)) width else 8
    h <- if (!is.null(height)) height else max(3, 0.8 * n_sigmas + 1)
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}


#' Compare pinned-trait experiments across conditions
#' Shows v and s time series for: base run, fixed-host, fixed-pathogen
#' @param model_name Fitness model
#' @param conditions Which conditions to show (default: all 4)
#' @param filename   Output filename (NULL = display only)
fig_pinned_comparison <- function(model_name = "acute",
                                  conditions = c("EThost_ETpath", "EThost_ERpath",
                                                 "ERhost_ETpath", "ERhost_ERpath"),
                                  sigma = 0.1,
                                  diploid = NULL,
                                  max_pts = 100,
                                  width = NULL,
                                  height = NULL,
                                  filename = NULL) {
  
  cat <- discover_experiments() %>%
    filter(fitness == model_name,
           abs(std_dev_move - sigma) < 1e-6)
  if (!is.null(diploid)) {
    cat <- cat %>% filter(diploid == !!diploid)
  } else if (n_distinct(cat$diploid) > 1) {
    message("  Both diploid & haploid found for ", model_name,
            " — defaulting to haploid. Set diploid=TRUE to override.")
    cat <- cat %>% filter(!diploid)
  }
  
  if (nrow(cat) == 0) {
    warning("No experiments found for ", model_name)
  }
  
  # Build labels for each experiment type
  label_run <- function(row) {
    if (!is.na(row$fix_host) && !is.na(row$fix_path))
      return(sprintf("fix H=%.1f, P=%.1f", row$fix_host, row$fix_path))
    if (!is.na(row$fix_host))
      return(sprintf("fix host=%.1f", row$fix_host))
    if (!is.na(row$fix_path))
      return(sprintf("fix path=%.1f", row$fix_path))
    "base"
  }
  
  # For each condition, gather base + pinned runs
  all_frames <- list()
  for (cond in conditions) {
    sub <- cat %>% filter(condition == cond)
    if (nrow(sub) == 0) next
    
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      tryCatch({
        d <- load_sim_by_row(sub, i)
        if (nrow(d) > 0) {
          td <- thin_for_plot(d, max_pts = max_pts)
          td$pin_label <- label_run(row)
          td$cond_label <- cond
          all_frames[[length(all_frames) + 1]] <- td
        }
      }, error = function(e) NULL)
    }
  }
  
  if (length(all_frames) == 0) {
    warning("No data loaded")
    return(NULL)
  }
  
  all_df <- bind_rows(all_frames)
  
  # Order: base first, then pinned
  pin_levels <- sort(unique(all_df$pin_label))
  pin_levels <- c(pin_levels[grepl("base", pin_levels)],
                  pin_levels[!grepl("base", pin_levels)])
  all_df$pin_label <- factor(all_df$pin_label, levels = pin_levels)
  
  # Nice condition labels
  cond_labels <- c(
    "EThost_ETpath" = "ET / ET",
    "EThost_ERpath" = "ET host / ER path",
    "ERhost_ETpath" = "ER host / ET path",
    "ERhost_ERpath" = "ER / ER"
  )
  all_df$cond_label <- factor(
    cond_labels[all_df$cond_label],
    levels = cond_labels[conditions]
  )
  
  yax <- auto_trait_axis(model_name)
  tax <- auto_time_axis(all_df)
  
  compact_theme <- mytheme +
    theme(strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 9),
          axis.text  = element_text(size = 11),
          axis.title = element_text(size = 14),
          panel.spacing = unit(3, "pt"))
  
  p_v <- ggplot(all_df, aes(gen, v)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_grid(pin_label ~ cond_label) +
    scale_x_continuous(breaks = c(1, 5e5, 1e6), labels = c("1", "500K", "1M")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(y = expression(italic(v)), x = NULL) +
    compact_theme
  
  p_s <- ggplot(all_df, aes(gen, s)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_grid(pin_label ~ cond_label) +
    scale_x_continuous(breaks = c(1, 5e5, 1e6), labels = c("1", "500K", "1M")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(y = expression(italic(c)), x = "Evolutionary time") +
    compact_theme
  
  n_pins <- length(pin_levels)
  n_conds <- length(unique(all_df$cond_label))
  
  combined <- (p_v / p_s) +
    plot_annotation(
      title = paste0(model_name, " — pinned trait comparison")
    )
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else max(6, 2.5 * n_conds)
    h <- if (!is.null(height)) height else max(4, 1.2 * n_pins * 2 + 1)
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("  Saved:", filename, "\n")
  }
  combined
}

# ============================================================================
# §4  FIGURE 1 -- Fitness Landscape (3-panel: host, pathogen, joint)
# ============================================================================

#' Three-panel fitness landscape for ANY registered fitness model.
fig_landscape <- function(model_name = "acute",
                          filename = NULL,
                          width = 10, height = 4) {
  
  grid    <- make_fitness_grid(model_name)
  br      <- calc_best_responses(model_name)
  nash_pt <- find_nash(model_name)
  dom     <- if (model_name %in% names(TRAIT_DISPLAY))
    TRAIT_DISPLAY[[model_name]] else TRAIT_DOMAIN[[model_name]]
  
  # Normalize to [0,1]
  grid$fH_norm <- grid$fH / max(grid$fH)
  grid$fP_norm <- grid$fP / max(grid$fP)
  grid$joint_norm <- grid$fH_norm * grid$fP_norm
  
  # Nice breaks for axis
  ax_breaks <- if (dom[2] <= 1) c(0, 0.5, 1) else pretty(dom, n = 4)
  
  # Panel A: Host fitness
  pA <- ggplot(grid, aes(x = v, y = s, z = fH_norm)) +
    geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
    geom_line(data = br$host, aes(x = v, y = s),
              color = "steelblue", linewidth = 1.5,
              inherit.aes = FALSE, linetype = "dashed") +
    labs(x = NULL, y = "c (clearance)") +
    coord_fixed(xlim = dom, ylim = dom, expand = FALSE) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    scale_fill_viridis_d(option = "viridis") +
    mytheme
  
  # Panel B: Pathogen fitness
  pB <- ggplot(grid, aes(x = v, y = s, z = fP_norm)) +
    geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
    geom_line(data = br$path, aes(x = v, y = s),
              color = "lightcoral", linewidth = 1.5,
              inherit.aes = FALSE, linetype = "dashed") +
    labs(x = "v (virulence)") +
    coord_fixed(xlim = dom, ylim = dom, expand = FALSE) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    scale_fill_viridis_d(option = "viridis") +
    mytheme
  
  # Panel C: Joint fitness with BR curves + Nash
  pC <- ggplot(grid, aes(x = v, y = s, z = joint_norm)) +
    geom_contour(aes(z = fH_norm), bins = 10, color = "steelblue", alpha = 0.5) +
    geom_contour(aes(z = fP_norm), bins = 10, color = "lightcoral", alpha = 0.5) +
    geom_vline(xintercept = nash_pt$v, linewidth = 1.5,
               color = "firebrick", linetype = "solid") +
    geom_hline(yintercept = nash_pt$s, linewidth = 1.5,
               color = "darkblue", linetype = "solid") +
    geom_line(data = br$host, aes(x = v, y = s),
              color = "steelblue", inherit.aes = FALSE,
              linewidth = 1.5, linetype = "dashed") +
    geom_line(data = br$path, aes(x = v, y = s),
              color = "lightcoral", inherit.aes = FALSE,
              linewidth = 1.5, linetype = "dashed") +
    geom_point(data = nash_pt, aes(x = v, y = s),
               color = "grey20", size = 5, inherit.aes = FALSE) +
    labs(x = NULL) +
    coord_fixed(xlim = dom, ylim = dom, expand = FALSE) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    mytheme
  
  # Assemble
  final <- (pA | strip_y(pB) | strip_y(pC)) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"),
          legend.position = "right",
          text = element_text(size = 14))
  
  if (!is.null(filename)) {
    ggsave(paste0("figures/", filename, ".pdf"), final, width = width, height = height)
    ggsave(paste0("figures/", filename, ".png"), final, width = width, height = height)
    cat("Saved:", filename, "\n")
  }
  
  cat(sprintf("  %s Nash: v*=%.3f, c*=%.3f\n",
              FITNESS_MODELS[[model_name]]$label, nash_pt$v, nash_pt$s))
  
  final
}

# Generate for any model:
# fig_landscape("acute",   "Figure1_acute")
# fig_landscape("minimal", "Figure1_minimal")
# fig_landscape("chronic", "Figure1_chronic")
# fig_landscape("taylor",  "Figure1_taylor")


# ============================================================================
# §5  FIGURE 2 -- Strategy lines, slopes, and destabilization
# ============================================================================

fig_strategy_panels <- function(model_name = "acute",
                                width = NULL, height = NULL,
                                filename = "Figure2") {
  grid    <- make_fitness_grid(model_name, resolution = 300)
  br      <- calc_best_responses(model_name)
  nash_pt <- find_nash(model_name)
  dom     <- TRAIT_DOMAIN[[model_name]]
  v_star  <- nash_pt$v
  s_star  <- nash_pt$s
  
  col_host <- "darkblue"
  col_path <- "firebrick"
  
  ax_breaks <- if (dom[2] <= 1) c(0, 0.5, 1) else pretty(dom, n = 4)
  
  host_line_fn <- function(v, bS, mS) clamp_trait(bS + mS * v, model_name)
  path_line_fn <- function(s, bV, mV) clamp_trait(bV + mV * s, model_name)
  
  pad <- (dom[2] - dom[1]) * 0.1
  v_seq <- seq(dom[1] - pad, dom[2] + pad, length.out = 500)
  s_seq <- seq(dom[1] - pad, dom[2] + pad, length.out = 500)
  
  # Panel builder
  make_panel <- function(bV, mV, bS, mS,
                         show_stable = TRUE,
                         boundary_pts = NULL,
                         show_yaxis = TRUE) {
    
    host_data <- data.frame(
      v = v_seq, s = host_line_fn(v_seq, bS, mS)
    ) %>% filter(v >= dom[1], v <= dom[2], s >= dom[1], s <= dom[2])
    path_data <- data.frame(
      s = s_seq, v = path_line_fn(s_seq, bV, mV)
    ) %>% filter(v >= dom[1], v <= dom[2], s >= dom[1], s <= dom[2])
    
    # Interior intersection
    den <- 1 - mV * mS
    v_int <- if (abs(den) > 1e-9) (bV + mV * bS) / den else v_star
    s_int <- bS + mS * v_int
    
    p <- ggplot(grid, aes(v, s)) +
      geom_contour(aes(z = fP), color = "lightcoral", bins = 10,
                   linewidth = 0.3, alpha = 0.7) +
      geom_contour(aes(z = fH), color = "steelblue", bins = 10,
                   linewidth = 0.3, alpha = 0.7) +
      geom_line(data = br$host, aes(v, s),
                color = "lightcoral", linewidth = 2, linetype = "dashed") +
      geom_line(data = br$path, aes(v, s),
                color = "steelblue", linewidth = 2, linetype = "dashed") +
      geom_line(data = host_data, aes(v, s),
                linetype = "solid", linewidth = 1.5, color = col_host) +
      geom_line(data = path_data, aes(v, s),
                linetype = "solid", linewidth = 1.5, color = col_path) +
      coord_fixed(xlim = dom, ylim = dom) +
      scale_x_continuous(breaks = ax_breaks) +
      scale_y_continuous(breaks = ax_breaks) +
      mytheme
    
    if (show_yaxis) {
      p <- p + labs(x = "v (virulence)", y = "c (clearance)")
    } else {
      p <- p + labs(x = NULL, y = NULL) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    if (show_stable) {
      p <- p + geom_point(aes(x = v_int, y = s_int), size = 5, colour = "black")
    } else {
      p <- p + geom_point(aes(x = v_int, y = s_int), size = 5, shape = 21,
                          fill = "gray70", colour = "black", stroke = 1)
    }
    
    if (!is.null(boundary_pts)) {
      p <- p + geom_point(data = boundary_pts, aes(x = v, y = s),
                          size = 5, colour = "black")
    }
    p
  }
  
  # Panel A: ET (flat strategies at Nash)
  pA <- make_panel(v_star, 0, s_star, 0, TRUE, NULL, TRUE)
  
  # Panel B: Host gains slope
  mS_B <- 0.5
  pB <- make_panel(v_star, 0, s_star - mS_B * v_star, mS_B, TRUE, NULL, FALSE)
  
  # Panel C: Pathogen gains slope
  mV_C <- 0.8
  pC <- make_panel(v_star - mV_C * s_star, mV_C, s_star, 0, TRUE, NULL, FALSE)
  
  # Panel D: Both large slopes -- destabilization
  mS_D <- 2.5; mV_D <- 2.5
  bS_D <- s_star - mS_D * v_star
  bV_D <- v_star - mV_D * s_star
  if (bS_D > dom[1] - 1e-3) bS_D <- dom[1] - 1e-3
  if (bV_D > dom[1] - 1e-3) bV_D <- dom[1] - 1e-3
  
  boundary_D <- data.frame(
    v = c(clamp_trait(bV_D, model_name), clamp_trait(bV_D + mV_D * dom[2], model_name)),
    s = c(clamp_trait(bS_D, model_name), clamp_trait(bS_D + mS_D * dom[2], model_name))
  )
  pD <- make_panel(bV_D, mV_D, bS_D, mS_D, FALSE, boundary_D, FALSE)
  
  # Combine — only panel A has axis titles
  final <- (pA | pB | pC | pD) +
    plot_annotation(tag_levels = "A") +
    plot_layout(widths = rep(1, 4)) &
    theme(plot.tag = element_text(face = "bold", size = 16),
          plot.tag.position = c(-0.01, 0.76))
  
  w <- if (!is.null(width)) width else 9
  h <- if (!is.null(height)) height else 6
  ggsave(paste0("figures/", filename, ".pdf"), final, width = w, height = h)
  ggsave(paste0("figures/", filename, ".png"), final, width = w, height = h)
  
  cat("\nStability conditions:\n")
  cat(sprintf("  A: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (stable)\n", 0, 0, 0))
  cat(sprintf("  B: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (stable)\n", mS_B, 0, 0))
  cat(sprintf("  C: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (stable)\n", 0, mV_C, 0))
  cat(sprintf("  D: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (UNSTABLE)\n",
              mS_D, mV_D, abs(mS_D * mV_D)))
  final
}


# ============================================================================
# §6  FIGURE 3 -- Time Series (v, c, W_H, W_P, omega)
# ============================================================================

# Shared axis settings — DEFAULTS for 1M-gen runs.
# Auto-overridden by auto_time_axis() when data is passed.
X_LIMS_LIN   <- c(1e4, 1e6)
X_BREAKS_LIN <- c(10000, 505000, 1e6)
X_LABS_LIN   <- c("10K", "505K", "1M")
W_LIMS_LOG   <- c(1e-2, 1e6)
W_BREAKS_LOG <- c(1e-2, 1e2, 1e6)
W_LABS_LOG   <- trans_format("log10", math_format(10^.x))


#' Compute sensible time-axis settings from data
#' Returns list(lims, breaks, labels) that can be passed to line_panel/omega_panel
auto_time_axis <- function(df, n_breaks = 3) {
  gen_range <- range(df$gen, na.rm = TRUE)
  lo <- gen_range[1]
  hi <- gen_range[2]
  
  # Nice labels
  fmt_label <- function(x) {
    if (x >= 1e6) sprintf("%.0fM", x / 1e6)
    else if (x >= 1e3) sprintf("%.0fK", x / 1e3)
    else as.character(x)
  }
  
  # Use clean breaks for common run lengths
  if (hi >= 9e5 && hi <= 1.1e6) {
    brk <- c(1, 5e5, 1e6)
    brk <- brk[brk >= lo & brk <= hi * 1.01]
  } else if (hi >= 4.5e5 && hi < 9e5) {
    brk <- c(1, 2.5e5, 5e5)
    brk <- brk[brk >= lo & brk <= hi * 1.01]
  } else {
    brk <- pretty(c(lo, hi), n = n_breaks)
    brk <- brk[brk >= lo & brk <= hi]
    if (length(brk) == 0) brk <- c(lo, hi)
    if (length(brk) > n_breaks + 1) {
      idx <- round(seq(1, length(brk), length.out = n_breaks + 1))
      brk <- brk[idx]
    }
  }
  
  labs <- sapply(brk, fmt_label)
  list(lims = c(lo, hi), breaks = brk, labels = labs)
}


#' Compute trait-axis limits from model name (or from data if model unknown)
auto_trait_axis <- function(model_name = NULL, df = NULL, y_var = NULL) {
  # Try model-specific display domain first
  if (!is.null(model_name) && model_name %in% names(TRAIT_DISPLAY)) {
    dom <- TRAIT_DISPLAY[[model_name]]
    # Clean breaks: 0, 0.5, 1 for [0,1] models; pretty() for wider domains (taylor)
    brk <- if (dom[2] <= 1) c(0, 0.5, 1) else pretty(dom, n = 4)
    return(list(lims = dom, breaks = brk))
  }
  # Fall back to data range
  if (!is.null(df) && !is.null(y_var) && y_var %in% names(df)) {
    rng <- range(df[[y_var]], na.rm = TRUE)
    pad <- (rng[2] - rng[1]) * 0.05
    dom <- c(max(0, rng[1] - pad), rng[2] + pad)
    brk <- pretty(dom, n = 4)
    return(list(lims = dom, breaks = brk))
  }
  # Default
  list(lims = c(0, 1), breaks = c(0, 0.5, 1))
}


# --- Line panel (v, c, W) ---
# model_name: if provided, uses TRAIT_DOMAIN for y-limits on trait variables
# x_lims/x_breaks/x_labels: if NULL, auto-detected from data
# use_step: if TRUE, uses geom_step instead of geom_line (better for SSWM data)
line_panel <- function(df, y_var, ylab = NULL,
                       show_xlab = FALSE, show_ylab = TRUE,
                       model_name = NULL,
                       x_lims = NULL, x_breaks = NULL, x_labels = NULL,
                       use_step = FALSE) {
  
  # Auto-detect time axis from data if not specified
  if (is.null(x_lims)) {
    tax <- auto_time_axis(df)
    x_lims <- tax$lims; x_breaks <- tax$breaks; x_labels <- tax$labels
  }
  
  # Auto-detect trait axis: use model domain for v/s, [0,1] for fitness
  is_trait <- y_var %in% c("v", "s")
  if (is_trait) {
    yax <- auto_trait_axis(model_name, df, y_var)
    y_lims <- yax$lims; y_breaks <- yax$breaks
  } else {
    y_lims <- c(0, 1); y_breaks <- c(0, 0.5, 1)
  }
  
  geom_fn <- if (use_step) geom_step else geom_line
  
  p <- ggplot(df, aes(x = gen, y = .data[[y_var]])) +
    geom_fn(linewidth = 0.5, alpha = 0.85) +
    scale_x_continuous(limits = x_lims, breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = y_lims, breaks = y_breaks) +
    coord_cartesian(xlim = x_lims) +
    mytheme
  
  if (show_ylab && !is.null(ylab)) p <- p + labs(y = ylab)
  else p <- p + labs(y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  if (show_xlab) p <- p + labs(x = "Evolutionary time")
  else p <- p + labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

# --- Omega spike panel ---
omega_panel <- function(df, who = c("Path", "Host"),
                        cap = 1e6, bins = 2000,
                        show_xlab = FALSE, show_ylab = TRUE, ylab = NULL,
                        x_lims = NULL, x_breaks = NULL, x_labels = NULL) {
  who <- match.arg(who)
  omega_col <- if (who == "Path") "omegaPath" else "omegaHost"
  
  # Auto-detect time axis from data if not specified
  if (is.null(x_lims)) {
    tax <- auto_time_axis(df)
    x_lims <- tax$lims; x_breaks <- tax$breaks; x_labels <- tax$labels
  }
  
  rng <- range(x_lims)
  edges <- seq(rng[1], rng[2], length.out = bins + 1)
  
  dat <- df %>%
    filter(event == "post", gen >= rng[1], gen <= rng[2]) %>%
    transmute(
      x = gen,
      y = pmin(cap, pmax(1e-12, as.numeric(.data[[omega_col]]))),
      lbin = cut(gen, breaks = edges, include.lowest = TRUE)
    ) %>%
    group_by(lbin) %>%
    slice_max(order_by = y, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  p <- ggplot(dat) +
    geom_segment(aes(x = x, xend = x, y = 1e-2, yend = pmax(1e-2, y)),
                 linewidth = 0.35, alpha = 0.9) +
    scale_x_continuous(limits = x_lims, breaks = x_breaks, labels = x_labels) +
    scale_y_log10(limits = W_LIMS_LOG, breaks = W_BREAKS_LOG,
                  labels = W_LABS_LOG, minor_breaks = NULL) +
    coord_cartesian(xlim = x_lims, ylim = W_LIMS_LOG) +
    mytheme
  
  if (show_ylab && !is.null(ylab)) p <- p + labs(y = ylab)
  else p <- p + labs(y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  if (show_xlab) p <- p + labs(x = "Evolutionary time")
  else p <- p + labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

#' Build the full 6x4 time series figure for one fitness model
#' Columns: ET-ET | ERpath-EThost | ERhost-ETpath | ER-ER
#' Rows: v, c, W_P, W_H, omega_P, omega_H
#' @param sigma  Filter by step size (e.g. 0.01). NULL = default/first match.
#' @param diploid  Filter by diploid flag. NULL = any.
#' @param max_pts  Max points per panel after thinning (default 2000). Use Inf for no thinning.
#' @param smooth   Rolling average window size (in number of points). NULL = no smoothing.
#'                 Try smooth = 50 for gentle smoothing, 200 for heavy.
#' @param step     If TRUE, uses geom_step (flat between events, vertical jumps).
#'                 More accurate for SSWM data and looks better when thinned.
#' @param x_lims   Override time axis limits, e.g. c(0, 1e6). NULL = auto-detect.
#' @param x_breaks Override time axis breaks, e.g. c(0, 5e5, 1e6). NULL = auto-detect.
#' @param x_labels Override time axis labels, e.g. c("0", "500K", "1M"). NULL = auto-detect.
fig_timeseries <- function(model_name = "acute", filename = NULL,
                           sigma = NULL, diploid = NULL,
                           max_pts = 2000, smooth = NULL,
                           step = FALSE,
                           width = NULL, height = NULL,
                           x_lims = NULL, x_breaks = NULL, x_labels = NULL) {
  
  scenarios <- c("ET-ET", "ERpath-EThost", "ERhost-ETpath", "ER-ER")
  col_titles <- c("ET / ET", "ET host / ER path",
                  "ER host / ET path", "ER / ER")
  
  tag <- model_name
  if (!is.null(diploid) && diploid) tag <- paste0(tag, " (diploid)")
  if (!is.null(sigma)) tag <- paste0(tag, " σ=", sigma)
  cat("\n  Loading time series for:", tag, "\n")
  
  dfs <- setNames(
    lapply(scenarios, function(sc) {
      d <- load_sim(model_name, sc, sigma = sigma, diploid_filter = diploid)
      if (is.null(d) || nrow(d) == 0) return(NULL)
      td <- thin_for_plot(d, max_pts = max_pts)
      if (!is.null(smooth) && smooth > 1) {
        k <- min(smooth, nrow(td))
        for (col in c("v", "s", "pathFit", "hostFit")) {
          if (col %in% names(td))
            td[[col]] <- rollmean(td[[col]], k = k, fill = NA, align = "center")
        }
        td <- td %>% filter(!is.na(v))
      }
      td
    }),
    scenarios
  )
  
  # Which scenarios actually loaded?
  available <- !vapply(dfs, is.null, logical(1))
  if (sum(available) == 0) {
    warning("No data loaded for any condition")
    return(NULL)
  }
  
  # Use only available scenarios
  active_scenarios <- scenarios[available]
  active_titles    <- col_titles[available]
  active_dfs       <- dfs[available]
  n_cols           <- length(active_scenarios)
  
  if (sum(available) < 4) {
    cat("  Note: only", sum(available), "of 4 conditions available:",
        paste(active_scenarios, collapse = ", "), "\n")
  }
  
  # Auto-detect shared time axis from all data, or use overrides
  all_gens <- unlist(lapply(active_dfs, function(d) d$gen))
  tax <- auto_time_axis(data.frame(gen = all_gens))
  if (!is.null(x_lims))   tax$lims   <- x_lims
  if (!is.null(x_breaks)) tax$breaks <- x_breaks
  if (!is.null(x_labels)) tax$labels <- x_labels
  
  rows <- list(
    list(var = "v",       ylab = expression(italic(v))),
    list(var = "s",       ylab = expression(italic(c))),
    list(var = "pathFit", ylab = expression(W[P])),
    list(var = "hostFit", ylab = expression(W[H]))
  )
  
  panels <- list()
  
  # Trait/fitness rows
  for (ri in seq_along(rows)) {
    row <- rows[[ri]]
    for (ci in seq_along(active_scenarios)) {
      p <- line_panel(
        active_dfs[[ci]], row$var,
        ylab = row$ylab,
        show_ylab = (ci == 1),
        show_xlab = FALSE,
        model_name = model_name,
        x_lims = tax$lims, x_breaks = tax$breaks, x_labels = tax$labels,
        use_step = step
      )
      # Column title on first row
      if (ri == 1) {
        p <- p + labs(title = active_titles[ci]) +
          theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
      }
      panels[[length(panels) + 1]] <- p
    }
  }
  
  # Omega rows
  for (who in c("Path", "Host")) {
    ylab <- if (who == "Path") expression(omega[P]) else expression(omega[H])
    is_last <- (who == "Host")
    for (ci in seq_along(active_scenarios)) {
      panels[[length(panels) + 1]] <- omega_panel(
        active_dfs[[ci]], who,
        ylab = ylab,
        show_ylab = (ci == 1),
        show_xlab = is_last,
        x_lims = tax$lims, x_breaks = tax$breaks, x_labels = tax$labels
      )
    }
  }
  
  total <- wrap_plots(panels, ncol = n_cols, byrow = TRUE) +
    plot_annotation(
      title = tag,
      tag_levels = "A"
    ) &
    theme(plot.tag.position = "topleft",
          plot.tag = element_text(face = "bold", size = 12))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 2.5 * n_cols + 1
    h <- if (!is.null(height)) height else 9
    ggsave(paste0("figures/", filename, ".pdf"), total, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), total, width = w, height = h)
    cat("  Saved:", filename, "\n")
  }
  total
}


# ============================================================================
# §7  FIGURE 4 -- Phase-space hex density on fitness landscape
# ============================================================================

hex_landscape_panel <- function(data, model_name = "acute",
                                nbins = 100, count_limits = NULL,
                                show_x = TRUE, show_y = TRUE,
                                show_nash = TRUE) {
  mod <- FITNESS_MODELS[[model_name]]
  dom <- TRAIT_DOMAIN[[model_name]]
  
  fgrid <- expand.grid(
    v = seq(dom[1], dom[2], length.out = 150),
    s = seq(dom[1], dom[2], length.out = 150)
  ) %>% mutate(fH = mod$fH(v, s), fP = mod$fP(v, s))
  
  br <- calc_best_responses(model_name, n = 300)
  nash_pt <- if (show_nash) find_nash(model_name) else NULL
  
  ax_breaks <- if (dom[2] <= 1) c(0, 0.5, 1) else pretty(dom, n = 4)
  
  p <- ggplot() +
    geom_contour(data = fgrid, aes(v, s, z = fP),
                 color = "lightcoral", alpha = 0.3, bins = 12, linewidth = 0.5) +
    geom_contour(data = fgrid, aes(v, s, z = fH),
                 color = "steelblue", alpha = 0.3, bins = 12, linewidth = 0.5) +
    geom_line(data = br$host, aes(v, s),
              color = "steelblue", linewidth = 1, linetype = "dashed") +
    geom_line(data = br$path, aes(v, s),
              color = "lightcoral", linewidth = 1, linetype = "dashed") +
    geom_hex(data = data, aes(x = v, y = s), bins = nbins, alpha = 0.7)
  
  fill_args <- list(option = "plasma", name = "Count",
                    trans = "log10", oob = scales::squish)
  if (!is.null(count_limits)) fill_args$limits <- count_limits
  p <- p + do.call(scale_fill_viridis_c, fill_args)
  
  if (!is.null(nash_pt)) {
    p <- p + geom_point(data = nash_pt, aes(x = v, y = s),
                        size = 4, color = "black", shape = 16)
  }
  
  # Realized mean as yellow cross
  mean_pt <- data.frame(v = mean(data$v, na.rm = TRUE),
                        s = mean(data$s, na.rm = TRUE))
  p <- p + geom_point(data = mean_pt, aes(x = v, y = s),
                      size = 4, color = "#FFD700", shape = 4, stroke = 1.5)
  
  p <- p +
    scale_x_continuous(breaks = ax_breaks, limits = dom) +
    scale_y_continuous(breaks = ax_breaks, limits = dom) +
    coord_fixed() + mytheme + theme(legend.position = "none")
  
  if (!show_x) p <- strip_x(p) else p <- p + labs(x = "v (virulence)")
  if (!show_y) p <- strip_y(p) else p <- p + labs(y = "c (clearance)")
  p
}

#' Full hex-density figure: models × scenarios grid
#' @param all_data Combined data from load_all_sims()
#' @param models Character vector of model names to include
#' @param filename Output file (or NULL for display only)
fig_hex_combined <- function(all_data = NULL, models = c("acute", "minimal"),
                             filename = "Figure4_hex",
                             nbins = 100,
                             sigma = NULL, diploid = NULL,
                             width = NULL, height = NULL) {
  
  # Load data if not provided
  if (is.null(all_data)) {
    all_data <- load_all_sims(models, sigma = sigma, diploid_filter = diploid)
  }
  
  scenarios <- c("ET-ET", "ERpath-EThost", "ERhost-ETpath", "ER-ER")
  col_titles <- c("ET / ET", "ET host / ER path",
                  "ER host / ET path", "ER / ER")
  
  panels <- list()
  n_mod <- length(models)
  
  for (mi in seq_along(models)) {
    mod <- models[mi]
    is_bottom <- (mi == n_mod)
    for (ci in seq_along(scenarios)) {
      sc <- scenarios[ci]
      d <- all_data %>% filter(fitness == mod, scenario == sc)
      if (nrow(d) == 0) {
        panels[[length(panels) + 1]] <- ggplot() + theme_void()
        next
      }
      
      p <- hex_landscape_panel(d, model_name = mod, nbins = nbins,
                               show_x = is_bottom,
                               show_y = (ci == 1))
      
      # Add column title on top row
      if (mi == 1) {
        p <- p + labs(title = col_titles[ci]) +
          theme(plot.title = element_text(hjust = 0.5, size = 11,
                                          face = "bold"))
      }
      # Add row label on left column
      if (ci == 1) {
        mod_label <- FITNESS_MODELS[[mod]]$label
        p <- p + labs(y = paste0(mod_label, "\nc (clearance)"))
      }
      
      panels[[length(panels) + 1]] <- p
    }
  }
  
  n_sc <- length(scenarios)
  combined <- wrap_plots(panels, ncol = n_sc, byrow = TRUE) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 12))
  
  if (!is.null(filename)) {
    h <- 3.5 * n_mod + 0.5
    w_out <- if (!is.null(width)) width else 3.5 * n_sc
    h_out <- if (!is.null(height)) height else h
    ggsave(paste0("figures/", filename, ".pdf"), combined,
           width = w_out, height = h_out)
    ggsave(paste0("figures/", filename, ".png"), combined,
           width = w_out, height = h_out)
    cat("Saved:", filename, "\n")
  }
  combined
}

# ============================================================================
# §8  FIGURE 5 -- Strategy snapshots
# ============================================================================

make_snapshot_panel <- function(es_data, gen_num, grid,
                                title_label = "",
                                lag_gens = 100,
                                show_x = TRUE, show_y = TRUE,
                                show_x_label = FALSE, show_y_label = FALSE,
                                model_name = "acute") {
  
  snapshot <- es_data %>% filter(gen == gen_num)
  if (nrow(snapshot) == 0) return(NULL)
  
  prev <- es_data %>%
    filter(gen <= gen_num - lag_gens) %>%
    arrange(desc(gen)) %>% slice(1)
  if (nrow(prev) == 0) prev <- snapshot
  
  dom <- TRAIT_DOMAIN[[model_name]]
  pad <- (dom[2] - dom[1]) * 0.05
  lo <- dom[1] - pad
  hi <- dom[2] + pad
  ax_breaks <- if (dom[2] <= 1) c(0, 0.5, 1) else pretty(dom, n = 4)
  
  v_seq <- seq(lo, hi, length.out = 500)
  s_seq <- seq(lo, hi, length.out = 500)
  clip <- function(df) df %>% filter(v >= dom[1], v <= dom[2],
                                     s >= dom[1], s <= dom[2])
  
  host_now  <- clip(tibble(v = v_seq, s = snapshot$bS + snapshot$mS * v_seq))
  path_now  <- clip(tibble(s = s_seq, v = snapshot$bV + snapshot$mV * s_seq))
  host_prev <- clip(tibble(v = v_seq, s = prev$bS + prev$mS * v_seq))
  path_prev <- clip(tibble(s = s_seq, v = prev$bV + prev$mV * s_seq))
  
  p <- ggplot() +
    geom_contour(data = grid, aes(v, s, z = fP),
                 color = "lightcoral", bins = 10, linewidth = 0.3, alpha = 0.7) +
    geom_contour(data = grid, aes(v, s, z = fH),
                 color = "steelblue", bins = 10, linewidth = 0.3, alpha = 0.7) +
    geom_line(data = host_prev, aes(v, s),
              color = "darkblue", linetype = "dotted", linewidth = 1.5, alpha = 0.85) +
    geom_line(data = path_prev, aes(v, s),
              color = "firebrick", linetype = "dotted", linewidth = 1.5, alpha = 0.85) +
    geom_line(data = host_now, aes(v, s),
              color = "darkblue", linetype = "solid", linewidth = 1.5) +
    geom_line(data = path_now, aes(v, s),
              color = "firebrick", linetype = "solid", linewidth = 1.5) +
    geom_point(aes(x = snapshot$v, y = snapshot$s), size = 5, color = "black") +
    coord_fixed(xlim = dom, ylim = dom) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    labs(title = title_label, x = NULL, y = NULL) +
    mytheme +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.02, vjust = -1, size = 14))
  
  if (!show_x)
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  
  if (!show_y)
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  
  if (show_x_label) p <- p + labs(x = "v (virulence)")
  if (show_y_label) p <- p + labs(y = "c (clearance)")
  
  p
}


#' Full strategy-snapshot figure: pick N evenly-spaced generations from an
#' ER-ER run, show how host & pathogen strategy lines evolve.
#' @param es_data  Data frame from an ES/ER-ER run (must have bS, mS, bV, mV)
#' @param model_name  Fitness model for background contours
#' @param n_panels  How many snapshots (default 6, arranged in 2 rows)
#' @param gens  Optional: explicit generation numbers to snapshot
#' @param filename  Output filename (or NULL for display only)
fig_snapshots <- function(es_data = NULL, model_name = "acute",
                          n_panels = 6, gens = NULL,
                          condition = "ERhost_ERpath", sigma = 0.1,
                          diploid = NULL,
                          width = NULL, height = NULL,
                          filename = "Figure5_snapshots") {
  
  if (is.null(es_data)) {
    es_data <- load_sim(model_name, condition, sigma = sigma,
                        diploid_filter = diploid)
    if (is.null(es_data) || nrow(es_data) == 0) {
      warning("No data found for ", model_name, "/", condition)
      return(invisible(NULL))
    }
  }
  
  grid <- make_fitness_grid(model_name, resolution = 200)
  
  # Pick snapshot generations
  available_gens <- sort(unique(es_data$gen))
  if (is.null(gens)) {
    idx <- round(seq(1, length(available_gens), length.out = n_panels))
    gens <- available_gens[idx]
  }
  n_panels <- length(gens)
  
  ncol <- min(n_panels, 3)
  nrow <- ceiling(n_panels / ncol)
  mid_col <- ceiling(ncol / 2)    # middle column for x-label
  mid_row <- ceiling(nrow / 2)    # middle row for y-label
  
  panels <- list()
  for (i in seq_along(gens)) {
    g <- gens[i]
    ri <- ceiling(i / ncol)
    ci <- ((i - 1) %% ncol) + 1
    
    panels[[i]] <- make_snapshot_panel(
      es_data, gen_num = g, grid = grid,
      title_label = paste0("gen ", format(g, big.mark = ",")),
      show_x = (ri == nrow),
      show_y = (ci == 1),
      show_x_label = (ri == nrow && ci == mid_col),
      show_y_label = (ci == 1 && ri == mid_row),
      model_name = model_name
    )
  }
  
  # Drop NULLs (if a generation wasn't found)
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) {
    warning("No panels could be created — check generation numbers")
    return(invisible(NULL))
  }
  
  combined <- wrap_plots(panels, ncol = ncol) +
    plot_annotation(
      title = paste0("Strategy snapshots (", model_name, ")"),
      tag_levels = "A"
    ) &
    theme(plot.tag = element_text(face = "bold", size = 14))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 3.5 * ncol
    h <- if (!is.null(height)) height else 3.8 * nrow
    ggsave(paste0("figures/", filename, ".pdf"), combined,
           width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined,
           width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}


# ============================================================================
# §9  FIGURE 6 -- Nash region violations & occupancy
# ============================================================================

calc_violation_grid <- function(model_name = "acute", resolution = 80) {
  # -----------------------------------------------------------------------
  # Nash stability via BEST-RESPONSE SLOPES (second derivatives)
  #
  # Host best response c*(v) satisfies  dW_H/ds = 0.
  #   Slope:  dc*/dv = -W_H,sv / W_H,ss   (implicit function theorem)
  #
  # Pathogen best response v*(c) satisfies  dW_P/dv = 0.
  #   Slope:  dv*/dc = -W_P,vs / W_P,vv
  #
  # Stability product:  dc*/dv * dv*/dc
  #   |product| < 1  =>  compatible (stable Nash)
  #   product  >= 1  =>  violates (red)   — both slopes same sign, too steep
  #   product  <= -1 =>  violates (blue)  — slopes opposite sign, too steep
  # -----------------------------------------------------------------------
  mod <- FITNESS_MODELS[[model_name]]
  dom <- TRAIT_DOMAIN[[model_name]]
  h <- (dom[2] - dom[1]) * 1e-3   # finite-difference step
  lo <- dom[1] + h * 2
  hi <- dom[2] - h * 2
  
  grid <- expand.grid(
    v = seq(lo, hi, length.out = resolution),
    s = seq(lo, hi, length.out = resolution)
  )
  
  # Vectorised second-derivative helpers (central differences)
  # W_ss  = d²W/ds²     = [W(v, s+h) - 2W(v, s) + W(v, s-h)] / h²
  # W_vv  = d²W/dv²     = [W(v+h, s) - 2W(v, s) + W(v-h, s)] / h²
  # W_sv  = d²W/(ds dv)  = [W(v+h,s+h) - W(v+h,s-h) - W(v-h,s+h) + W(v-h,s-h)] / (4h²)
  
  grid %>%
    rowwise() %>%
    mutate(
      # --- Host second partials (needed: W_H,ss and W_H,sv) ---
      fH_ss = (mod$fH(v, min(s + h, hi)) - 2 * mod$fH(v, s) +
                 mod$fH(v, max(s - h, lo))) / h^2,
      fH_sv = (mod$fH(min(v + h, hi), min(s + h, hi)) -
                 mod$fH(min(v + h, hi), max(s - h, lo)) -
                 mod$fH(max(v - h, lo), min(s + h, hi)) +
                 mod$fH(max(v - h, lo), max(s - h, lo))) / (4 * h^2),
      
      # --- Pathogen second partials (needed: W_P,vv and W_P,vs) ---
      fP_vv = (mod$fP(min(v + h, hi), s) - 2 * mod$fP(v, s) +
                 mod$fP(max(v - h, lo), s)) / h^2,
      fP_vs = (mod$fP(min(v + h, hi), min(s + h, hi)) -
                 mod$fP(min(v + h, hi), max(s - h, lo)) -
                 mod$fP(max(v - h, lo), min(s + h, hi)) +
                 mod$fP(max(v - h, lo), max(s - h, lo))) / (4 * h^2),
      
      # --- Best-response slopes ---
      # Host:    dc*/dv = -W_H,sv / W_H,ss
      # Pathogen: dv*/dc = -W_P,vs / W_P,vv
      br_host = -fH_sv / (fH_ss + 1e-12),   # dc*/dv
      br_path = -fP_vs / (fP_vv + 1e-12),   # dv*/dc
      
      # --- Stability product ---
      prod_mv = br_host * br_path,
      zone = case_when(
        prod_mv >= 1  ~ "violates (>=1)",
        prod_mv <= -1 ~ "violates (<=-1)",
        TRUE          ~ "compatible"
      )
    ) %>%
    ungroup()
}


#' Figure 6A: Nash-stability violation map — shows where |mS·mV| > 1 in
#' trait space, overlaid with simulation trajectory density.
fig_nash_violation_map <- function(model_name = "acute",
                                   es_data = NULL,
                                   condition = "ERhost_ERpath", sigma = 0.1,
                                   diploid = NULL,
                                   resolution = 80,
                                   width = NULL, height = NULL,
                                   filename = NULL) {
  if (is.null(es_data)) {
    es_data <- load_sim(model_name, condition, sigma = sigma,
                        diploid_filter = diploid)
  }
  
  vgrid <- calc_violation_grid(model_name, resolution)
  nash_pt <- find_nash(model_name)
  dom <- TRAIT_DOMAIN[[model_name]]
  
  p <- ggplot(vgrid, aes(v, s)) +
    geom_tile(aes(fill = zone), alpha = 0.7) +
    scale_fill_manual(values = c("compatible" = "#E8E8E8",
                                 "violates (>=1)" = "#FFAAAA",
                                 "violates (<=-1)" = "#AAD4FF"),
                      name = "Stability") +
    geom_point(data = nash_pt, aes(x = v, y = s),
               size = 5, color = "black", shape = 16) +
    coord_fixed(xlim = dom, ylim = dom) +
    labs(x = "v (virulence)", y = "c (clearance)",
         title = paste0(model_name, " — Nash stability regions")) +
    mytheme
  
  # Overlay simulation trajectory if provided
  if (!is.null(es_data) && nrow(es_data) > 0) {
    thin <- thin_for_plot(es_data)
    p <- p + geom_path(data = thin, aes(x = v, y = s),
                       color = "grey40", alpha = 0.15, linewidth = 0.2)
  }
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 6
    h <- if (!is.null(height)) height else 5.5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


#' Figure 6B: Slope distribution — scatter of (mS, mV) from ER-ER runs,
#' with stability hyperbolas at mS·mV = ±1.
#' This maps to manuscript "Figure 2: Phase Space: Strategy Slopes."
fig_slope_distribution <- function(es_data = NULL,
                                   model_name = "acute",
                                   condition = "ERhost_ERpath", sigma = 0.1,
                                   diploid = NULL,
                                   filename = NULL,
                                   width = NULL, height = NULL,
                                   title_label = "Phase Space: Strategy Slopes") {
  
  if (is.null(es_data)) {
    es_data <- load_sim(model_name, condition, sigma = sigma,
                        diploid_filter = diploid)
    if (is.null(es_data) || nrow(es_data) == 0) {
      warning("No data found for ", model_name, "/", condition)
      return(invisible(NULL))
    }
  }
  
  if (!"mS" %in% names(es_data) || !"mV" %in% names(es_data)) {
    warning("Data must include mS and mV columns (ER run)")
    return(invisible(NULL))
  }
  
  thin <- thin_for_plot(es_data)
  
  # Stability hyperbolas: mS * mV = ±1
  ms_seq <- seq(-5, 5, length.out = 500)
  hyp_pos <- data.frame(mS = ms_seq, mV =  1 / ms_seq)
  hyp_neg <- data.frame(mS = ms_seq, mV = -1 / ms_seq)
  
  # Compute axis limits from data (clip extreme outliers with quantiles)
  q_mS <- quantile(thin$mS, c(0.01, 0.99), na.rm = TRUE)
  q_mV <- quantile(thin$mV, c(0.01, 0.99), na.rm = TRUE)
  pad <- 0.15  # 15% padding
  xlim <- q_mS + c(-1, 1) * diff(q_mS) * pad
  ylim <- q_mV + c(-1, 1) * diff(q_mV) * pad
  
  # Classify interior vs boundary
  thin <- thin %>%
    mutate(
      prod_slopes = mS * mV,
      interior = abs(prod_slopes) < 1
    )
  
  pct_interior <- mean(thin$interior, na.rm = TRUE) * 100
  
  # Plot unstable first, then stable on top so blue is visible
  p <- ggplot(thin, aes(x = mS, y = mV)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray70", linewidth = 0.3) +
    geom_point(data = thin %>% filter(!interior),
               aes(color = interior), size = 0.8, alpha = 0.4) +
    geom_point(data = thin %>% filter(interior),
               aes(color = interior), size = 0.8, alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "#2171B5", "FALSE" = "#CB181D"),
                       labels = c("TRUE" = "stable", "FALSE" = "unstable"),
                       name = "Region") +
    geom_line(data = hyp_pos %>% filter(abs(mV) < max(abs(ylim))),
              aes(mS, mV), color = "red", linetype = "dashed",
              linewidth = 0.8, inherit.aes = FALSE) +
    geom_line(data = hyp_neg %>% filter(abs(mV) < max(abs(ylim))),
              aes(mS, mV), color = "red", linetype = "dashed",
              linewidth = 0.8, inherit.aes = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = expression(m[S]~"(host slope)"),
         y = expression(m[V]~"(pathogen slope)"),
         title = title_label,
         subtitle = sprintf("%.0f%% of time in stable region", pct_interior)) +
    mytheme +
    theme(legend.position = "right")
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 7
    h <- if (!is.null(height)) height else 6
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


#' Combined Figure 6: Nash violations + slope distribution
fig_nash_combined <- function(model_name = "acute",
                              es_data = NULL,
                              condition = "ERhost_ERpath", sigma = 0.1,
                              diploid = NULL,
                              width = NULL, height = NULL,
                              filename = "Figure6_Nash") {
  
  if (is.null(es_data)) {
    es_data <- load_sim(model_name, condition, sigma = sigma,
                        diploid_filter = diploid)
  }
  
  pA <- fig_nash_violation_map(model_name, es_data)
  
  pB <- if (!is.null(es_data) && "mS" %in% names(es_data))
    fig_slope_distribution(es_data)
  else
    ggplot() + theme_void() + labs(title = "(no ER data)")
  
  combined <- (pA | pB) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 13
    h <- if (!is.null(height)) height else 6
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}

# ============================================================================
# §10  FIGURE 7 -- Strategy parameter evolution & stability
# ============================================================================

fig_strategy_evolution <- function(es_data = NULL, model_name = "acute",
                                   condition = "ERhost_ERpath", sigma = 0.1,
                                   diploid = NULL,
                                   filename = "Figure_strategy_params",
                                   width = NULL, height = NULL) {
  
  if (is.null(es_data)) {
    es_data <- load_sim(model_name, condition, sigma = sigma,
                        diploid_filter = diploid)
    if (is.null(es_data) || nrow(es_data) == 0) {
      warning("No data found for ", model_name, "/", condition)
      return(invisible(NULL))
    }
  }
  
  thin <- es_data %>%
    mutate(row_num = row_number()) %>%
    filter(row_num %% 10 == 0)
  
  # Auto-detect appropriate x-axis
  gen_range <- range(thin$gen, na.rm = TRUE)
  use_log <- gen_range[2] > 1e5  # log scale for long runs
  
  if (use_log) {
    # Build nice log10 breaks within data range
    log_lo <- floor(log10(max(gen_range[1], 1)))
    log_hi <- ceiling(log10(gen_range[2]))
    log_breaks <- 10^seq(log_lo, log_hi)
    log_x <- scale_x_log10(
      breaks = log_breaks,
      labels = trans_format("log10", math_format(10^.x))
    )
  } else {
    tax <- auto_time_axis(thin)
    log_x <- scale_x_continuous(
      limits = tax$lims, breaks = tax$breaks, labels = tax$labels
    )
  }
  
  p_bS <- ggplot(thin, aes(gen, bS)) +
    geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
    log_x + labs(x = "Generation", y = "Host intercept") + mytheme
  
  p_mS <- ggplot(thin, aes(gen, mS)) +
    geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    log_x + labs(x = "Generation", y = "Host slope") + mytheme
  
  p_bV <- ggplot(thin, aes(gen, bV)) +
    geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
    log_x + labs(x = "Generation", y = "Pathogen intercept") + mytheme
  
  p_mV <- ggplot(thin, aes(gen, mV)) +
    geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    log_x + labs(x = "Generation", y = "Pathogen slope") + mytheme
  
  combined <- (p_bS | p_mS) / (p_bV | p_mV) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  
  w <- if (!is.null(width)) width else 7.5
  h <- if (!is.null(height)) height else 9
  ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
  combined
}


# ============================================================================
# §11  TIME SERIES STATISTICS
# ============================================================================

calc_cv <- function(x) {
  x <- na.omit(x)
  if (length(x) < 2) return(NA_real_)
  sd(x) / mean(x)
}

calc_spectral_slope <- function(x) {
  x <- na.omit(x)
  if (length(x) < 50) return(NA_real_)
  tryCatch({
    x_dt <- residuals(lm(x ~ seq_along(x)))
    spec <- spectrum(x_dt, plot = FALSE)
    freq <- spec$freq[-1]; power <- spec$spec[-1]
    valid <- freq > 0 & power > 0
    if (sum(valid) < 10) return(NA_real_)
    -coef(lm(log10(power[valid]) ~ log10(freq[valid])))[2]
  }, error = function(e) NA_real_)
}

calc_acf_at_lag <- function(x, lag = 1) {
  x <- na.omit(x)
  if (length(x) < lag + 10) return(NA_real_)
  tryCatch(acf(x, lag.max = lag, plot = FALSE)$acf[lag + 1],
           error = function(e) NA_real_)
}

calc_correlation_length <- function(x, threshold = 0.1, max_lag = NULL) {
  x <- na.omit(x)
  if (length(x) < 50) return(NA_real_)
  # Only skip truly zero-variance signals (numerical noise)
  if (sd(x) < .Machine$double.eps * 100) return(NA_real_)
  if (is.null(max_lag)) max_lag <- min(5000, floor(length(x) / 2))
  tryCatch({
    acf_vals <- as.numeric(acf(x, lag.max = max_lag, plot = FALSE)$acf[-1])
    below <- which(abs(acf_vals) < threshold)
    if (length(below) > 0) below[1] else max_lag
  }, error = function(e) NA_real_)
}


# ============================================================================
# §12  STEP SIZE ANALYSIS
# ============================================================================

calc_step_sizes <- function(df) {
  df %>%
    filter(event == "post") %>%
    arrange(gen) %>%
    mutate(
      delta_v = abs(v - lag(v)),
      delta_s = abs(s - lag(s))
    ) %>%
    filter(!is.na(delta_v))
}


# ============================================================================
# §13  NEUTRAL DRIFT ANALYSIS
# ============================================================================

identify_neutral_events <- function(df,
                                    geno_threshold = 0.01,
                                    pheno_threshold = 0.001) {
  df %>%
    arrange(gen) %>%
    mutate(
      delta_v   = abs(v - lag(v)),
      delta_s   = abs(s - lag(s)),
      pheno     = sqrt(delta_v^2 + delta_s^2),
      delta_bV  = abs(bV - lag(bV)),
      delta_mV  = abs(mV - lag(mV)),
      delta_bS  = abs(bS - lag(bS)),
      delta_mS  = abs(mS - lag(mS)),
      geno      = sqrt(delta_bV^2 + delta_mV^2 + delta_bS^2 + delta_mS^2),
      is_neutral = (geno > geno_threshold) & (pheno < pheno_threshold),
      decoupling = geno / (pheno + 1e-10)
    ) %>%
    filter(!is.na(pheno))
}


# ============================================================================
# §14  BOUNDARY ANALYSIS
# ============================================================================

calc_boundary_occupancy <- function(df, model_name = "acute",
                                    threshold = NULL) {
  dom <- TRAIT_DOMAIN[[model_name]]
  if (is.null(threshold)) {
    threshold <- (dom[2] - dom[1]) * 0.02  # 2% of domain
  }
  lo <- dom[1] + threshold
  hi <- dom[2] - threshold
  
  df %>%
    summarise(
      n = n(),
      host_at_boundary = mean(s < lo | s > hi),
      path_at_boundary = mean(v < lo | v > hi),
      any_at_boundary  = mean((s < lo | s > hi) | (v < lo | v > hi)),
      host_lower = mean(s < lo),
      host_upper = mean(s > hi),
      path_lower = mean(v < lo),
      path_upper = mean(v > hi),
      mean_v     = mean(v),
      mean_s     = mean(s),
      .groups = "drop"
    )
}

# ============================================================================
# §15  CROSS-CONDITION COMPARISON FIGURES
# ============================================================================
#
# Compare time-series statistics, step sizes, neutral drift, and dwell times
# across conditions, sigma values, and pinned-partner experiments.
#
# Shared color palette — high contrast, ordered ET/ET -> ER/ER
CONDITION_COLORS <- c(
  "ET / ET"            = "#1B9E77",   # teal
  "ET host / ER path"  = "#D95F02",   # orange
  "ER host / ET path"  = "#7570B3",   # purple
  "ER / ER"            = "#E7298A"    # magenta
)
scale_fill_condition <- function(...)
  scale_fill_manual(values = CONDITION_COLORS, ...)
scale_color_condition <- function(...)
  scale_color_manual(values = CONDITION_COLORS, ...)
#
# All functions accept:
#   sigma = 0.1       — one sigma (compare conditions)
#   sigma = NULL       — all sigmas (faceted by sigma)
#   include_pinned     — also load fixed-host / fixed-path runs
#

# --- Shared helper: load experiments from catalog ---
load_all_conditions <- function(model_name, sigma = 0.1, diploid = NULL,
                                include_pinned = FALSE, gamma_filter = 0.01,
                                tag_filter = NA) {
  cat <- discover_experiments()
  
  # Filter to model
  sub <- cat %>% filter(fitness == model_name)
  
  # Tag filter: NA (default) = exclude tagged runs; NULL = all; string = match
  if (is.na(tag_filter)) {
    sub <- sub %>% filter(is.na(tag))
  } else if (!is.null(tag_filter)) {
    sub <- sub %>% filter(!is.na(tag) & tag == tag_filter)
  }
  
  # Diploid filter
  if (!is.null(diploid)) sub <- sub %>% filter(diploid == !!diploid)
  
  # Sigma filter
  if (!is.null(sigma)) {
    sub <- sub %>% filter(abs(std_dev_move - sigma) < 1e-6)
  }
  
  # Gamma filter — default to 0.01 to exclude gamma-sweep runs
  # NA gamma means legacy config (default 0.01), so include those too
  if (!is.null(gamma_filter)) {
    sub <- sub %>% filter(is.na(gamma) | abs(gamma - gamma_filter) < 1e-6)
  }
  
  # Exclude or include pinned runs
  
  if (!include_pinned) {
    sub <- sub %>% filter(is.na(fix_host), is.na(fix_path))
  }
  
  if (nrow(sub) == 0) return(tibble())
  
  # Nice condition labels — ordered ET/ET -> mixed -> ER/ER
  cond_labels <- c(
    "EThost_ETpath" = "ET / ET",
    "EThost_ERpath" = "ET host / ER path",
    "ERhost_ETpath" = "ER host / ET path",
    "ERhost_ERpath" = "ER / ER"
  )
  cond_order <- c("ET / ET", "ET host / ER path",
                  "ER host / ET path", "ER / ER")
  
  all_df <- load_sim_set(sub) %>%
    mutate(
      scenario = factor(cond_labels[condition], levels = cond_order),
      sigma_label = sprintf("\u03c3 = %g", sigma)
    )
  
  # Label pinned runs
  if (include_pinned) {
    all_df <- all_df %>%
      mutate(
        run_type = case_when(
          !is.na(fix_host) & !is.na(fix_path) ~ "both fixed",
          !is.na(fix_host)  ~ "host fixed",
          !is.na(fix_path)  ~ "path fixed",
          TRUE              ~ "coevolving"
        )
      )
  }
  
  all_df
}


# =============================================================================
# REPLICATE-AWARE LOADING AND PLOTTING
# =============================================================================

#' Load all replicates for a given model and condition(s)
#'
#' Discovers runs tagged with rep1, rep2, ... and loads them into a single
#' data frame with a `rep` column. Runs without a rep tag are treated as rep=0
#' (the original / baseline run).
#'
#' @param model_name Fitness model name
#' @param sigma Step size filter (default 0.1)
#' @param diploid Diploid filter (NULL = any)
#' @param gamma_filter Gamma filter (default 0.01; NULL = all)
#' @param conditions Character vector of conditions to load (NULL = all 4)
#' @return tibble with columns: gen, v, s, condition, scenario, rep, ...
load_replicates <- function(model_name, sigma = 0.1, diploid = NULL,
                            gamma_filter = 0.01, conditions = NULL,
                            tag_filter = NA) {
  cat <- discover_experiments()
  
  sub <- cat %>% filter(fitness == model_name)
  
  # Tag filter: NA (default) = exclude tagged; NULL = all; string = match
  if (is.na(tag_filter)) {
    sub <- sub %>% filter(is.na(tag))
  } else if (!is.null(tag_filter)) {
    sub <- sub %>% filter(!is.na(tag) & tag == tag_filter)
  }
  
  if (!is.null(diploid)) sub <- sub %>% filter(diploid == !!diploid)
  if (!is.null(sigma))   sub <- sub %>% filter(abs(std_dev_move - sigma) < 1e-6)
  if (!is.null(gamma_filter)) {
    sub <- sub %>% filter(is.na(gamma) | abs(gamma - gamma_filter) < 1e-6)
  }
  # Exclude pinned runs
  sub <- sub %>% filter(is.na(fix_host), is.na(fix_path))
  
  if (!is.null(conditions)) {
    sub <- sub %>% filter(condition %in% conditions)
  }
  
  if (nrow(sub) == 0) {
    warning("No experiments found for ", model_name)
    return(tibble())
  }
  
  # Nice labels
  cond_labels <- c(
    "EThost_ETpath" = "ET / ET",
    "EThost_ERpath" = "ET host / ER path",
    "ERhost_ETpath" = "ER host / ET path",
    "ERhost_ERpath" = "ER / ER"
  )
  cond_order <- c("ET / ET", "ET host / ER path",
                  "ER host / ET path", "ER / ER")
  
  # Load each row, tagging with rep
  all_dfs <- lapply(seq_len(nrow(sub)), function(i) {
    row <- sub[i, ]
    df <- read_csv(row$csv, show_col_types = FALSE)
    rep_val <- if (is.na(row$rep)) 0L else as.integer(row$rep)
    df %>% mutate(
      condition = row$condition,
      scenario  = factor(cond_labels[row$condition], levels = cond_order),
      rep       = rep_val,
      rep_label = paste0("rep ", rep_val)
    )
  })
  
  bind_rows(all_dfs)
}


#' Figure: Overlay replicate time series
#'
#' Plots v and s time series with replicates overlaid as semi-transparent lines.
#' Faceted by condition (columns).
#'
#' @param model_name Fitness model name
#' @param sigma Step size (default 0.1)
#' @param diploid Diploid filter
#' @param conditions Which conditions to plot (NULL = all)
#' @param alpha Transparency for replicate lines (default 0.4)
#' @param log_time Use log10 x-axis (default TRUE)
fig_replicate_timeseries <- function(model_name = "taylor",
                                     sigma = 0.1, diploid = NULL,
                                     conditions = NULL,
                                     alpha = 0.4, log_time = TRUE,
                                     width = NULL, height = NULL,
                                     filename = NULL) {
  
  df <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                        conditions = conditions)
  if (nrow(df) == 0) {
    warning("No replicate data found"); return(invisible(NULL))
  }
  
  n_reps <- length(unique(df$rep))
  n_conds <- length(unique(df$scenario))
  cat("Plotting", n_reps, "replicates across", n_conds, "conditions\n")
  
  # Virulence panel
  p_v <- ggplot(df, aes(x = gen, y = v, color = factor(rep), group = rep)) +
    geom_line(alpha = alpha, linewidth = 0.3) +
    facet_wrap(~scenario, nrow = 1, scales = "free_y") +
    labs(y = "Virulence (v)", x = NULL, color = "Replicate") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 10))
  
  # Clearance panel
  p_s <- ggplot(df, aes(x = gen, y = s, color = factor(rep), group = rep)) +
    geom_line(alpha = alpha, linewidth = 0.3) +
    facet_wrap(~scenario, nrow = 1, scales = "free_y") +
    labs(y = "Clearance (c)", x = "Generation", color = "Replicate") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_blank())
  
  if (log_time) {
    p_v <- p_v + scale_x_log10(labels = scales::label_comma())
    p_s <- p_s + scale_x_log10(labels = scales::label_comma())
  }
  
  # Use distinguishable colors for replicates
  rep_cols <- c("0" = "#2D3748", "1" = "#1B9E77", "2" = "#D95F02",
                "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E",
                "6" = "#E6AB02", "7" = "#A6761D", "8" = "#666666",
                "9" = "#1F78B4")
  p_v <- p_v + scale_color_manual(values = rep_cols, na.value = "#999999")
  p_s <- p_s + scale_color_manual(values = rep_cols, na.value = "#999999")
  
  p <- p_v / p_s + plot_annotation(
    title = paste0(str_to_title(model_name), " — Replicate Overlay"),
    subtitle = paste0(n_reps, " replicates, \u03c3 = ", sigma)
  )
  
  # Save
  if (is.null(filename)) {
    filename <- paste0("Replicate_timeseries_", model_name)
  }
  w <- width  %||% max(8, n_conds * 3)
  h <- height %||% 6
  ggsave(paste0(filename, ".pdf"), p, width = w, height = h)
  ggsave(paste0(filename, ".png"), p, width = w, height = h, dpi = 200)
  cat("Saved:", filename, ".pdf/.png\n")
  p
}


#' Figure: Replicate trait density overlay
#'
#' Shows density distributions of v and s across replicates, 
#' faceted by condition. Good for checking whether replicates
#' converge to similar distributions.
fig_replicate_density <- function(model_name = "taylor",
                                  sigma = 0.1, diploid = NULL,
                                  conditions = NULL,
                                  width = NULL, height = NULL,
                                  filename = NULL) {
  
  df <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                        conditions = conditions)
  if (nrow(df) == 0) {
    warning("No replicate data found"); return(invisible(NULL))
  }
  
  n_reps <- length(unique(df$rep))
  n_conds <- length(unique(df$scenario))
  
  long_df <- df %>%
    pivot_longer(cols = c(v, s), names_to = "trait",
                 values_to = "value") %>%
    mutate(trait = ifelse(trait == "v", "Virulence (v)", "Clearance (c)"))
  
  p <- ggplot(long_df, aes(x = value, fill = factor(rep))) +
    geom_density(alpha = 0.3, linewidth = 0.3) +
    facet_grid(trait ~ scenario, scales = "free") +
    labs(x = "Trait value", y = "Density", fill = "Replicate",
         title = paste0(str_to_title(model_name), " — Replicate Trait Distributions"),
         subtitle = paste0(n_reps, " replicates")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  
  rep_cols <- c("0" = "#2D3748", "1" = "#1B9E77", "2" = "#D95F02",
                "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E",
                "6" = "#E6AB02", "7" = "#A6761D", "8" = "#666666",
                "9" = "#1F78B4")
  p <- p + scale_fill_manual(values = rep_cols, na.value = "#999999")
  
  if (is.null(filename)) {
    filename <- paste0("Replicate_density_", model_name)
  }
  w <- width  %||% max(8, n_conds * 2.5)
  h <- height %||% 6
  ggsave(paste0(filename, ".pdf"), p, width = w, height = h)
  ggsave(paste0(filename, ".png"), p, width = w, height = h, dpi = 200)
  cat("Saved:", filename, ".pdf/.png\n")
  p
}


#' Figure: Time-series statistics across conditions
#' Computes CV, spectral slope, and correlation length for v and s
#' in sliding windows, then plots distributions as violin + box plots.
#' sigma = NULL facets by sigma; include_pinned adds fixed-partner runs.
fig_ts_stats <- function(model_name = "acute",
                         sigma = 0.1, diploid = NULL,
                         include_pinned = FALSE,
                         window = 2000, step = 500,
                         width = NULL, height = NULL,
                         filename = "Figure_ts_stats") {
  
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  
  # Grouping columns depend on what varies
  grp_cols <- "scenario"
  if (is.null(sigma)) grp_cols <- c(grp_cols, "sigma_label")
  if (include_pinned) grp_cols <- c(grp_cols, "run_type")
  
  # --- CV & spectral slope in sliding windows ---
  win_stats <- all_df %>%
    group_by(across(all_of(grp_cols))) %>%
    group_map(function(df, key) {
      n <- nrow(df)
      if (n < window) return(NULL)
      starts <- seq(1, n - window + 1, by = step)
      lapply(starts, function(s) {
        chunk <- df[s:(s + window - 1), ]
        bind_cols(key, tibble(
          cv_v    = calc_cv(chunk$v),
          cv_s    = calc_cv(chunk$s),
          slope_v = calc_spectral_slope(chunk$v),
          slope_s = calc_spectral_slope(chunk$s)
        ))
      }) %>% bind_rows()
    }) %>% bind_rows()
  
  if (nrow(win_stats) == 0) {
    warning("Windows too large for data"); return(invisible(NULL))
  }
  
  # --- Correlation length on the FULL time series per group ---
  corr_stats <- all_df %>%
    group_by(across(all_of(grp_cols))) %>%
    summarise(
      corr_len_v = calc_correlation_length(v),
      corr_len_s = calc_correlation_length(s),
      .groups = "drop"
    )
  
  # Combine into long format
  win_long <- win_stats %>%
    pivot_longer(-all_of(grp_cols),
                 names_to = "metric", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      trait = ifelse(grepl("_v$", metric), "virulence (v)", "clearance (c)"),
      stat  = case_when(
        grepl("^cv",    metric) ~ "CV",
        grepl("^slope", metric) ~ "Spectral slope"
      ),
      stat = factor(stat, levels = c("CV", "Spectral slope"))
    )
  
  corr_long <- corr_stats %>%
    pivot_longer(-all_of(grp_cols),
                 names_to = "metric", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      trait = ifelse(grepl("_v$", metric), "virulence (v)", "clearance (c)")
    )
  
  # Build x-axis labels
  if (include_pinned) {
    win_long <- win_long %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
    corr_long <- corr_long %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    win_long  <- win_long  %>% mutate(x_label = scenario)
    corr_long <- corr_long %>% mutate(x_label = scenario)
  }
  
  # Top panels: CV & spectral slope (windowed, violin+box)
  p_top <- ggplot(win_long, aes(x = x_label, y = value, fill = scenario)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
    facet_grid(stat ~ trait, scales = "free_y") +
    scale_fill_condition() +
    labs(x = NULL, y = NULL) +
    mytheme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none",
          strip.text = element_text(size = 12))
  
  # Bottom panel: correlation length (full series, bar chart)
  p_bot <- ggplot(corr_long, aes(x = x_label, y = value, fill = scenario)) +
    geom_col(alpha = 0.7, width = 0.6) +
    facet_wrap(~ trait) +
    scale_fill_condition() +
    labs(x = NULL, y = "Correlation length\n(generations)") +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none",
          strip.text = element_blank())
  
  p <- (p_top / p_bot) +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(title = paste0(model_name, " — time-series statistics"))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 9
    h <- if (!is.null(height)) height else 9
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


#' Figure: ACF decay curves — shows when autocorrelation is lost
#' Plots full ACF(lag) curves for v and s per condition, so you can see
#' the memory timescale rather than just ACF(1).
#' @param max_lag  Maximum lag in generations (default 10000)
#' @param thin_acf Plot every Nth lag point (avoids over-plotting)
fig_acf_decay <- function(model_name = "acute",
                          sigma = 0.1, diploid = NULL,
                          include_pinned = FALSE,
                          max_lag = 10000, thin_acf = 10,
                          width = NULL, height = NULL,
                          filename = "Figure_acf_decay") {
  
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  
  # Grouping
  grp_cols <- "scenario"
  if (is.null(sigma)) grp_cols <- c(grp_cols, "sigma_label")
  if (include_pinned) grp_cols <- c(grp_cols, "run_type")
  
  # Compute full ACF per group
  acf_list <- all_df %>%
    group_by(across(all_of(grp_cols))) %>%
    group_map(function(df, key) {
      n <- nrow(df)
      ml <- min(max_lag, n - 10)
      if (ml < 10) return(NULL)
      
      acf_v <- tryCatch(
        as.numeric(acf(df$v, lag.max = ml, plot = FALSE)$acf[-1]),
        error = function(e) NULL)
      acf_s <- tryCatch(
        as.numeric(acf(df$s, lag.max = ml, plot = FALSE)$acf[-1]),
        error = function(e) NULL)
      
      if (is.null(acf_v) && is.null(acf_s)) return(NULL)
      
      lags <- seq_len(ml)
      # Thin for plotting
      keep <- seq(1, ml, by = thin_acf)
      
      rows <- list()
      if (!is.null(acf_v))
        rows[[1]] <- bind_cols(key, tibble(
          lag = lags[keep], acf = acf_v[keep], trait = "virulence (v)"))
      if (!is.null(acf_s))
        rows[[2]] <- bind_cols(key, tibble(
          lag = lags[keep], acf = acf_s[keep], trait = "clearance (c)"))
      bind_rows(rows)
    }) %>% bind_rows()
  
  if (nrow(acf_list) == 0) {
    warning("Could not compute ACF"); return(invisible(NULL))
  }
  
  # Build line group label
  if (include_pinned) {
    acf_list <- acf_list %>%
      mutate(group_label = paste0(scenario, " (", run_type, ")"))
  } else {
    acf_list <- acf_list %>%
      mutate(group_label = as.character(scenario))
  }
  
  p <- ggplot(acf_list, aes(x = lag, y = acf, color = scenario)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
    geom_line(aes(linetype = if (include_pinned) run_type else NULL,
                  group = group_label),
              alpha = 0.8, linewidth = 0.6) +
    facet_wrap(~ trait) +
    scale_color_condition() +
    labs(x = "Lag (generations)", y = "Autocorrelation",
         title = paste0(model_name, " — ACF decay"),
         color = "Condition") +
    mytheme +
    theme(legend.position = "bottom")
  
  # Facet by sigma if multiple
  if (is.null(sigma) && n_distinct(acf_list$sigma_label) > 1) {
    p <- p + facet_grid(sigma_label ~ trait)
  }
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 10
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


#' Figure: Mutational step-size distributions across conditions
fig_step_sizes <- function(model_name = "acute",
                           sigma = 0.1, diploid = NULL,
                           include_pinned = FALSE,
                           width = NULL, height = NULL,
                           filename = "Figure_step_sizes") {
  
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  
  grp_cols <- "scenario"
  if (is.null(sigma)) grp_cols <- c(grp_cols, "sigma_label")
  if (include_pinned) grp_cols <- c(grp_cols, "run_type")
  
  steps <- all_df %>%
    group_by(across(all_of(grp_cols))) %>%
    group_modify(~ calc_step_sizes(.x)) %>%
    ungroup() %>%
    pivot_longer(c(delta_v, delta_s),
                 names_to = "trait", values_to = "step") %>%
    mutate(trait = ifelse(trait == "delta_v",
                          "|delta v|", "|delta c|"))
  
  if (include_pinned) {
    steps <- steps %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    steps <- steps %>% mutate(x_label = scenario)
  }
  
  p <- ggplot(steps, aes(x = x_label, y = step, fill = scenario)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_y_log10() +
    scale_fill_condition() +
    labs(x = NULL, y = "Step size (log scale)",
         title = paste0(model_name, " — mutational step sizes")) +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  
  if (is.null(sigma) && n_distinct(steps$sigma_label) > 1) {
    p <- p + facet_grid(sigma_label ~ trait, scales = "free_y")
  }
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 8
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


#' Figure: Neutral drift analysis across conditions
fig_neutral_drift <- function(model_name = "acute",
                              sigma = 0.1, diploid = NULL,
                              include_pinned = FALSE,
                              width = NULL, height = NULL,
                              filename = "Figure_neutral_drift") {
  
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  
  er_cols <- c("bS", "mS", "bV", "mV")
  has_er <- all(er_cols %in% names(all_df))
  if (!has_er) {
    warning("Need ER data (bS, mS, bV, mV) for neutral drift analysis")
    return(invisible(NULL))
  }
  
  grp_cols <- "scenario"
  if (is.null(sigma)) grp_cols <- c(grp_cols, "sigma_label")
  if (include_pinned) grp_cols <- c(grp_cols, "run_type")
  
  drift <- all_df %>%
    filter(!is.na(bS) & !is.na(mS) & !is.na(bV) & !is.na(mV)) %>%
    group_by(across(all_of(grp_cols))) %>%
    group_modify(~ identify_neutral_events(.x)) %>%
    ungroup()
  
  if (nrow(drift) == 0) {
    warning("No drift data computed"); return(invisible(NULL))
  }
  
  # Panel A: Fraction of neutral events
  frac_df <- drift %>%
    group_by(across(all_of(grp_cols))) %>%
    summarise(neutral_frac = mean(is_neutral, na.rm = TRUE),
              n = n(), .groups = "drop")
  
  if (include_pinned) {
    frac_df <- frac_df %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    frac_df <- frac_df %>% mutate(x_label = scenario)
  }
  
  pA <- ggplot(frac_df, aes(x = x_label, y = neutral_frac, fill = scenario)) +
    geom_col(alpha = 0.7, width = 0.6) +
    scale_fill_condition() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Neutral fraction",
         title = "Fraction of neutral mutations") +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  
  # Panel B: Decoupling ratio distribution
  if (include_pinned) {
    drift <- drift %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    drift <- drift %>% mutate(x_label = scenario)
  }
  
  pB <- ggplot(drift %>% filter(decoupling < quantile(decoupling, 0.99,
                                                      na.rm = TRUE)),
               aes(x = x_label, y = decoupling, fill = scenario)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8) +
    scale_y_log10() +
    scale_fill_condition() +
    labs(x = NULL, y = "Geno / pheno ratio (log)",
         title = "Genotype-phenotype decoupling") +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  
  # Facet by sigma if multiple
  if (is.null(sigma) && n_distinct(frac_df$sigma_label) > 1) {
    pA <- pA + facet_wrap(~ sigma_label)
    pB <- pB + facet_wrap(~ sigma_label)
  }
  
  combined <- (pA | pB) +
    plot_annotation(
      title = paste0(model_name, " — neutral drift analysis"),
      tag_levels = "A"
    ) &
    theme(plot.tag = element_text(face = "bold", size = 16))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 10
    h <- if (!is.null(height)) height else 5.5
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}


#' Figure: Dwell-time distributions across conditions
#' Regions: "nash", "stable", "boundary"
fig_dwell_times <- function(model_name = "acute",
                            region = "nash",
                            sigma = 0.1, diploid = NULL,
                            include_pinned = FALSE,
                            nash_radius = NULL,
                            width = NULL, height = NULL,
                            filename = "Figure_dwell_times") {
  
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  
  nash_pt <- find_nash(model_name)
  dom <- TRAIT_DOMAIN[[model_name]]
  if (is.null(nash_radius)) nash_radius <- (dom[2] - dom[1]) * 0.1
  
  grp_cols <- "scenario"
  if (is.null(sigma)) grp_cols <- c(grp_cols, "sigma_label")
  if (include_pinned) grp_cols <- c(grp_cols, "run_type")
  
  all_df <- all_df %>%
    mutate(in_region = switch(region,
                              nash = {
                                sqrt((v - nash_pt$v)^2 + (s - nash_pt$s)^2) < nash_radius
                              },
                              stable = {
                                # Only meaningful for runs where both players have ER (nonzero slopes)
                                has_slopes <- "mS" %in% names(all_df) & "mV" %in% names(all_df)
                                if (has_slopes)
                                  ifelse(is.na(mS) | is.na(mV) | (mS == 0 & mV == 0),
                                         NA, abs(mS * mV) < 1)
                                else
                                  rep(NA, n())
                              },
                              boundary = {
                                thresh <- (dom[2] - dom[1]) * 0.02
                                (v < dom[1] + thresh) | (v > dom[2] - thresh) |
                                  (s < dom[1] + thresh) | (s > dom[2] - thresh)
                              },
                              stop("Unknown region: ", region)
    ))
  
  all_df <- all_df %>% filter(!is.na(in_region))
  if (nrow(all_df) == 0) {
    warning("No valid data for region '", region, "'")
    return(invisible(NULL))
  }
  
  # Compute run lengths per group
  dwell_df <- all_df %>%
    arrange(across(all_of(grp_cols)), gen) %>%
    group_by(across(all_of(grp_cols))) %>%
    mutate(run_id = cumsum(in_region != lag(in_region, default = !in_region[1]))) %>%
    group_by(across(all_of(c(grp_cols, "run_id", "in_region")))) %>%
    summarise(dwell = n(), .groups = "drop") %>%
    filter(in_region)
  
  if (nrow(dwell_df) == 0) {
    warning("No dwell events for region '", region, "'")
    return(invisible(NULL))
  }
  
  if (include_pinned) {
    dwell_df <- dwell_df %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    dwell_df <- dwell_df %>% mutate(x_label = scenario)
  }
  
  region_labels <- c(nash = "near Nash", stable = "stable (|mS mV| < 1)",
                     boundary = "at boundary")
  
  p <- ggplot(dwell_df, aes(x = x_label, y = dwell, fill = scenario)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8) +
    scale_y_log10() +
    scale_fill_condition() +
    labs(x = NULL, y = "Dwell time (generations, log scale)",
         title = paste0(model_name, " — dwell times: ",
                        region_labels[region])) +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  
  if (is.null(sigma) && n_distinct(dwell_df$sigma_label) > 1) {
    p <- p + facet_wrap(~ sigma_label)
  }
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 7
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


#' Figure: Boundary occupancy across conditions
#' Bar chart showing fraction of time each trait (host clearance, pathogen
#' virulence) spends at trait-space boundaries, broken down by upper/lower.
fig_boundary_occupancy <- function(model_name = "acute",
                                   sigma = 0.1, diploid = NULL,
                                   include_pinned = FALSE,
                                   width = NULL, height = NULL,
                                   filename = "Figure_boundary_occupancy") {
  
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  
  grp_cols <- "scenario"
  if (is.null(sigma)) grp_cols <- c(grp_cols, "sigma_label")
  if (include_pinned) grp_cols <- c(grp_cols, "run_type")
  
  occ <- all_df %>%
    group_by(across(all_of(grp_cols))) %>%
    group_modify(~ calc_boundary_occupancy(.x, model_name)) %>%
    ungroup()
  
  # Pivot to long for plotting
  occ_long <- occ %>%
    dplyr::select(all_of(grp_cols), host_lower, host_upper,
                  path_lower, path_upper) %>%
    pivot_longer(-all_of(grp_cols),
                 names_to = "boundary", values_to = "fraction") %>%
    mutate(
      player = ifelse(grepl("^host", boundary),
                      "clearance (c)", "virulence (v)"),
      side = ifelse(grepl("lower$", boundary), "lower", "upper")
    )
  
  if (include_pinned) {
    occ_long <- occ_long %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    occ_long <- occ_long %>% mutate(x_label = scenario)
  }
  
  p <- ggplot(occ_long, aes(x = x_label, y = fraction, fill = side)) +
    geom_col(position = "stack", alpha = 0.8, width = 0.6) +
    facet_wrap(~ player) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("lower" = "#4393C3", "upper" = "#D6604D"),
                      name = "Boundary") +
    labs(x = NULL, y = "Fraction of time at boundary",
         title = paste0(model_name, " — boundary occupancy")) +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  
  if (is.null(sigma) && "sigma_label" %in% names(occ_long) &&
      n_distinct(occ_long$sigma_label) > 1) {
    p <- p + facet_grid(sigma_label ~ player)
  }
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 8
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


#' Figure: Trait density distributions across conditions
#' Overlaid density curves for v and s, one per condition, with Nash
#' equilibrium marked as a vertical line.
fig_trait_density <- function(model_name = "acute",
                              sigma = 0.1, diploid = NULL,
                              include_pinned = FALSE,
                              width = NULL, height = NULL,
                              filename = "Figure_trait_density") {
  
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  
  nash_pt <- find_nash(model_name)
  disp <- if (model_name %in% names(TRAIT_DISPLAY))
    TRAIT_DISPLAY[[model_name]] else TRAIT_DOMAIN[[model_name]]
  
  # Clean axis breaks
  ax_breaks <- if (disp[2] <= 1) c(0, 0.25, 0.5, 0.75, 1) else pretty(disp, n = 5)
  
  scenarios <- levels(all_df$scenario)
  if (is.null(scenarios)) scenarios <- sort(unique(all_df$scenario))
  n_scen <- length(scenarios)
  
  # Thin scatter points for context (especially useful for concentrated conditions)
  set.seed(42)
  thin_n <- 2000
  all_thin <- all_df %>%
    group_by(scenario) %>%
    filter(row_number() %in% sample(seq_len(dplyr::n()),
                                    size = min(thin_n, dplyr::n()))) %>%
    ungroup()
  
  panels <- lapply(seq_along(scenarios), function(i) {
    sc <- scenarios[i]
    df_sc <- all_df %>% filter(scenario == sc)
    df_thin <- all_thin %>% filter(scenario == sc)
    
    show_y <- (i == 1)
    
    is_mid <- (i == ceiling(n_scen / 2))
    
    p <- ggplot(df_sc, aes(x = v, y = s)) +
      geom_point(data = df_thin, aes(x = v, y = s),
                 color = "grey60", size = 0.1, alpha = 0.3,
                 inherit.aes = FALSE) +
      stat_density_2d(aes(fill = after_stat(density)),
                      geom = "raster", contour = FALSE, alpha = 0.85) +
      scale_fill_viridis_c(option = "viridis") +
      geom_point(data = nash_pt, aes(x = v, y = s),
                 color = "red", size = 3, shape = 4, stroke = 1.5,
                 inherit.aes = FALSE) +
      coord_fixed(xlim = disp, ylim = disp) +
      scale_x_continuous(breaks = ax_breaks) +
      scale_y_continuous(breaks = ax_breaks) +
      labs(title = sc,
           x = if (is_mid) "v (virulence)" else NULL,
           y = if (show_y) "c (clearance)" else NULL) +
      mytheme +
      theme(plot.title = element_text(size = 11, hjust = 0.5),
            legend.position = "none")
    
    if (!show_y) {
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    }
    p
  })
  
  combined <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = paste0(model_name, " — trait distributions"),
      tag_levels = "A"
    ) &
    theme(plot.tag = element_text(face = "bold", size = 16))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 3.5 * n_scen
    h <- if (!is.null(height)) height else 4.5
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}

#-------
# ============================================================================
# §11  Gamma Sweep -- Mutation Rate Asymmetry Comparison
# ============================================================================

#' Load all gamma-sweep experiments for a model
#' Returns a data frame with gamma as an additional column
load_gamma_sweep <- function(model_name, sigma = 0.1, diploid = TRUE,
                             conditions = NULL, max_pts = 2000) {
  cat <- discover_experiments() %>%
    filter(fitness == model_name,
           is.na(fix_host), is.na(fix_path))
  
  if (!is.null(sigma))   cat <- cat %>% filter(abs(std_dev_move - sigma) < 1e-6)
  if (!is.null(diploid)) cat <- cat %>% filter(diploid == !!diploid)
  if (!is.null(conditions)) cat <- cat %>% filter(condition %in% conditions)
  
  # Need multiple gamma values
  if (n_distinct(cat$gamma) < 2) {
    warning("Only ", n_distinct(cat$gamma), " gamma value(s) found — ",
            "run a gamma sweep first.\n",
            "  Available gammas: ", paste(unique(cat$gamma), collapse = ", "))
    return(tibble())
  }
  
  cat("  Loading gamma sweep: ", model_name, 
      " (", n_distinct(cat$gamma), " gamma values × ",
      n_distinct(cat$condition), " conditions)\n")
  
  cond_labels <- c(
    "EThost_ETpath" = "ET / ET",
    "EThost_ERpath" = "ET host / ER path",
    "ERhost_ETpath" = "ER host / ET path",
    "ERhost_ERpath" = "ER / ER"
  )
  cond_order <- c("ET / ET", "ET host / ER path",
                  "ER host / ET path", "ER / ER")
  
  all_df <- load_sim_set(cat) %>%
    mutate(
      scenario = factor(cond_labels[condition], levels = cond_order),
      gamma_label = sprintf("\u03b3 = %g", gamma),
      gamma_desc = case_when(
        gamma < 0.1  ~ "path-fast",
        gamma > 0.9  ~ "host-fast",
        abs(gamma - 0.5) < 0.05 ~ "equal",
        gamma < 0.5  ~ "path-biased",
        TRUE         ~ "host-biased"
      )
    )
  
  # Order gamma labels by value
  gamma_order <- sort(unique(all_df$gamma))
  all_df$gamma_label <- factor(
    all_df$gamma_label,
    levels = sprintf("\u03b3 = %g", gamma_order)
  )
  
  all_df
}


#' Figure: Gamma sweep time series grid
#' Rows = gamma values, Columns = conditions
#' Shows v and s in separate grids
fig_gamma_timeseries <- function(model_name = "acute",
                                 sigma = 0.1, diploid = TRUE,
                                 conditions = NULL,
                                 max_pts = 100,
                                 width = NULL, height = NULL,
                                 filename = "Gamma_sweep_timeseries") {
  
  all_df <- load_gamma_sweep(model_name, sigma, diploid, conditions, max_pts)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  thin_df <- all_df %>%
    group_by(gamma_label, scenario) %>%
    group_modify(~ thin_for_plot(.x, max_pts = max_pts)) %>%
    ungroup()
  
  yax <- auto_trait_axis(model_name)
  
  p_v <- ggplot(thin_df, aes(gen, v, color = scenario)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_grid(gamma_label ~ scenario) +
    scale_y_continuous(limits = yax$lims, breaks = yax$breaks) +
    scale_color_condition() +
    labs(y = expression(italic(v) ~ "(virulence)"), x = NULL,
         title = paste0(model_name, " — virulence across \u03b3")) +
    mytheme +
    theme(legend.position = "none",
          strip.text = element_text(size = 9))
  
  p_s <- ggplot(thin_df, aes(gen, s, color = scenario)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_grid(gamma_label ~ scenario) +
    scale_y_continuous(limits = yax$lims, breaks = yax$breaks) +
    scale_color_condition() +
    labs(y = expression(italic(c) ~ "(clearance)"), x = "Evolutionary time",
         title = paste0(model_name, " — clearance across \u03b3")) +
    mytheme +
    theme(legend.position = "none",
          strip.text = element_text(size = 9))
  
  combined <- p_v / p_s
  
  if (!is.null(filename)) {
    n_gammas <- n_distinct(thin_df$gamma_label)
    n_conds  <- n_distinct(thin_df$scenario)
    w <- if (!is.null(width)) width else max(8, 2.5 * n_conds)
    h <- if (!is.null(height)) height else max(8, 1.5 * n_gammas * 2 + 2)
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}


#' Figure: Gamma sweep summary statistics
#' Compares CV, realized step sizes, and neutral drift fraction across gamma values
fig_gamma_summary <- function(model_name = "acute",
                              sigma = 0.1, diploid = TRUE,
                              conditions = NULL,
                              width = NULL, height = NULL,
                              filename = "Gamma_sweep_summary") {
  
  all_df <- load_gamma_sweep(model_name, sigma, diploid, conditions)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  # --- Panel 1: Realized step sizes by gamma and condition ---
  steps <- all_df %>%
    group_by(gamma_label, scenario) %>%
    group_modify(~ calc_step_sizes(.x)) %>%
    ungroup() %>%
    pivot_longer(c(delta_v, delta_s),
                 names_to = "trait", values_to = "step") %>%
    mutate(trait = ifelse(trait == "delta_v", "|Δv|", "|Δc|"))
  
  p_steps <- ggplot(steps, aes(x = gamma_label, y = step, fill = scenario)) +
    geom_violin(alpha = 0.4, scale = "width", position = position_dodge(0.8)) +
    geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8,
                 position = position_dodge(0.8)) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_y_log10() +
    scale_fill_condition() +
    labs(x = NULL, y = "Realized step size (log)",
         title = "Realized step sizes") +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  
  # --- Panel 2: CV of traits by gamma and condition ---
  cv_df <- all_df %>%
    group_by(gamma_label, scenario) %>%
    summarise(
      cv_v = sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE),
      cv_s = sd(s, na.rm = TRUE) / mean(s, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(c(cv_v, cv_s), names_to = "trait", values_to = "cv") %>%
    mutate(trait = ifelse(trait == "cv_v", "CV(v)", "CV(c)"))
  
  p_cv <- ggplot(cv_df, aes(x = gamma_label, y = cv, 
                            fill = scenario, group = scenario)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    facet_wrap(~ trait) +
    scale_fill_condition() +
    labs(x = NULL, y = "Coefficient of variation",
         title = "Trait variability") +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  
  # --- Panel 3: Neutral drift fraction (ER conditions only) ---
  er_cols <- c("bS", "mS", "bV", "mV")
  has_er <- all(er_cols %in% names(all_df))
  
  if (has_er) {
    drift <- all_df %>%
      filter(!is.na(bS) & !is.na(mS) & !is.na(bV) & !is.na(mV)) %>%
      group_by(gamma_label, scenario) %>%
      group_modify(~ {
        events <- identify_neutral_events(.x)
        if (nrow(events) == 0) return(tibble(neutral_frac = NA_real_))
        tibble(neutral_frac = mean(events$is_neutral, na.rm = TRUE))
      }) %>%
      ungroup() %>%
      filter(!is.na(neutral_frac))
    
    if (nrow(drift) > 0) {
      p_drift <- ggplot(drift, aes(x = gamma_label, y = neutral_frac, 
                                   fill = scenario)) +
        geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
        scale_fill_condition() +
        labs(x = NULL, y = "Fraction neutral",
             title = "Neutral drift events") +
        mytheme +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              legend.position = "none")
      
      combined <- (p_steps / (p_cv | p_drift)) +
        plot_annotation(
          title = paste0(model_name, " — mutation rate asymmetry (\u03b3) sweep"),
          tag_levels = "A"
        )
    } else {
      combined <- (p_steps / p_cv) +
        plot_annotation(
          title = paste0(model_name, " — mutation rate asymmetry (\u03b3) sweep"),
          tag_levels = "A"
        )
    }
  } else {
    combined <- (p_steps / p_cv) +
      plot_annotation(
        title = paste0(model_name, " — mutation rate asymmetry (\u03b3) sweep"),
        tag_levels = "A"
      )
  }
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 12
    h <- if (!is.null(height)) height else 10
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}


# ============================================================================
# GAMMA SWEEP FIGURES (mutation rate asymmetry)
# ============================================================================
# Run AFTER completing gamma sweep experiments:
#   python run_experiments.py -f acute --gamma-sweep 0.01,0.1,0.5,0.9,0.99 --diploid
#
# Then refresh catalog and generate figures:
#   refresh_catalog()
#   fig_gamma_timeseries("acute", filename = "Gamma_timeseries_acute")
#   fig_gamma_summary("acute", filename = "Gamma_summary_acute")
#
# For focused ER-ER only comparison:
#   fig_gamma_timeseries("acute", conditions = c("ERhost_ERpath"),
#                        filename = "Gamma_timeseries_acute_ERER")
# ============================================================================
# GAMMA SWEEP — Additional Analysis Functions
# ============================================================================
# Add these to Plots.R (after the existing fig_gamma_summary function)
#
# PURPOSE: These functions address the core scientific question of whether
# volatility in coevolutionary dynamics tracks ER status or mutation rate
# asymmetry (gamma). Together with the existing fig_gamma_timeseries() and
# fig_gamma_summary(), they provide a comprehensive test.
#
# RUNNING THE EXPERIMENTS — execute from the project root:
# -------------------------------------------------------
# Focused run (recommended first — ER-ER only, ~5 sims):
#   python run_experiments.py -f acute -c ERhost_ERpath \
#     --gamma-sweep 0.01,0.1,0.5,0.9,0.99 --diploid
#
# Full run (all conditions × 5 gammas = 20 sims):
#   python run_experiments.py -f acute \
#     --gamma-sweep 0.01,0.1,0.5,0.9,0.99 --diploid
#
# Minimal 3-point sweep (fastest, captures the essentials):
#   python run_experiments.py -f acute \
#     --gamma-sweep 0.01,0.5,0.99 --diploid
#
# Quick test (10K gens, just to check it works):
#   python run_experiments.py -f acute -c ERhost_ERpath \
#     --gamma-sweep 0.01,0.5,0.99 --diploid --quick
#
# Additional models:
#   python run_experiments.py -f minimal --gamma-sweep 0.01,0.5,0.99 --diploid
#   python run_experiments.py -f taylor  --gamma-sweep 0.01,0.5,0.99 --diploid
#
# After running, refresh the catalog in R:
#   refresh_catalog()
#   list_experiments()   # verify gamma runs appear
# -------------------------------------------------------


# ============================================================================
# 1. WHO-MUTATES FRACTION
# ============================================================================
# Shows what fraction of substitution events are host vs pathogen mutations
# at each gamma. Validates that gamma actually shifts the substitution balance
# as expected, and reveals whether ER status modifies the host/path ratio.

fig_gamma_who_mutates <- function(model_name = "acute",
                                  sigma = 0.1, diploid = TRUE,
                                  conditions = NULL,
                                  width = NULL, height = NULL,
                                  filename = "Gamma_who_mutates") {
  
  all_df <- load_gamma_sweep(model_name, sigma, diploid, conditions)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  # The mutator column records "host" or "path" for each substitution event
  if (!"mutator" %in% names(all_df)) {
    warning("No 'mutator' column found — need full CSV with mutator info")
    return(invisible(NULL))
  }
  
  # Compute host-mutation fraction per (gamma, condition)
  frac_df <- all_df %>%
    filter(mutator %in% c("host", "path")) %>%
    group_by(gamma_label, gamma, scenario) %>%
    summarise(
      n_host = sum(mutator == "host"),
      n_path = sum(mutator == "path"),
      n_total = n(),
      frac_host = n_host / n_total,
      .groups = "drop"
    )
  
  # Expected line: frac_host = gamma (if rates scale linearly)
  expected <- tibble(
    gamma = seq(0, 1, 0.01),
    expected_frac = gamma  # naive expectation
  )
  
  p <- ggplot(frac_df, aes(x = gamma, y = frac_host, 
                           color = scenario, shape = scenario)) +
    geom_line(data = expected, aes(x = gamma, y = expected_frac),
              inherit.aes = FALSE,
              color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 3, alpha = 0.9) +
    geom_line(aes(group = scenario), alpha = 0.5) +
    scale_color_condition() +
    annotate("text", x = 0.85, y = 0.15, label = "γ = frac(host)",
             color = "gray50", size = 3, fontface = "italic") +
    labs(x = expression(gamma ~ "(prob host mutates)"),
         y = "Fraction of substitutions that are host",
         title = paste0(model_name, " — who mutates?"),
         subtitle = "Dashed = naïve expectation (frac_host = γ)") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    mytheme +
    theme(legend.position = "bottom")
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 7
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


# ============================================================================
# 2. ER-PAIR SYMMETRY TEST
# ============================================================================
# The key test: is (ER-host / ET-path, γ=0.99) dynamically equivalent to
# (ET-host / ER-path, γ=0.01)?
#
# If yes → rate asymmetry × ER status interaction drives dynamics
# If no  → host-path biological asymmetry matters independently
#
# Compares trait distributions (KDE) and summary stats for "matched pairs"
# where the fast player is always the ER player vs always the ET player.

fig_gamma_symmetry_test <- function(model_name = "acute",
                                    sigma = 0.1, diploid = TRUE,
                                    width = NULL, height = NULL,
                                    filename = "Gamma_symmetry_test") {
  
  all_df <- load_gamma_sweep(model_name, sigma, diploid)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  # Define the matched pairs
  # Pair A: ER player is fast
  #   - ET-host / ER-path with γ=0.01 (path=ER is fast)
  #   - ER-host / ET-path with γ=0.99 (host=ER is fast)
  # Pair B: ER player is slow
  #   - ET-host / ER-path with γ=0.99 (path=ER is slow)
  #   - ER-host / ET-path with γ=0.01 (host=ER is slow)
  
  pairs <- all_df %>%
    filter(condition %in% c("EThost_ERpath", "ERhost_ETpath")) %>%
    mutate(
      pair_label = case_when(
        condition == "EThost_ERpath" & gamma < 0.1 ~ "ER-fast (path ER, γ=0.01)",
        condition == "ERhost_ETpath" & gamma > 0.9 ~ "ER-fast (host ER, γ=0.99)",
        condition == "EThost_ERpath" & gamma > 0.9 ~ "ER-slow (path ER, γ=0.99)",
        condition == "ERhost_ETpath" & gamma < 0.1 ~ "ER-slow (host ER, γ=0.01)",
        condition == "EThost_ERpath" & abs(gamma - 0.5) < 0.1 ~ "Equal (path ER, γ=0.5)",
        condition == "ERhost_ETpath" & abs(gamma - 0.5) < 0.1 ~ "Equal (host ER, γ=0.5)",
        TRUE ~ NA_character_
      ),
      speed_class = case_when(
        (condition == "EThost_ERpath" & gamma < 0.1) |
          (condition == "ERhost_ETpath" & gamma > 0.9) ~ "ER player fast",
        (condition == "EThost_ERpath" & gamma > 0.9) |
          (condition == "ERhost_ETpath" & gamma < 0.1) ~ "ER player slow",
        TRUE ~ "Equal rates"
      )
    ) %>%
    filter(!is.na(pair_label))
  
  if (nrow(pairs) == 0) {
    warning("Need gamma = 0.01 and 0.99 for both asymmetric conditions.\n",
            "  Run: python run_experiments.py -f ", model_name,
            " --gamma-sweep 0.01,0.5,0.99 --diploid")
    return(invisible(NULL))
  }
  
  yax <- auto_trait_axis(model_name)
  
  # Panel A: Virulence density by matched pair
  p_v <- ggplot(pairs, aes(x = v, fill = pair_label, color = pair_label)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    facet_wrap(~ speed_class, ncol = 1) +
    scale_x_continuous(limits = yax$lims) +
    labs(x = expression(italic(v) ~ "(virulence)"),
         y = "Density", fill = NULL, color = NULL) +
    mytheme +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8))
  
  # Panel B: Clearance density by matched pair
  p_s <- ggplot(pairs, aes(x = s, fill = pair_label, color = pair_label)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    facet_wrap(~ speed_class, ncol = 1) +
    scale_x_continuous(limits = yax$lims) +
    labs(x = expression(italic(c) ~ "(clearance)"),
         y = "Density", fill = NULL, color = NULL) +
    mytheme +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8))
  
  # Panel C: Summary stats comparison
  stats <- pairs %>%
    group_by(pair_label, speed_class) %>%
    summarise(
      cv_v = sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE),
      cv_s = sd(s, na.rm = TRUE) / mean(s, na.rm = TRUE),
      mean_v = mean(v, na.rm = TRUE),
      mean_s = mean(s, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(c(cv_v, cv_s), names_to = "metric", values_to = "value") %>%
    mutate(metric = ifelse(metric == "cv_v", "CV(v)", "CV(c)"))
  
  p_stats <- ggplot(stats, aes(x = pair_label, y = value, fill = speed_class)) +
    geom_col(alpha = 0.8, width = 0.7) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = NULL, y = "Value", fill = NULL) +
    mytheme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          legend.position = "none")
  
  combined <- (p_v | p_s) / p_stats +
    plot_annotation(
      title = paste0(model_name, " — ER speed symmetry test"),
      subtitle = "Do matched pairs (ER-fast vs ER-slow) show equivalent dynamics?",
      tag_levels = "A"
    ) +
    plot_layout(heights = c(2, 1))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 10
    h <- if (!is.null(height)) height else 10
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), combined, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  combined
}


# ============================================================================
# 3. TRAIT DENSITY × GAMMA HEATMAP
# ============================================================================
# For a single condition (e.g. ER-ER), shows how the 2D trait distribution
# shifts as gamma changes. Each panel = one gamma value, with hexbin or
# contour density in (v, c) space. Nash equilibrium marked.

fig_gamma_trait_density <- function(model_name = "acute",
                                    condition_filter = "ERhost_ERpath",
                                    sigma = 0.1, diploid = TRUE,
                                    width = NULL, height = NULL,
                                    filename = "Gamma_trait_density") {
  
  all_df <- load_gamma_sweep(model_name, sigma, diploid,
                             conditions = condition_filter)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  yax <- auto_trait_axis(model_name)
  
  # Get Nash equilibrium for reference
  nash <- tryCatch({
    mod <- FITNESS_MODELS[[model_name]]
    nash_eq(mod$fH, mod$fP, mod$params,
            TRAIT_DOMAIN[[model_name]][1], TRAIT_DOMAIN[[model_name]][2])
  }, error = function(e) list(v = NA, s = NA))
  
  p <- ggplot(all_df, aes(x = v, y = s)) +
    geom_hex(bins = 40, alpha = 0.9) +
    scale_fill_viridis_c(option = "magma", trans = "log10",
                         name = "Count") +
    facet_wrap(~ gamma_label, nrow = 1) +
    coord_fixed(xlim = yax$lims, ylim = yax$lims) +
    labs(x = expression(italic(v) ~ "(virulence)"),
         y = expression(italic(c) ~ "(clearance)"),
         title = paste0(model_name, " / ",
                        condition_filter, " — trait density across \u03b3")) +
    mytheme +
    theme(strip.text = element_text(size = 10))
  
  # Add Nash point if found
  if (!is.na(nash$v)) {
    p <- p + geom_point(data = data.frame(v = nash$v, s = nash$s),
                        aes(v, s), color = "white", shape = 4,
                        size = 3, stroke = 1.5)
  }
  
  if (!is.null(filename)) {
    n_gammas <- n_distinct(all_df$gamma_label)
    w <- if (!is.null(width)) width else max(8, 3 * n_gammas)
    h <- if (!is.null(height)) height else 4
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


# ============================================================================
# 4. SUBSTITUTION TEMPO COMPARISON
# ============================================================================
# Shows the rate of evolutionary change (substitutions per unit time)
# as a function of gamma. Separates host and pathogen substitution rates.
# This reveals whether gamma primarily controls WHO mutates or also
# changes the TOTAL rate of evolution.

fig_gamma_tempo <- function(model_name = "acute",
                            sigma = 0.1, diploid = TRUE,
                            conditions = NULL,
                            width = NULL, height = NULL,
                            filename = "Gamma_tempo") {
  
  all_df <- load_gamma_sweep(model_name, sigma, diploid, conditions)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  if (!"mutator" %in% names(all_df)) {
    warning("No 'mutator' column — need full CSV")
    return(invisible(NULL))
  }
  
  # Count substitution events per unit evolutionary time
  tempo_df <- all_df %>%
    filter(mutator %in% c("host", "path")) %>%
    group_by(gamma, gamma_label, scenario) %>%
    summarise(
      time_span = max(gen, na.rm = TRUE) - min(gen, na.rm = TRUE),
      n_host_subs = sum(mutator == "host"),
      n_path_subs = sum(mutator == "path"),
      n_total     = n(),
      .groups = "drop"
    ) %>%
    mutate(
      rate_host  = n_host_subs / time_span,
      rate_path  = n_path_subs / time_span,
      rate_total = n_total / time_span
    ) %>%
    pivot_longer(c(rate_host, rate_path, rate_total),
                 names_to = "rate_type", values_to = "rate") %>%
    mutate(rate_type = case_when(
      rate_type == "rate_host"  ~ "Host subs / gen",
      rate_type == "rate_path"  ~ "Pathogen subs / gen",
      rate_type == "rate_total" ~ "Total subs / gen"
    ))
  
  p <- ggplot(tempo_df, aes(x = gamma, y = rate, 
                            color = scenario, shape = rate_type)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_line(aes(group = interaction(scenario, rate_type)), alpha = 0.4) +
    facet_wrap(~ rate_type, scales = "free_y") +
    scale_color_condition() +
    labs(x = expression(gamma),
         y = "Substitution rate (events / generation)",
         title = paste0(model_name, " — evolutionary tempo across \u03b3"),
         color = "Condition") +
    mytheme +
    theme(legend.position = "bottom")
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 12
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


# ============================================================================
# 5. ACF COMPARISON ACROSS GAMMA
# ============================================================================
# Autocorrelation decay for v and s traits across gamma values.
# If ER drives long-range temporal correlations (punctuated equilibrium),
# ACF shape should be similar regardless of gamma.
# If gamma matters, fast-player traits should decorrelate faster.

fig_gamma_acf <- function(model_name = "acute",
                          condition_filter = "ERhost_ERpath",
                          sigma = 0.1, diploid = TRUE,
                          max_lag = 5000,
                          width = NULL, height = NULL,
                          filename = "Gamma_acf") {
  
  all_df <- load_gamma_sweep(model_name, sigma, diploid,
                             conditions = condition_filter)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  acf_list <- all_df %>%
    group_by(gamma_label) %>%
    group_modify(~ {
      acf_v <- acf(.x$v, lag.max = max_lag, plot = FALSE)
      acf_s <- acf(.x$s, lag.max = max_lag, plot = FALSE)
      bind_rows(
        tibble(lag = acf_v$lag[-1], acf = acf_v$acf[-1], trait = "v (virulence)"),
        tibble(lag = acf_s$lag[-1], acf = acf_s$acf[-1], trait = "c (clearance)")
      )
    }) %>%
    ungroup()
  
  p <- ggplot(acf_list, aes(x = lag, y = acf, color = gamma_label)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
    geom_line(alpha = 0.8, linewidth = 0.6) +
    facet_wrap(~ trait) +
    labs(x = "Lag (generations)", y = "Autocorrelation",
         title = paste0(model_name, " / ", condition_filter,
                        " — ACF across \u03b3"),
         color = NULL) +
    mytheme +
    theme(legend.position = "bottom")
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 10
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}


# ============================================================================
# CALL BLOCK 
# ============================================================================

# # --- Gamma sweep figures ---
# # (Requires gamma sweep data — run experiments first)
# #
# # Recommended execution sequence:
# #   1. Run experiments:
# #      python run_experiments.py -f acute --gamma-sweep 0.01,0.5,0.99 --diploid
# #   2. Refresh catalog:
# #      refresh_catalog()
# #   3. Generate figures:
#
# # Core figures (already in Plots.R):
# fig_gamma_timeseries("acute", filename = "Gamma_timeseries_acute")
# fig_gamma_summary("acute", filename = "Gamma_summary_acute")
#
# # New figures:
# fig_gamma_who_mutates("acute", filename = "Gamma_who_mutates_acute")
# fig_gamma_symmetry_test("acute", filename = "Gamma_symmetry_acute")
# fig_gamma_trait_density("acute", filename = "Gamma_trait_density_ERER_acute")
# fig_gamma_tempo("acute", filename = "Gamma_tempo_acute")
# fig_gamma_acf("acute", filename = "Gamma_acf_ERER_acute")
#
# # ER-ER only (focused):
# fig_gamma_timeseries("acute", conditions = "ERhost_ERpath",
#                       filename = "Gamma_ERER_timeseries_acute")
# fig_gamma_trait_density("acute", condition_filter = "ERhost_ERpath",
#                          filename = "Gamma_ERER_density_acute")
#
# # Asymmetric conditions only (for symmetry test):
# fig_gamma_timeseries("acute",
#                       conditions = c("EThost_ERpath", "ERhost_ETpath"),
#                       filename = "Gamma_asymmetric_timeseries_acute")
#
# # Other models:
# fig_gamma_timeseries("minimal", filename = "Gamma_timeseries_minimal")
# fig_gamma_summary("minimal", filename = "Gamma_summary_minimal")

#-------
#PLOTS#
#-------
# Overwrite default trait domains with model-specific ones if provided

# Trait domain per model.  Taylor uses rates (unbounded); others use [0,1].
TRAIT_DOMAIN <- list(
  acute   = c(0.001, 0.999),
  chronic = c(0.001, 0.999),
  minimal = c(0.001, 0.999),
  taylor  = c(0.01,  20.0)     # Nash ≈ (v*=9, c*=3)
)

# Clean axis limits for plotting (not the simulation clamp bounds)
TRAIT_DISPLAY <- list(
  acute   = c(0, 1),
  chronic = c(0, 1),
  minimal = c(0, 1),
  taylor  = c(0, 20)
)

# These do not need simulation data, so we can call them directly.
# landscapes for all simulations
fig_landscape(model_name = "acute", filename = "Fitness_landscapes_acute")
fig_landscape(model_name = "minimal", filename = "Fitness_landscapes_minimal")
fig_landscape(model_name = "taylor", filename = "Fitness_landscapes_taylor")

# strategy panels (analytical — no simulation data needed)
fig_strategy_panels(model_name = "acute",   filename = "Strategies_acute")
fig_strategy_panels(model_name = "minimal", filename = "Strategies_minimal")
fig_strategy_panels(model_name = "taylor",  filename = "Strategies_taylor")

#SIMULATION PLOTS#

# See everything at a glance
list_experiments()

# Base runs (auto-picks the untagged ones)
# Max_pts will thin the data
# there is step and smoothing functions too but they do not work well right now

fig_timeseries("acute", "Time_series_acute", diploid = TRUE, max_pts = 100,
               sigma = 0.1, width = 9, height = 10) 

fig_timeseries("minimal", "Time_series_minimal", diploid = TRUE, max_pts = 100, 
               sigma = 0.1, width = 9, height = 10)

fig_timeseries("taylor", "Time_series_taylor", diploid = TRUE, max_pts = 100,
               sigma = 0.1, width = 9, height = 10)

fig_timeseries("taylor", "Time_series_taylor_not_thinned", diploid = TRUE, max_pts = Inf,
               sigma = 0.1, width = 9, height = 10)

fig_timeseries("chronic", "Time_series_chronic", diploid = TRUE, max_pts = 1000,
               sigma = 0.1, width = 9, height = 10)

# Specific variants
# Works even if only some conditions exist for that variant
# (Only for acute)
# Step size experiments
fig_step_size_comparison(model_name = "acute", diploid = TRUE, 
                         filename = "Stepsize_sweeps_acute_EThost_ETpath", 
                         width = 8, height = 10)

# Host and Path fixed experiments
fig_pinned_comparison("acute", diploid = TRUE, sigma = 0.1, 
                      filename = "Pinned_partner_comparison_acute", 
                      width = 8, height = 10)


# Other figures

# Nash violation map — shows where in trait space the Nash condition is violated, 
# colored by magnitude of violation. Useful for understanding the geometry of the 
#landscape and why certain strategies are stable or not.

fig_nash_violation_map(model_name = "acute", diploid = T, resolution = 100,
                       width = 7, height = 5, 
                       filename = "Nash_violation_map_acute")

fig_nash_violation_map(model_name = "minimal", diploid = T, resolution = 100,
                       width = 7, height = 5, 
                       filename = "Nash_violation_map_minimal")

fig_nash_violation_map(model_name = "taylor", diploid = T, resolution = 100,
                       width = 7, height = 5, 
                       filename = "Nash_violation_map_taylor")


# Snapshots of evolutionary trajectories in trait space, colored by time, with Nash
fig_snapshots(model_name = "acute", condition = "ERhost_ERpath", n_panels = 6, 
              diploid = TRUE, sigma = 0.1, width = 7, height = 5, 
              filename = "ER-ER_Snapshots_acute")

fig_snapshots(model_name = "minimal", condition = "ERhost_ERpath", n_panels = 6, 
              diploid = TRUE, sigma = 0.1, width = 7, height = 5, 
              filename = "ER-ER_Snapshots_minimal")

fig_snapshots(model_name = "taylor", condition = "ERhost_ERpath", n_panels = 6, 
              diploid = TRUE, sigma = 0.1, width = 7, height = 5, 
              filename = "ER-ER_Snapshots_taylor")


# Combined hexbin of all trajectories in trait space, colored by time, with Nash
fig_hex_combined(models = "acute", diploid = TRUE, sigma = 0.1, 
                 width = 8, height = 4, 
                 filename = "ER-ER_strategy_points_acute") 

fig_hex_combined(models = "minimal", diploid = TRUE, sigma = 0.1, 
                 width = 8, height = 4, 
                 filename = "ER-ER_strategy_points_minimal") 

fig_hex_combined(models = "taylor", diploid = TRUE, sigma = 0.1, 
                 width = 8, height = 4, 
                 filename = "ER-ER_strategy_points_taylor") 

# Evolution of strategies over time, with Nash equilibrium marked, faceted by condition
fig_strategy_evolution(model_name = "acute", diploid = TRUE, sigma = 0.1, 
                       width = 8, height = 6, 
                       filename = "ER-ER_strategy_evolution_acute")

fig_strategy_evolution(model_name = "minimal", diploid = TRUE, sigma = 0.1, 
                       width = 8, height = 6, 
                       filename = "ER-ER_strategy_evolution_minimal")

fig_strategy_evolution(model_name = "taylor", diploid = TRUE, sigma = 0.1, 
                       width = 8, height = 6, 
                       filename = "ER-ER_strategy_evolution_taylor")

# Distribution of spectral slopes in the time series, which can indicate stability 
#and memory effects.
fig_slope_distribution(model_name = "acute", diploid = TRUE, sigma = 0.1, 
                       width = 6, height = 5, 
                       filename = "ER-ER_stability_acute")  

fig_slope_distribution(model_name = "minimal", diploid = TRUE, sigma = 0.1, 
                       width = 6, height = 5, 
                       filename = "ER-ER_stability_minimal")  

fig_slope_distribution(model_name = "taylor", diploid = TRUE, sigma = 0.1, 
                       width = 6, height = 5, 
                       filename = "ER-ER_stability_taylor")  


# Other TS stat, diagnostic figures: 
# CV of traits in sliding windows, spectral slope distribution, correlation length
fig_ts_stats("acute", diploid = TRUE, sigma = 0.1, filename = "TS_stats_acute")
fig_ts_stats("minimal", diploid = TRUE, sigma = 0.1, filename = "TS_stats_minimal")
fig_ts_stats("taylor", diploid = TRUE, sigma = 0.1, filename = "TS_stats_taylor")

# Step size distribution across conditions, which can indicate how the effective mutation

fig_step_sizes("acute", diploid = TRUE, sigma = 0.1, filename = "Realized_step_sizes_acute")
fig_step_sizes("minimal", diploid = TRUE, sigma = 0.1, filename = "Realized_step_sizes_minimal")
fig_step_sizes("taylor", diploid = TRUE, sigma = 0.1, filename = "Realized_step_sizes_taylor")

# Neutral drift analysis: fraction of mutations that are effectively neutral, and the
# distribution of their decoupling ratios (genotype vs phenotype change), which can
# indicate how much of the evolutionary dynamics is driven by drift vs selection, and
# how much genotypic change is decoupled from phenotypic change.
fig_neutral_drift("acute", diploid = TRUE, sigma = 0.1, filename = "Neutral_drift_acute")
fig_neutral_drift("minimal", diploid = TRUE, sigma = 0.1, filename = "Neutral_drift_minimal")
fig_neutral_drift("taylor", diploid = TRUE, sigma = 0.1, filename = "Neutral_drift_taylor")

# Dwell time distributions in different regions of trait space (near Nash, stable
# regions, boundaries), which can indicate how long populations tend to stay in these
# regions and how that differs across conditions.
fig_dwell_times("acute", region = "nash", diploid = TRUE, filename = "Dwell_nash_acute")
fig_dwell_times("acute", region = "stable", diploid = TRUE, filename = "Dwell_stable_acute")
fig_dwell_times("acute", region = "boundary", diploid = TRUE, filename = "Dwell_boundary_acute")

fig_dwell_times("minimal", region = "nash", diploid = TRUE, filename = "Dwell_nash_minimal")
fig_dwell_times("minimal", region = "stable", diploid = TRUE, filename = "Dwell_stable_minimal")
fig_dwell_times("minimal", region = "boundary", diploid = TRUE, filename = "Dwell_boundary_minimal")

fig_dwell_times("taylor", region = "nash", diploid = TRUE, filename = "Dwell_nash_taylor")
fig_dwell_times("taylor", region = "stable", diploid = TRUE, filename = "Dwell_stable_taylor")
fig_dwell_times("taylor", region = "boundary", diploid = TRUE, filename = "Dwell_boundary_taylor")

# Trait density distributions across conditions, with Nash equilibrium marked, which can
# indicate how the population is distributed in trait space and how close it is to the Nash
# equilibrium under different conditions.
fig_trait_density("acute", sigma = 0.1, diploid = TRUE, filename = "Trait_density_acute")
fig_trait_density("minimal", sigma = 0.1, diploid = TRUE, filename = "Trait_density_minimal")
fig_trait_density("taylor", sigma = 0.1, diploid = TRUE, filename = "Trait_density_taylor")

# Boundary occupancy: fraction of time spent at trait-space boundaries, which can indicate
# how much the population is pushed against the limits of trait space under different conditions.
fig_boundary_occupancy("acute", sigma = 0.1, diploid = TRUE, filename = "Boundary_occupancy_acute")
fig_boundary_occupancy("minimal", sigma = 0.1, diploid = TRUE, filename = "Boundary_occupancy_minimal")
fig_boundary_occupancy("taylor", sigma = 0.1, diploid = TRUE, filename = "Boundary_occupancy_taylor")