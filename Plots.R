# ============================================================================
# Plots.R -- plotting and analysis helpers for the Goldstein host-pathogen model
# Canan Karakoc
#
# Definitions only: sourcing this file draws nothing. The manuscript figures
# are assembled from these helpers in manuscript_figures.R. Run R from the
# repository root so the relative paths (results/, figures/) resolve.
#
# STRUCTURE:
#   §0  Setup (libraries, theme, saving)
#   §1  Fitness models (defined once, registered by name)
#   §2  Utilities
#   §3  Data loading
#   §4  Fig 1      model and response-rule geometry
#   §5  Fig 2, 4C  time series (also S2, S6, S10)
#   §6  Fig 3      mechanism of destabilisation (also S11)
#   §7  Fig 4      selection and tempo (also S7)
#   §8  S1, S3-S5  supplementary diagnostics
#   §9  S8-S9      tracking-strength model
#
# ADDING A NEW FITNESS MODEL:
#   1. Add its fitness functions in §1
#   2. Register it in FITNESS_MODELS (and TRAIT_DOMAIN / TRAIT_DISPLAY)
#   3. The loaders in §3 pick up its runs from config.json automatically
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
library(jsonlite)

dir.create("figures", showWarnings = FALSE)

# --- Global theme ---
mytheme <- theme_bw() +
  theme(
    axis.ticks.length   = unit(0.2, "cm"),
    legend.text         = element_text(size = 16),
    axis.text           = element_text(size = 18, color = "black"),
    axis.title          = element_text(size = 19),
    plot.title          = element_text(size = 18),
    panel.border        = element_rect(fill = NA, colour = "black", linewidth = 1),
    strip.text.x        = element_text(size = 18),
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


# Bold-tag theme reused across grouped figures
.tag_theme <- theme(plot.tag = element_text(face = "bold", size = 18),
                    plot.tag.position = c(-0.02, 1.04))

# iCloud Drive in ~/Documents intermittently triggers
# "Error in grDevices::dev.off() : write failed" when ggsave closes a large
# PDF: the sync daemon grabs the inode before R finishes flushing.  Workaround
# is to render to /tmp (outside iCloud), then copy in.
safe_ggsave <- function(filename, plot, ...) {
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(filename)))
  ggsave(tmp, plot, ...)
  file.copy(tmp, filename, overwrite = TRUE)
  file.remove(tmp)
  invisible(filename)
}

# ============================================================================
# §1  FITNESS MODELS
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


# ----- Tracking model (minimal + best-response "tracking" term, scaled by k) -----
# Same shape as the minimal model with an added term that rewards each player
# for tracking the opponent's trait.  Reduces EXACTLY to `minimal` at k = 0.
#   W_H = c(1-c)(1-v) + k c^2 (1-c) v      (c == s, host clearance)
#   W_P = v(1-v)(1-c) + k v^2 (1-v) c
# The tracking strength k lives in the params list so a single fitness function
# serves every k; register_tracking_k() (below) makes one registry entry per k.
PAR_TRACKING <- list(k = 1.0)

fH_tracking <- function(v, s, p = PAR_TRACKING) {
  s * (1 - s) * (1 - v) + p$k * s^2 * (1 - s) * v
}

fP_tracking <- function(v, s, p = PAR_TRACKING) {
  v * (1 - v) * (1 - s) + p$k * v^2 * (1 - v) * s
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
  ),
  tracking = list(               # base tracking model (k = 1); per-k variants
    fH     = fH_tracking,        # are added by register_tracking_k() below.
    fP     = fP_tracking,
    params = PAR_TRACKING,
    label  = "Tracking"
  )
)

# Trait domain per model.  Taylor uses rates (unbounded); others use [0,1].
TRAIT_DOMAIN <- list(
  acute    = c(0.001, 0.999),
  chronic  = c(0.001, 0.999),
  minimal  = c(0.001, 0.999),
  tracking = c(0.001, 0.999),
  taylor   = c(0.01,  30.0)     # Nash ≈ (v*=9, c*=3)
)

# Clean axis limits for plotting (not the simulation clamp bounds)
TRAIT_DISPLAY <- list(
  acute    = c(0, 1),
  chronic  = c(0, 1),
  minimal  = c(0, 1),
  tracking = c(0, 1),
  taylor   = c(0, 30)
)

# ----- Tracking k-variants: one registry entry per tracking strength ----------
# The whole plotting framework is keyed on a model-name string, so we register
# a distinct model per k (e.g. "tracking_k2").  Passing "tracking_k2" to any
# model-name-driven figure (e.g. fig_snapshots) uses the k=2 fitness;
# load_sim() maps the same name back to the on-disk fitness=="tracking" runs
# with TRACKING_K == 2.  So a k-sweep is just a loop over these names.
tracking_model_name <- function(k) sprintf("tracking_k%g", k)

# Extract k from a "tracking_k<k>" name; NA for the bare "tracking".
tracking_k_of <- function(model_name) {
  m <- regmatches(model_name, regexpr("(?<=_k)[0-9.]+$", model_name, perl = TRUE))
  if (length(m) == 0) NA_real_ else as.numeric(m)
}

register_tracking_k <- function(ks) {
  for (k in ks) {
    nm <- tracking_model_name(k)
    if (!is.null(FITNESS_MODELS[[nm]])) next          # already registered
    local({
      kk <- k
      FITNESS_MODELS[[nm]] <<- list(
        fH     = function(v, s, ...) s * (1 - s) * (1 - v) + kk * s^2 * (1 - s) * v,
        fP     = function(v, s, ...) v * (1 - v) * (1 - s) + kk * v^2 * (1 - v) * s,
        params = list(k = kk),
        label  = sprintf("Tracking (k=%g)", kk)
      )
    })
    TRAIT_DOMAIN[[nm]]  <<- c(0.001, 0.999)
    TRAIT_DISPLAY[[nm]] <<- c(0, 1)
  }
  invisible(NULL)
}

# k values of the tracking sweep (k = 0 is the minimal model)
TRACKING_KS <- c(0, 0.5, 1, 1.5, 2, 3, 4)
register_tracking_k(TRACKING_KS)


# ============================================================================
# §2  UTILITIES
# ============================================================================


# Model-aware clamping (defaults to [0,1])
clamp_trait <- function(x, model_name = NULL) {
  if (is.null(model_name)) {
    return(pmin(1, pmax(0, x)))
  }
  dom <- TRAIT_DOMAIN[[model_name]]
  pmin(dom[2], pmax(dom[1], x))
}

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
# §3  DATA LOADING (auto-discovery from config.json)
# ============================================================================
# Every run directory holds config.json + simulation.csv. The catalog is built
# from the JSON, so no paths are hard-coded. Point all loaders at one tree with
#   options(ggt.results_root = "results")   # the default


# Global: minimum generation to include when loading data. Burn-in is already
# excluded from simulation.csv, so the default keeps every recorded row.
MIN_GEN_CUTOFF <- 0

# Cached catalog — avoid re-scanning filesystem on every load
.catalog_cache <- new.env(parent = emptyenv())
.catalog_cache$data <- NULL
.catalog_cache$root <- NULL

#' Scan results/ tree and build a catalog of all experiments.
#' Results are cached; call discover_experiments(refresh = TRUE) to re-scan.
#' The root can also be set globally, so every figure reads the same tree:
#'   options(ggt.results_root = "results")
discover_experiments <- function(results_root = getOption("ggt.results_root", "results"),
                                 refresh = FALSE) {

  # Lazily create the cache if it wasn't defined (e.g. when sourcing only part
  # of this file interactively, so lines 438-440 never ran).
  if (!exists(".catalog_cache", envir = globalenv())) {
    .catalog_cache <<- new.env(parent = emptyenv())
    .catalog_cache$data <- NULL
    .catalog_cache$root <- NULL
  }

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
        tracking_k    = if (is.null(cfg$TRACKING_K)) NA_real_ else as.numeric(cfg$TRACKING_K),
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

#' Human-readable label for an experiment row
experiment_label <- function(row) {
  parts <- c(row$condition)
  if (identical(row$fitness, "tracking") &&
      !is.null(row$tracking_k) && !is.na(row$tracking_k))
    parts <- c(parts, sprintf("k=%g", row$tracking_k))
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
                     sigma = NULL, diploid_filter = NULL,
                     tag_filter = NA, tracking_k = NULL) {

  cat <- discover_experiments()
  condition <- SCENARIO_TO_CONDITION[scenario]
  if (is.na(condition)) condition <- scenario

  # Tracking k-variants ("tracking_k2") map to the on-disk fitness=="tracking"
  # runs, restricted to the matching TRACKING_K.  An explicit tracking_k arg
  # overrides the value encoded in the model name.
  fitness_model <- model
  if (grepl("^tracking", model)) {
    fitness_model <- "tracking"
    if (is.null(tracking_k)) {
      kn <- tracking_k_of(model)
      if (!is.na(kn)) tracking_k <- kn
    }
  }

  subset <- cat %>% filter(fitness == fitness_model, condition == !!condition)

  if (!is.null(tracking_k) && "tracking_k" %in% names(subset)) {
    subset <- subset %>%
      filter(!is.na(tracking_k) & abs(tracking_k - !!tracking_k) < 1e-6)
  }
  
  # Tag filter: NA (default) = exclude tagged; NULL = all; string = match
  # NOTE: check is.null() BEFORE is.na() — is.na(NULL) is logical(0), which
  # makes `if (is.na(tag_filter))` error with "argument is of length zero".
  if (is.null(tag_filter)) {
    # include all runs regardless of tag
  } else if (is.na(tag_filter)) {
    subset <- subset %>% filter(is.na(tag))
  } else {
    subset <- subset %>% filter(!is.na(tag) & tag == tag_filter)
  }
  
  if (nrow(subset) == 0) {
    warning(paste("No data for", model, "/", condition,
                  if (!is.null(tracking_k)) paste0(" (k=", tracking_k, ")") else "",
                  "\n  Available:", paste(unique(cat$condition[cat$fitness == fitness_model]),
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
    available <- cat %>% filter(fitness == fitness_model, condition == !!condition)
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

# --- Shared helper: load experiments from catalog ---
load_all_conditions <- function(model_name, sigma = 0.1, diploid = NULL,
                                include_pinned = FALSE, gamma_filter = 0.01,
                                tag_filter = NA, tag_prefix = NULL,
                                keep_pre = FALSE) {

  # Replicate-aware path. load_replicates handles both tagged runs (tag_prefix)
  # and untagged ones (rep taken from config), and always returns a 'rep'
  # column, so statistics computed downstream never mix lineages. Only pinned
  # runs, which load_replicates excludes, still take the old path below.
  if (!include_pinned) {
    df <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                          gamma_filter = gamma_filter, tag_filter = tag_filter,
                          tag_prefix = tag_prefix, keep_pre = keep_pre)
    if (nrow(df) == 0) return(tibble())
    if (!is.null(sigma))
      df <- df %>% mutate(sigma_label = sprintf("\u03c3 = %g", sigma))
    return(df)
  }

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
                            tag_filter = NA, tag_prefix = NULL,
                            reps = NULL,
                            min_gen = MIN_GEN_CUTOFF,
                            results_root = getOption("ggt.results_root", "results"),
                            keep_pre = FALSE) {
  cat <- discover_experiments(results_root)

  # Tracking models ("tracking_k2") live on disk as fitness == "tracking"
  # restricted to a TRACKING_K value — mirror load_sim()'s handling.
  fitness_model <- model_name
  tk <- NULL
  if (grepl("^tracking", model_name)) {
    fitness_model <- "tracking"
    tk <- tracking_k_of(model_name)
  }

  sub <- cat %>% filter(fitness == fitness_model)

  if (!is.null(tk) && "tracking_k" %in% names(sub)) {
    sub <- sub %>% filter(!is.na(tracking_k) & abs(tracking_k - !!tk) < 1e-6)
  }

  # Tag matching: tag_prefix takes precedence over tag_filter
  # tag_prefix matches tags starting with a prefix (e.g. "zoom" matches the
  # zoom runs' tag) and extracts rep number from suffix
  if (!is.null(tag_prefix)) {
    sub <- sub %>% filter(!is.na(tag) & grepl(paste0("^", tag_prefix), tag))
  } else if (is.na(tag_filter)) {
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

    # Skip empty or malformed CSVs
    if (nrow(df) == 0 || !"event" %in% names(df)) {
      warning("Skipping empty/malformed CSV: ", row$csv)
      return(NULL)
    }

    # Extract rep number: prefer config rep field, else parse from tag suffix
    if (!is.na(row$rep)) {
      rep_val <- as.integer(row$rep)
    } else if (!is.null(tag_prefix) && !is.na(row$tag)) {
      suffix <- sub(paste0("^", tag_prefix), "", row$tag)
      rep_val <- suppressWarnings(as.integer(suffix))
      if (is.na(rep_val)) rep_val <- i  # fallback: use row index
    } else {
      rep_val <- 0L
    }

    df %>%
      filter(keep_pre | event == "post", gen > min_gen) %>%
      mutate(
        condition = row$condition,
        scenario  = factor(cond_labels[row$condition], levels = cond_order),
        rep       = rep_val,
        rep_label = paste0("rep ", rep_val),
        host_pop  = row$host_pop,   # carried for omega / 2N normalisation
        path_pop  = row$path_pop,
        gamma     = row$gamma,      # prob_host_mutate, for the tempo ratio R
        omegaPath = suppressWarnings(as.numeric(omegaPath)),
        omegaHost = suppressWarnings(as.numeric(omegaHost))
      )
  })

  result <- bind_rows(all_dfs)

  # Optional: keep only a subset of replicates (e.g. reps = c(1, 2, 3))
  if (!is.null(reps)) {
    result <- result %>% filter(rep %in% reps)
    if (nrow(result) == 0)
      warning("No replicates matched reps = ", paste(reps, collapse = ", "))
  }

  n_reps <- length(unique(result$rep))
  n_conds <- length(unique(result$condition))
  cat("Loaded", n_reps, "replicates across", n_conds, "conditions for", model_name, "\n")
  result
}

#' Pair each recorded "pre" row with the "post" row of the same substitution.
#' The writer emits pre (state before step_generation) then post (state after)
#' for every recorded generation, so post - pre is exactly one substitution,
#' and the post row's dwell is the evolutionary time spent in the PRE state.
#' Needs data loaded with keep_pre = TRUE.
event_pairs <- function(df) {
  df %>%
    # one group per run; replicate numbers repeat across gamma in the tempo sweep
    group_by(across(any_of(c("scenario", "rep", "gamma")))) %>%
    mutate(nx_event = lead(event), nx_gen = lead(gen),
           v1 = lead(v), c1 = lead(s),
           bS1 = lead(bS), mS1 = lead(mS), bV1 = lead(bV), mV1 = lead(mV),
           mutator1 = lead(mutator),
           s_coef = suppressWarnings(as.numeric(lead(mutSelCoeff))),
           dwell1 = suppressWarnings(as.numeric(lead(dwell)))) %>%
    ungroup() %>%
    filter(event == "pre", nx_event == "post", nx_gen == gen) %>%
    transmute(scenario, condition, rep, gen, across(any_of("gamma")),
              mutator = mutator1, s_coef, dwell = dwell1,
              v0 = v, c0 = s, v1, c1,
              dbS = bS1 - bS, dmS = mS1 - mS, dbV = bV1 - bV, dmV = mV1 - mV,
              host_pop, path_pop)
}

#' Keep only finished runs: every completed run records the same number of rows,
#' so a run still being written (or killed) has fewer. Reports what it skips.
drop_unfinished_runs <- function(d) {
  keys <- intersect(c("condition", "gamma", "rep"), names(d))
  n <- d %>% filter(event == "post") %>% count(across(all_of(keys)), name = "n_post")
  full <- max(n$n_post)
  partial <- n %>% filter(n_post < full)
  if (nrow(partial))
    cat(sprintf("  skipping %d unfinished run(s): %s\n", nrow(partial),
                paste(do.call(paste, c(partial[keys], sep = " ")), collapse = "; ")))
  d %>% semi_join(filter(n, n_post == full), by = keys)
}

# ============================================================================
# §4  FIG 1 -- model and response-rule geometry (analytical)
# ============================================================================


#' Three landscape panels: host, pathogen, joint with Nash.
#' Returns list(host = ggplot, path = ggplot, joint = ggplot).
panels_landscape <- function(model_name = "minimal") {
  grid    <- make_fitness_grid(model_name)
  br      <- calc_best_responses(model_name)
  nash_pt <- find_nash(model_name)
  dom     <- if (model_name %in% names(TRAIT_DISPLAY))
    TRAIT_DISPLAY[[model_name]] else TRAIT_DOMAIN[[model_name]]
  
  grid$fH_norm    <- grid$fH / max(grid$fH)
  grid$fP_norm    <- grid$fP / max(grid$fP)
  grid$joint_norm <- grid$fH_norm * grid$fP_norm
  
  ax_breaks <- if (dom[2] <= 1) c(0, 0.5, 1) else pretty(dom, n = 4)
  
  pHost <- ggplot(grid, aes(x = v, y = s, z = fH_norm)) +
    geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
    geom_line(data = br$host, aes(x = v, y = s),
              color = "steelblue", linewidth = 1.5,
              inherit.aes = FALSE, linetype = "dashed") +
    labs(x = "v (virulence)", y = "c (clearance)") +
    coord_fixed(xlim = dom, ylim = dom, expand = FALSE) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    scale_fill_viridis_d(option = "viridis", guide = "none") +
    mytheme

  pPath <- ggplot(grid, aes(x = v, y = s, z = fP_norm)) +
    geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
    geom_line(data = br$path, aes(x = v, y = s),
              color = "lightcoral", linewidth = 1.5,
              inherit.aes = FALSE, linetype = "dashed") +
    labs(x = "v (virulence)", y = NULL) +
    coord_fixed(xlim = dom, ylim = dom, expand = FALSE) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    scale_fill_viridis_d(option = "viridis", guide = "none") +
    mytheme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  pJoint <- ggplot(grid, aes(x = v, y = s, z = joint_norm)) +
    geom_contour(aes(z = fH_norm), bins = 10, color = "steelblue", alpha = 0.5) +
    geom_contour(aes(z = fP_norm), bins = 10, color = "lightcoral", alpha = 0.5) +
    geom_line(data = br$host, aes(x = v, y = s),
              color = "steelblue", inherit.aes = FALSE,
              linewidth = 1.5, linetype = "dashed") +
    geom_line(data = br$path, aes(x = v, y = s),
              color = "lightcoral", inherit.aes = FALSE,
              linewidth = 1.5, linetype = "dashed") +
    geom_point(data = nash_pt, aes(x = v, y = s),
               color = "grey20", size = 4, inherit.aes = FALSE) +
    labs(x = "v (virulence)", y = NULL) +
    coord_fixed(xlim = dom, ylim = dom, expand = FALSE) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    mytheme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  list(host = pHost, path = pPath, joint = pJoint)
}

#' Four strategy-geometry panels: ET, host-slope, path-slope, both-slope.
#' Returns list(et, host_slope, path_slope, both_slope).
panels_strategy <- function(model_name = "minimal") {
  grid    <- make_fitness_grid(model_name, resolution = 300)
  br      <- calc_best_responses(model_name)
  nash_pt <- find_nash(model_name)
  dom     <- TRAIT_DOMAIN[[model_name]]
  v_star  <- nash_pt$v
  s_star  <- nash_pt$s
  
  col_host <- "darkblue"
  col_path <- "firebrick"
  
  ax_breaks <- if (dom[2] <= 1) c(0, 0.5, 1) else pretty(dom, n = 4)
  
  pad   <- (dom[2] - dom[1]) * 0.1
  v_seq <- seq(dom[1] - pad, dom[2] + pad, length.out = 500)
  s_seq <- seq(dom[1] - pad, dom[2] + pad, length.out = 500)
  
  host_line_fn <- function(v, bS, mS) clamp_trait(bS + mS * v, model_name)
  path_line_fn <- function(s, bV, mV) clamp_trait(bV + mV * s, model_name)
  
  make_panel <- function(bV, mV, bS, mS,
                         show_stable = TRUE,
                         boundary_pts = NULL,
                         show_yaxis = TRUE) {
    
    host_data <- data.frame(v = v_seq, s = host_line_fn(v_seq, bS, mS)) %>%
      filter(v >= dom[1], v <= dom[2], s >= dom[1], s <= dom[2])
    path_data <- data.frame(s = s_seq, v = path_line_fn(s_seq, bV, mV)) %>%
      filter(v >= dom[1], v <= dom[2], s >= dom[1], s <= dom[2])
    
    den   <- 1 - mV * mS
    v_int <- if (abs(den) > 1e-9) (bV + mV * bS) / den else v_star
    s_int <- bS + mS * v_int
    
    p <- ggplot(grid, aes(v, s)) +
      geom_contour(aes(z = fH), color = "steelblue", bins = 10,
                   linewidth = 0.3, alpha = 0.7) +
      geom_contour(aes(z = fP), color = "lightcoral", bins = 10,
                   linewidth = 0.3, alpha = 0.7) +
      geom_line(data = host_data, aes(v, s),
                linetype = "solid", linewidth = 1.5, color = col_host) +
      geom_line(data = path_data, aes(v, s),
                linetype = "solid", linewidth = 1.5, color = col_path) +
      coord_fixed(xlim = dom, ylim = dom, expand = FALSE) +
      scale_x_continuous(breaks = ax_breaks) +
      scale_y_continuous(breaks = ax_breaks) +
      mytheme

    if (show_yaxis) {
      p <- p + labs(x = "v (virulence)", y = "c (clearance)")
    } else {
      p <- p + labs(x = "v (virulence)", y = NULL) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    if (show_stable) {
      p <- p + geom_point(aes(x = v_int, y = s_int), size = 4, colour = "black")
    } else {
      p <- p + geom_point(aes(x = v_int, y = s_int), size = 4, shape = 21,
                          fill = "gray70", colour = "black", stroke = 1)
    }
    if (!is.null(boundary_pts)) {
      p <- p + geom_point(data = boundary_pts, aes(x = v, y = s),
                          size = 4, colour = "black")
    }
    p
  }
  
  # Panel D in original: both slopes large -> destabilization
  mS_D <- 2.5; mV_D <- 2.5
  bS_D <- s_star - mS_D * v_star
  bV_D <- v_star - mV_D * s_star
  if (bS_D > dom[1] - 1e-3) bS_D <- dom[1] - 1e-3
  if (bV_D > dom[1] - 1e-3) bV_D <- dom[1] - 1e-3
  boundary_D <- data.frame(
    v = c(clamp_trait(bV_D, model_name),
          clamp_trait(bV_D + mV_D * dom[2], model_name)),
    s = c(clamp_trait(bS_D, model_name),
          clamp_trait(bS_D + mS_D * dom[2], model_name))
  )
  
  list(
    et         = make_panel(v_star, 0, s_star, 0, TRUE, NULL, TRUE),
    host_slope = make_panel(v_star, 0,
                            s_star - 0.5 * v_star, 0.5,
                            TRUE, NULL, FALSE),
    path_slope = make_panel(v_star - 0.8 * s_star, 0.8,
                            s_star, 0, TRUE, NULL, FALSE),
    both_slope = make_panel(bV_D, mV_D, bS_D, mS_D,
                            FALSE, boundary_D, FALSE)
  )
}

fig_grouped_1_setup <- function(model_name = "minimal",
                                filename = "Fig1_setup",
                                width = 13, height = 9) {

  L <- panels_landscape(model_name)
  S <- panels_strategy (model_name)

  top    <- L$host       | L$path       | L$joint
  bottom <- S$et         | S$host_slope | S$path_slope | S$both_slope

  out <- (top / bottom) +
    plot_layout(heights = c(1.33, 1)) +
    plot_annotation(tag_levels = "A") &
    .tag_theme &
    theme(plot.margin = margin(12, 14, 4, 18))
  
  if (!is.null(filename)) {
    safe_ggsave(paste0("figures/", filename, ".pdf"), out,
                width = width, height = height)
    safe_ggsave(paste0("figures/", filename, ".png"), out,
                width = width, height = height, dpi = 300)
    cat("Saved:", filename, "\n")
  }
  out
}

# ============================================================================
# §5  TIME SERIES -- Fig 2, Fig 4C, S2, S6, S10
# ============================================================================


#' Compute sensible time-axis settings from data
#' Returns list(lims, breaks, labels) that can be passed to line_panel/omega_panel
auto_time_axis <- function(df, n_breaks = 3) {
  gen_range <- range(df$gen, na.rm = TRUE)
  lo <- gen_range[1]
  hi <- gen_range[2]
  
  # Nice labels
  fmt_label <- function(x) {
    if (x >= 1e6) paste0(format(x / 1e6, trim = TRUE), "M")
    else if (x >= 1e3) paste0(format(x / 1e3, trim = TRUE), "K")
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
    # Allow ticks 1% past either end (runs stop at e.g. 99 900, not 100 000),
    # and coarsen until pretty() gives <= n_breaks + 1 evenly spaced ticks
    tol <- 0.01 * (hi - lo)
    for (n in n_breaks:1) {
      brk <- pretty(c(lo, hi), n = n)
      brk <- brk[brk >= lo - tol & brk <= hi + tol]
      if (length(brk) <= n_breaks + 1) break
    }
    if (length(brk) == 0) brk <- c(lo, hi)
  }
  
  labs <- sapply(brk, fmt_label)
  # Pad the lower limit by 2% of the span so the first label isn't clipped,
  # and extend the upper limit to the last tick so it isn't dropped
  lo <- min(lo, brk); hi <- max(hi, brk)
  list(lims = c(lo - 0.02 * (hi - lo), hi), breaks = brk, labels = labs)
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

# Replicate color palette (colorblind-friendly, up to 9 replicates)
REP_COLORS <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3",
                "4" = "#E7298A", "5" = "#2D3748", "6" = "#E6AB02",
                "7" = "#A6761D", "8" = "#666666", "9" = "#1F78B4",
                "0" = "#66A61E")

# --- Line panel (v, c, W) ---
# model_name: if provided, uses TRAIT_DOMAIN for y-limits on trait variables
# x_lims/x_breaks/x_labels: if NULL, auto-detected from data
# use_step: if TRUE, uses geom_step instead of geom_line (better for SSWM data)
# highlight_rep: with has_reps, draw this rep in black over the others in grey
line_panel <- function(df, y_var, ylab = NULL,
                       show_xlab = FALSE, show_ylab = TRUE,
                       model_name = NULL,
                       x_lims = NULL, x_breaks = NULL, x_labels = NULL,
                       use_step = FALSE, has_reps = FALSE,
                       show_legend = FALSE, highlight_rep = NULL) {

  rep_guide <- if (show_legend) guide_legend(title = "rep", override.aes = list(alpha = 1, linewidth = 1)) else "none"

  # Auto-detect time axis from data if not specified
  if (is.null(x_lims)) {
    tax <- auto_time_axis(df)
    x_lims <- tax$lims; x_breaks <- tax$breaks; x_labels <- tax$labels
  }

  # Auto-detect trait axis: use model domain for v/s, auto-range for fitness
  is_trait <- y_var %in% c("v", "s")
  if (is_trait) {
    yax <- auto_trait_axis(model_name, df, y_var)
    y_lims <- yax$lims; y_breaks <- yax$breaks
  } else {
    # Auto-scale fitness panels from data (chronic W_H = 1/m can exceed 1)
    rng <- range(df[[y_var]], na.rm = TRUE)
    if (rng[2] <= 1.05) {
      y_lims <- c(0, 1); y_breaks <- c(0, 0.5, 1)
    } else {
      span <- rng[2] - rng[1]
      # Enforce minimum span (10% of midpoint) so near-constant series
      # don't zoom into numerical noise
      min_span <- max(0.1 * mean(rng), 0.1)
      if (span < min_span) {
        mid <- mean(rng)
        rng <- c(mid - min_span / 2, mid + min_span / 2)
      }
      pad <- (rng[2] - rng[1]) * 0.05
      y_lims <- c(max(0, rng[1] - pad), rng[2] + pad)
      y_breaks <- pretty(y_lims, n = 4)
    }
  }

  geom_fn <- if (use_step) geom_step else geom_line

  if (has_reps && "rep" %in% names(df) && !is.null(highlight_rep)) {
    # Other replicates as a grey backdrop, the highlighted one in black on top
    p <- ggplot(mapping = aes(x = gen, y = .data[[y_var]], group = rep)) +
      geom_fn(data = filter(df, rep != highlight_rep),
              colour = "grey80", linewidth = 0.25) +
      geom_fn(data = filter(df, rep == highlight_rep),
              colour = "black", linewidth = 0.35)
  } else if (has_reps && "rep" %in% names(df)) {
    n_reps <- length(unique(df$rep))
    lw <- if (n_reps <= 3) 0.4 else 0.3
    al <- if (n_reps <= 3) 0.7 else 0.5
    p <- ggplot(df, aes(x = gen, y = .data[[y_var]],
                        color = factor(rep, levels = names(REP_COLORS)),
                        group = rep)) +
      geom_fn(linewidth = lw, alpha = al) +
      scale_color_manual(values = REP_COLORS, guide = rep_guide, drop = TRUE)
  } else {
    p <- ggplot(df, aes(x = gen, y = .data[[y_var]])) +
      geom_fn(linewidth = 0.5, alpha = 0.85)
  }

  p <- p +
    scale_x_continuous(limits = x_lims, breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = y_lims, breaks = y_breaks) +
    coord_cartesian(xlim = x_lims) +
    mytheme

  if (show_ylab && !is.null(ylab)) p <- p + labs(y = ylab)
  else p <- p + labs(y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  if (show_xlab) p <- p + labs(x = NULL)
  else p <- p + labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

# --- Omega spike panel ---
# Omega = cumulative substitution rate per player per generation.
# Neutral rate = 1 (haploid) or 0.5 (diploid).
# omega > 1 means positive selection is accelerating substitutions;
# omega < 1 means most mutations are deleterious or nearly neutral.
omega_panel <- function(df, who = c("Path", "Host"),
                        show_xlab = FALSE, show_ylab = TRUE, ylab = NULL,
                        x_lims = NULL, x_breaks = NULL, x_labels = NULL,
                        has_reps = FALSE, show_legend = FALSE,
                        highlight_rep = NULL, normalize = FALSE) {
  who <- match.arg(who)
  # normalize: plot omega / 2N (each player's ceiling), so host and pathogen,
  # whose N differ 100-fold, share one 0-1 scale
  pop_col <- if (who == "Path") "path_pop" else "host_pop"
  normalize <- normalize && pop_col %in% names(df)
  rep_guide <- if (show_legend) guide_legend(title = "rep", override.aes = list(alpha = 1, linewidth = 1)) else "none"
  omega_col <- if (who == "Path") "omegaPath" else "omegaHost"

  # Auto-detect time axis from data if not specified
  if (is.null(x_lims)) {
    tax <- auto_time_axis(df)
    x_lims <- tax$lims; x_breaks <- tax$breaks; x_labels <- tax$labels
  }

  rng <- range(x_lims)

  # Prepare data: clamp zero/NA omega to tiny value so lines stay connected
  prep_omega <- function(d) {
    d %>%
      filter(gen >= rng[1], gen <= rng[2]) %>%
      mutate(y = suppressWarnings(as.numeric(.data[[omega_col]])) /
                 (if (normalize) 2 * .data[[pop_col]] else 1)) %>%
      filter(!is.na(y)) %>%
      mutate(y = pmax(y, 1e-10))
  }

  # Fixed y-axis limits: 10^-2 to 10^6
  all_vals <- prep_omega(df)$y
  if (length(all_vals) == 0) {
    # No valid data — return empty panel
    p <- ggplot() + theme_void()
    if (show_ylab && !is.null(ylab)) p <- p + labs(y = ylab)
    return(p)
  }
  y_lims   <- if (normalize) c(1e-8, 1.5) else c(1e-2, 1e8)
  y_breaks <- if (normalize) c(1e-8, 1e-4, 1) else c(1e-2, 1e2, 1e6)

  # Plot omega directly (no binning — keeps lines connected)
  if (has_reps && "rep" %in% names(df) && !is.null(highlight_rep)) {
    # Other replicates as a grey backdrop, the highlighted one in black on top
    dat <- prep_omega(df) %>% mutate(x = gen)
    p <- ggplot(mapping = aes(x = x, y = y, group = rep)) +
      geom_line(data = filter(dat, rep != highlight_rep),
                colour = "grey80", linewidth = 0.25) +
      geom_line(data = filter(dat, rep == highlight_rep),
                colour = "black", linewidth = 0.3)
  } else if (has_reps && "rep" %in% names(df)) {
    dat <- prep_omega(df) %>%
      mutate(x = gen, rep = factor(rep, levels = names(REP_COLORS)))

    n_reps <- length(unique(dat$rep))
    al <- if (n_reps <= 3) 0.6 else 0.4
    p <- ggplot(dat) +
      geom_line(aes(x = x, y = y, color = rep, group = rep),
                linewidth = 0.3, alpha = al) +
      scale_color_manual(values = REP_COLORS, guide = rep_guide, drop = TRUE)
  } else {
    dat <- prep_omega(df) %>%
      mutate(x = gen)

    p <- ggplot(dat) +
      geom_line(aes(x = x, y = y),
                linewidth = 0.35, alpha = 0.9)
  }

  p <- p +
    scale_x_continuous(limits = x_lims, breaks = x_breaks, labels = x_labels) +
    scale_y_log10(breaks = y_breaks,
                  labels = trans_format("log10", math_format(10^.x)),
                  minor_breaks = NULL) +
    coord_cartesian(xlim = x_lims, ylim = y_lims) +
    mytheme

  if (show_ylab && !is.null(ylab)) p <- p + labs(y = ylab)
  else p <- p + labs(y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  if (show_xlab) p <- p + labs(x = NULL)
  else p <- p + labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

#' Build the full 6xN time series figure for one fitness model
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
#' @param tag_prefix  When set (e.g. "zoom"), loads all replicates matching
#'   this tag prefix and overlays them as colored semi-transparent lines.
#'   Overrides tag_filter.
#' @param highlight_rep  In replicate mode, draw this rep in black over the
#'   other reps in grey. NULL = coloured overlay.
#' @param show_omega  FALSE drops the two omega rows (v, c, W_P, W_H only).
#' @param window   c(start, end) generations to plot, e.g. c(0, 2000).
#' @param results_root  Where to discover runs, e.g. "results_zoom" for the
#'   write_every = 1 zoom runs (pass tag_prefix = "zoom" with it).
fig_timeseries <- function(model_name = "acute", filename = NULL,
                           sigma = NULL, diploid = NULL,
                           max_pts = 2000, smooth = NULL,
                           step = FALSE,
                           width = NULL, height = NULL,
                           x_lims = NULL, x_breaks = NULL, x_labels = NULL,
                           tag_filter = NA, tag_prefix = NULL,
                           reps = NULL, rep_legend = FALSE, replicates = FALSE,
                           conditions = NULL,
                           highlight_rep = NULL, show_omega = TRUE,
                           window = NULL,
                           results_root = getOption("ggt.results_root", "results"),
                           show_title = TRUE) {

  cond_names  <- c("EThost_ETpath", "EThost_ERpath", "ERhost_ETpath", "ERhost_ERpath")
  col_titles  <- c("ET / ET", "ET host / ER path",
                    "ER host / ET path", "ER / ER")
  # Backward-compat: old scenario names used by load_sim
  scenario_map <- c("EThost_ETpath" = "ET-ET", "EThost_ERpath" = "ERpath-EThost",
                     "ERhost_ETpath" = "ERhost-ETpath", "ERhost_ERpath" = "ER-ER")

  # Filter to requested conditions
  if (!is.null(conditions)) {
    keep <- cond_names %in% conditions | scenario_map %in% conditions
    cond_names <- cond_names[keep]
    col_titles <- col_titles[keep]
  }

  # Overlay replicates when the caller asks for it in any of three ways:
  #   - tag_prefix set    (tagged reps, e.g. the zoom runs)
  #   - reps = c(...)      (explicit rep subset)
  #   - replicates = TRUE  (untagged reps identified by config 'rep', e.g. acute)
  has_reps <- !is.null(tag_prefix) || !is.null(reps) || isTRUE(replicates)

  tag <- model_name
  if (!is.null(diploid) && diploid) tag <- paste0(tag, " (diploid)")
  if (!is.null(sigma)) tag <- paste0(tag, " \u03c3=", sigma)
  if (!is.null(window)) tag <- paste0(tag, ", gen ", window[1], "\u2013", window[2])
  cat("\n  Loading time series for:", tag,
      if (has_reps) paste0(" [replicates: ",
                           if (!is.null(tag_prefix)) paste0(tag_prefix, "*") else "untagged",
                           "]") else "", "\n")

  if (has_reps) {
    # --- Replicate mode: load all reps via load_replicates ---
    all_rep_data <- load_replicates(
      model_name, sigma = sigma, diploid = diploid,
      conditions = cond_names, tag_prefix = tag_prefix, reps = reps,
      results_root = results_root
    )
    if (nrow(all_rep_data) == 0) {
      warning("No replicate data loaded"); return(NULL)
    }

    dfs <- setNames(
      lapply(cond_names, function(cn) {
        d <- all_rep_data %>% filter(condition == cn)
        if (!is.null(window)) d <- d %>% filter(gen >= window[1], gen <= window[2])
        if (nrow(d) == 0) return(NULL)
        # Thin per replicate to keep overlay readable
        d %>%
          group_by(rep) %>%
          group_modify(~thin_for_plot(.x, max_pts = max_pts)) %>%
          ungroup()
      }),
      cond_names
    )
  } else {
    # --- Single-run mode: load via load_sim (backward compatible) ---
    dfs <- setNames(
      lapply(cond_names, function(cn) {
        sc <- scenario_map[cn]
        d <- suppressWarnings(load_sim(model_name, sc, sigma = sigma,
                                       diploid_filter = diploid,
                                       tag_filter = tag_filter))
        if (is.null(d) || nrow(d) == 0) return(NULL)
        if (!is.null(window)) d <- d %>% filter(gen >= window[1], gen <= window[2])
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
      cond_names
    )
  }

  # Which conditions actually loaded?
  available <- !vapply(dfs, is.null, logical(1))
  if (sum(available) == 0) {
    warning("No data loaded for any condition")
    return(NULL)
  }

  active_conds  <- cond_names[available]
  active_titles <- col_titles[available]
  active_dfs    <- dfs[available]
  n_cols        <- length(active_conds)

  if (sum(available) < length(cond_names)) {
    cat("  Note: only", sum(available), "of", length(cond_names), "conditions available:",
        paste(active_titles, collapse = ", "), "\n")
  }
  if (has_reps) {
    n_reps <- length(unique(all_rep_data$rep))
    cat("  Overlaying", n_reps, "replicates per condition\n")
  }

  # Auto-detect shared time axis from all data (or the window), or use overrides
  all_gens <- if (!is.null(window)) window else unlist(lapply(active_dfs, function(d) d$gen))
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
    is_last <- !show_omega && ri == length(rows)
    for (ci in seq_along(active_conds)) {
      p <- line_panel(
        active_dfs[[ci]], row$var,
        ylab = row$ylab,
        show_ylab = (ci == 1),
        show_xlab = is_last,
        model_name = model_name,
        x_lims = tax$lims, x_breaks = tax$breaks, x_labels = tax$labels,
        use_step = step,
        has_reps = has_reps,
        show_legend = (rep_legend && has_reps),
        highlight_rep = highlight_rep
      )
      if (is_last && ci == 1) p <- p + labs(x = "Substitutions")
      # Column title on first row
      if (ri == 1) {
        p <- p + labs(title = active_titles[ci]) +
          theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
      }
      panels[[length(panels) + 1]] <- p
    }
  }

  # Omega rows
  for (who in if (show_omega) c("Path", "Host") else character(0)) {
    ylab <- if (who == "Path") expression(omega[P]) else expression(omega[H])
    is_last <- (who == "Host")
    for (ci in seq_along(active_conds)) {
      p <- omega_panel(
        active_dfs[[ci]], who,
        ylab = ylab,
        show_ylab = (ci == 1),
        show_xlab = is_last,
        x_lims = tax$lims, x_breaks = tax$breaks, x_labels = tax$labels,
        has_reps = has_reps,
        show_legend = (rep_legend && has_reps),
        highlight_rep = highlight_rep
      )
      # X-axis title only on bottom-left panel
      if (is_last && ci == 1) {
        p <- p + labs(x = "Substitutions")
      }
      panels[[length(panels) + 1]] <- p
    }
  }

  total <- wrap_plots(panels, ncol = n_cols, byrow = TRUE)
  # Collect the per-rep colour legend into a single shared legend
  if (rep_legend && has_reps) {
    total <- total + plot_layout(guides = "collect")
  }
  total <- total +
    plot_annotation(
      title = if (show_title) tag else NULL,   # off for manuscript figures
      tag_levels = "A"
    ) &
    theme(plot.tag.position = "topleft",
          plot.tag = element_text(face = "bold", size = 12))

  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 2.5 * n_cols + 1
    h <- if (!is.null(height)) height else 2.2 * (length(rows) + 2 * show_omega)
    # safe_ggsave: renders outside iCloud first (plain ggsave hits "write failed")
    safe_ggsave(paste0("figures/", filename, ".pdf"), total, width = w, height = h)
    safe_ggsave(paste0("figures/", filename, ".png"), total, width = w, height = h)
    cat("  Saved:", filename, "\n")
  }
  total
}

#' Omega time series, pathogen (top) and host (bottom) x four scenarios, all
#' replicates in grey with one highlighted in black, thinned per replicate.
#' normalize = TRUE plots omega / 2N with the cap (1) and neutral (1/2N) marked.
fig_omega_timeseries <- function(model_name = "minimal", sigma = 0.01,
                                 diploid = TRUE, tag_prefix = NULL,
                                 normalize = TRUE, highlight_rep = 1,
                                 max_pts = 400) {
  cond_names <- c("EThost_ETpath", "EThost_ERpath", "ERhost_ETpath", "ERhost_ERpath")
  col_titles <- c("ET / ET", "ET host / ER path", "ER host / ET path", "ER / ER")
  reps <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                          conditions = cond_names, tag_prefix = tag_prefix)
  if (is.null(reps) || nrow(reps) == 0) stop("No replicate data for omega time series")
  tax <- auto_time_axis(reps)
  reps_ts <- reps %>%
    group_by(condition, rep) %>%
    group_modify(~ thin_for_plot(.x, max_pts = max_pts)) %>%
    ungroup()

  panels <- list()
  for (who in c("Path", "Host")) {
    two_n <- 2 * first(reps[[if (who == "Path") "path_pop" else "host_pop"]])
    for (ci in seq_along(cond_names)) {
      ylab <- if (who == "Path") {
        if (normalize) expression(omega[P] / 2 * N[P]) else expression(omega[P])
      } else {
        if (normalize) expression(omega[H] / 2 * N[H]) else expression(omega[H])
      }
      p <- omega_panel(reps_ts %>% filter(condition == cond_names[ci]), who,
                       ylab = ylab, show_ylab = (ci == 1),
                       show_xlab = (who == "Host"),
                       x_lims = tax$lims, x_breaks = tax$breaks, x_labels = tax$labels,
                       has_reps = TRUE, highlight_rep = highlight_rep,
                       normalize = normalize) +
        # neutral omega = 1 (1/2N once normalised); the cap 2N becomes 1
        geom_hline(yintercept = if (normalize) 1 / two_n else 1,
                   linetype = "dashed", colour = "grey40")
      if (normalize)
        p <- p + geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40")
      if (who == "Path")
        p <- p + labs(title = col_titles[ci]) +
          theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
      if (who == "Host" && ci == 1) p <- p + labs(x = "Substitutions")
      panels[[length(panels) + 1]] <- p
    }
  }
  wrap_plots(panels, ncol = 4, byrow = TRUE)
}

#' Fig 4C: a zoomed window of one replicate, with host omega aligned under the
#' quantity that explains it.
#'   ER host / ET path: clearance c stuck at 0 -> omega_H at its cap
#'     (supply-limited: beneficial mutations everywhere, few substitutions)
#'   ET host / ER path: pathogen slope m_v keeps moving -> the host's landscape
#'     is non-stationary, so omega_H stays elevated while W_H is fine.
#' Rows: explaining quantity, host fitness W_H, omega_H / 2N_H.
fig_aligned_zoom <- function(model_name = "minimal", sigma = 0.01, diploid = TRUE,
                             results_root = getOption("ggt.results_root", "results"),
                             tag_prefix = NULL, rep = 1, window = c(0, 2000)) {
  cols <- list(
    list(cond = "ERhost_ETpath", title = "ER host / ET path",
         var = "s",  lab = expression(italic(c)), pseudo_log = FALSE),
    # m_v has rare spikes of several hundred: a pseudo-log axis keeps the
    # everyday wander near 0 visible alongside them
    list(cond = "EThost_ERpath", title = "ET host / ER path",
         var = "mV", lab = expression(m[v]), pseudo_log = TRUE))

  d <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                       conditions = vapply(cols, `[[`, "", "cond"),
                       tag_prefix = tag_prefix, reps = rep,
                       results_root = results_root) %>%
    filter(gen >= window[1], gen <= window[2]) %>%
    mutate(omega_H = suppressWarnings(as.numeric(omegaHost)) / (2 * host_pop))
  if (nrow(d) == 0) stop("No data in window for the aligned zoom panel")
  dt <- stats::median(diff(sort(unique(d$gen))))
  cat(sprintf("  aligned zoom: rep %d, substitutions %d-%d, one row every %g substitutions (%s)\n",
              rep, window[1], window[2], dt, results_root))

  xs <- scale_x_continuous(limits = window, expand = c(0.01, 0))
  panel <- function(dd, y, ylab, show_x, logy = FALSE, hlines = NULL) {
    p <- ggplot(dd, aes(gen, .data[[y]])) +
      geom_line(linewidth = 0.35, colour = "black") + xs + mytheme +
      labs(x = if (show_x) "Substitutions" else NULL, y = ylab)
    if (!is.null(hlines))
      p <- p + geom_hline(yintercept = hlines, linetype = c("dotted", "dashed")[seq_along(hlines)],
                          colour = "grey40")
    if (logy) p <- p + scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                                     limits = c(1e-6, 1.5))
    if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p
  }

  panels <- list()
  for (col in cols) {
    dd <- d %>% filter(condition == col$cond)
    two_n <- 2 * first(dd$host_pop)
    p1 <- panel(dd, col$var, col$lab, FALSE) + labs(title = col$title) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    if (col$pseudo_log)
      p1 <- p1 + scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                                    breaks = c(-100, -10, 0, 10, 100))
    p2 <- panel(dd, "hostFit", expression(W[H]), FALSE)
    p3 <- panel(dd, "omega_H", expression(omega[H] / 2 * N[H]), TRUE, logy = TRUE,
                hlines = c(1, 1 / two_n))   # cap (dotted), neutral (dashed)
    panels <- c(panels, list(p1, p2, p3))
  }
  # Column-major: column 1 = ER host / ET path, column 2 = ET host / ER path
  wrap_plots(panels, ncol = 2, byrow = FALSE)
}

# ============================================================================
# §6  MECHANISM OF DESTABILISATION -- Fig 3, S11
# ============================================================================


#' Neutral drift per substitution, for reactive (ER) players only.
#'
#' For each substitution made by an ER player:
#'   genotype change = RMS change of that player's reaction norm, trait(x) =
#'     b + m x, over the opponent's whole trait range x in [0, 1]:
#'     sqrt(db^2 + db*dm + dm^2/3). This is in trait units, so it can be
#'     compared with the phenotype change (intercepts and slopes alone are not);
#'   phenotype change = |change in the player's own realised trait|.
#' Neutral = the norm moved by more than one mutational step (sigma) while the
#' trait moved less than a tenth as much (ratio > ratio_threshold).
#' ET players are excluded: their genotype IS their phenotype, so the ratio is 1
#' by construction and a neutral step is impossible.
#' near_neutral flags substitutions the simulator itself treats as effectively
#' neutral (|s| N <= NEUTRAL_THRESH = 0.01); the ratio panel drops them, since
#' their near-zero phenotype steps would otherwise inflate its upper tail.
neutral_events <- function(pairs, sigma = 0.01, ratio_threshold = 10,
                           neutral_thresh = 0.01) {
  rms <- function(db, dm) sqrt(pmax(db^2 + db * dm + dm^2 / 3, 0))
  bind_rows(
    pairs %>% filter(mutator == "host", grepl("^ER", scenario)) %>%
      transmute(scenario, rep, player = "host", geno = rms(dbS, dmS),
                pheno = abs(c1 - c0), sN = abs(s_coef) * host_pop),
    pairs %>% filter(mutator == "path", grepl("ER path$|^ER / ER$", scenario)) %>%
      transmute(scenario, rep, player = "pathogen", geno = rms(dbV, dmV),
                pheno = abs(v1 - v0), sN = abs(s_coef) * path_pop)
  ) %>%
    mutate(ratio        = geno / pmax(pheno, 1e-12),
           is_neutral   = geno > sigma & ratio > ratio_threshold,
           near_neutral = !is.na(sN) & sN <= neutral_thresh)
}

load_neutral_events <- function(model_name, sigma = 0.01, diploid = TRUE,
                                tag_prefix = NULL) {
  d <- load_all_conditions(model_name, sigma, diploid, tag_prefix = tag_prefix,
                           keep_pre = TRUE)
  neutral_events(event_pairs(d), sigma = sigma)
}

# Label for the ET/ET slot, which has no reactive player to analyse
.no_er_label <- function(y) {
  annotate("text", x = "ET / ET", y = y, label = "no ER\nplayer",
           size = 5, colour = "grey45", lineheight = 0.9)
}

#' Fig 3A: fraction of an ER player's own substitutions that are neutral.
fig_neutral_fraction <- function(events) {
  frac <- events %>%
    group_by(scenario, player) %>%
    summarise(frac = mean(is_neutral), .groups = "drop")
  ggplot(frac, aes(x = scenario, y = frac, fill = player)) +
    geom_col(position = position_dodge(width = 0.75, preserve = "single"),
             width = 0.7, colour = "grey20", linewidth = 0.2) +
    .no_er_label(0.5) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c(host = "grey35", pathogen = "grey80"),
                      name = NULL) +
    labs(x = NULL, y = "Neutral substitutions") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
}

#' Fig 3B: genotype / phenotype change ratio of ER players' own substitutions,
#' excluding effectively neutral ones and steps with no phenotype change at all.
fig_decoupling <- function(events, ratio_threshold = 10) {
  dat <- events %>% filter(!near_neutral, pheno > 1e-9)
  dropped <- events %>% group_by(scenario, player) %>%
    summarise(near_neutral = mean(near_neutral),
              no_pheno = mean(!near_neutral & pheno <= 1e-9), .groups = "drop")
  cat("  decoupling panel: fraction of substitutions dropped\n")
  for (i in seq_len(nrow(dropped)))
    cat(sprintf("    %-18s %-9s near-neutral %.1f%%, zero phenotype change %.1f%%\n",
                dropped$scenario[i], dropped$player[i],
                100 * dropped$near_neutral[i], 100 * dropped$no_pheno[i]))
  dodge <- position_dodge(width = 0.85, preserve = "single")
  ggplot(dat, aes(x = scenario, y = ratio, fill = player)) +
    geom_violin(colour = "grey35", scale = "width", position = dodge) +
    geom_boxplot(aes(group = interaction(scenario, player)), width = 0.1,
                 outlier.size = 0.3, fill = "white", position = dodge,
                 show.legend = FALSE) +
    geom_hline(yintercept = ratio_threshold, linetype = "dashed", colour = "grey40") +
    .no_er_label(10^stats::quantile(log10(dat$ratio), 0.97)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_fill_manual(values = c(host = "grey55", pathogen = "grey90"),
                      name = NULL) +
    labs(x = NULL, y = "Genotype / phenotype change") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
}

calc_step_sizes <- function(df) {
  # Preferred: the realised step of ONE substitution, read straight off the
  # pre/post pair the writer emits for each recorded generation (post is the
  # state after exactly one event). Needs keep_pre = TRUE at load time.
  if ("event" %in% names(df) && any(df$event == "pre")) {
    g <- intersect(c("rep", "gen"), names(df))
    return(
      df %>%
        filter(event %in% c("pre", "post")) %>%
        group_by(across(all_of(g))) %>%
        filter(n() == 2) %>%
        summarise(delta_v = abs(v[event == "post"] - v[event == "pre"]),
                  delta_s = abs(s[event == "post"] - s[event == "pre"]),
                  mutator = mutator[event == "post"],   # who substituted
                  .groups = "drop")
    )
  }

  # Fallback (post rows only): change across the whole recording interval.
  # Differences must never span two runs: pooled replicates share gen values,
  # so sorting by gen alone interleaves independent lineages
  df %>%
    filter(event == "post") %>%
    mutate(.run = if ("rep" %in% names(df)) rep else 0L) %>%
    arrange(.run, gen) %>%
    mutate(
      delta_v = ifelse(.run == lag(.run), abs(v - lag(v)), NA_real_),
      delta_s = ifelse(.run == lag(.run), abs(s - lag(s)), NA_real_)
    ) %>%
    select(-.run) %>%
    filter(!is.na(delta_v))
}

#' Figure: Mutational step-size distributions across conditions
fig_step_sizes <- function(model_name = "acute",
                           sigma = 0.1, diploid = NULL,
                           include_pinned = FALSE,
                           width = NULL, height = NULL,
                           filename = "Figure_step_sizes",
                           tag_filter = NA, tag_prefix = NULL) {

  # keep_pre: per-event steps come from the pre/post pair of one generation
  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned,
                                tag_filter = tag_filter, tag_prefix = tag_prefix,
                                keep_pre = TRUE)
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
                          "abs(Delta*v)", "abs(Delta*c)"))   # plotmath, parsed in facets

  # Split per-event steps by who substituted: "own" when the trait's player
  # mutated, "opponent's" when the other player did and a reaction norm moved
  # the trait in response. Pooling the two hides the shape of each.
  has_source <- "mutator" %in% names(steps)
  steps <- steps %>%
    mutate(src = if (has_source)
                   ifelse((trait == "abs(Delta*c)" & mutator == "host") |
                          (trait == "abs(Delta*v)" & mutator == "path"),
                          "own substitution", "opponent's substitution")
                 else "all",
           src = factor(src, c("own substitution", "opponent's substitution", "all")))

  # A fixed-trait (ET) player should not respond to its opponent at all, but ET
  # is implemented with a residual reaction-norm slope of ~1e-4, so its trait
  # shifts by ~1e-6 on every opponent substitution. Drop those structural
  # near-zeros rather than plot them as induced steps.
  if (has_source) {
    host_et <- grepl("^ET", steps$scenario)
    path_et <- grepl("ET path$|^ET / ET$", steps$scenario)
    steps <- steps %>%
      filter(!(src == "opponent's substitution" &
               ((trait == "abs(Delta*c)" & host_et) |
                (trait == "abs(Delta*v)" & path_et))))
  }

  # A log axis cannot show events that left the trait exactly where it was, so
  # drop them explicitly and report how many: what is plotted is the step size
  # GIVEN that the trait moved.
  zero_frac <- steps %>%
    group_by(scenario, trait, src) %>%
    summarise(zero = mean(step == 0), .groups = "drop") %>%
    filter(zero > 0)
  if (nrow(zero_frac) > 0) {
    cat("  fraction of events with no change in the trait (omitted from log axis):\n")
    for (i in seq_len(nrow(zero_frac)))
      cat(sprintf("    %-18s %-2s %-24s %.0f%%\n", zero_frac$scenario[i],
                  gsub("abs\\(Delta\\*|\\)", "", zero_frac$trait[i]),
                  as.character(zero_frac$src[i]), 100 * zero_frac$zero[i]))
  }
  steps <- steps %>% filter(step > 0)
  
  if (include_pinned) {
    steps <- steps %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    steps <- steps %>% mutate(x_label = scenario)
  }
  
  # With pre/post rows each step is one substitution; otherwise it is the
  # change accumulated across the whole recording interval
  per_event <- "event" %in% names(all_df) && any(all_df$event == "pre")
  dt <- stats::median(diff(sort(unique(all_df$gen))))
  y_lab <- if (per_event) "Step size per substitution"
           else paste0("Change per ", format(dt), " substitutions")

  # Greyscale: the x axis already names the scenario, so fill carries nothing
  dodge <- position_dodge(width = 0.85)
  p <- ggplot(steps, aes(x = x_label, y = step, fill = src)) +
    geom_violin(colour = "grey35", scale = "width", position = dodge) +
    geom_boxplot(aes(group = interaction(x_label, src)), width = 0.1,
                 outlier.size = 0.3, fill = "white", position = dodge,
                 show.legend = FALSE) +
    scale_fill_manual(values = c("own substitution" = "grey55",
                                 "opponent's substitution" = "grey90",
                                 "all" = "grey80"), name = NULL) +
    facet_wrap(~ trait, scales = "free_y", labeller = label_parsed) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = NULL, y = y_lab,
         title = paste0(model_name, " — mutational step sizes")) +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = if (has_source) "top" else "none")
  
  if (is.null(sigma) && n_distinct(steps$sigma_label) > 1) {
    p <- p + facet_grid(sigma_label ~ trait, scales = "free_y",
                        labeller = labeller(trait = label_parsed))
  } else if (!is.null(sigma)) {
    # Genetic step size, the reference for realised-step amplification
    p <- p + geom_hline(yintercept = sigma, linetype = "dashed", color = "gray40")
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

#' Figure 6B: Slope distribution — scatter of (mS, mV) from ER-ER runs,
#' with stability hyperbolas at mS·mV = ±1.
#' This maps to manuscript "Figure 2: Phase Space: Strategy Slopes."
fig_slope_distribution <- function(es_data = NULL,
                                   model_name = "acute",
                                   condition = "ERhost_ERpath", sigma = 0.1,
                                   diploid = NULL,
                                   filename = NULL,
                                   width = NULL, height = NULL,
                                   title_label = "Phase Space: Strategy Slopes",
                                   pct_inline = FALSE,
                                   lim_q = 0.01,   # axis limits: [lim_q, 1 - lim_q] quantiles
                                   tag_filter = NA, tag_prefix = NULL) {

  if (is.null(es_data)) {
    # load_replicates covers tagged and untagged runs and keeps every replicate
    # (load_sim returned a single run whenever tag_prefix was NULL)
    es_data <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                               conditions = condition, tag_filter = tag_filter,
                               tag_prefix = tag_prefix)
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
  q_mS <- quantile(thin$mS, c(lim_q, 1 - lim_q), na.rm = TRUE)
  q_mV <- quantile(thin$mV, c(lim_q, 1 - lim_q), na.rm = TRUE)
  pad <- 0.15  # 15% padding
  xlim <- q_mS + c(-1, 1) * diff(q_mS) * pad
  ylim <- q_mV + c(-1, 1) * diff(q_mV) * pad
  
  # Classify interior vs boundary
  thin <- thin %>%
    mutate(
      prod_slopes = mS * mV,
      interior = abs(prod_slopes) < 1
    )
  
  # Percentage from every post row, not the thinned subset used for plotting
  post <- es_data %>% filter(event == "post")
  pct_interior <- mean(abs(post$mS * post$mV) < 1, na.rm = TRUE) * 100
  
  # Plot unstable first, then stable on top so blue is visible
  p <- ggplot(thin, aes(x = mS, y = mV)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray70", linewidth = 0.3) +
    geom_point(data = thin %>% filter(!interior),
               aes(color = interior), size = 0.8, alpha = 0.4) +
    geom_point(data = thin %>% filter(interior),
               aes(color = interior), size = 0.8, alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "grey70", "FALSE" = "black"),
                       labels = c("TRUE" = "stable", "FALSE" = "unstable"),
                       name = NULL) +
    geom_line(data = hyp_pos %>% filter(abs(mV) < max(abs(ylim))),
              aes(mS, mV), color = "grey25", linetype = "dashed",
              linewidth = 0.8, inherit.aes = FALSE) +
    geom_line(data = hyp_neg %>% filter(abs(mV) < max(abs(ylim))),
              aes(mS, mV), color = "grey25", linetype = "dashed",
              linewidth = 0.8, inherit.aes = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = expression(m[c]), y = expression(m[v]),
         title = title_label,
         subtitle = if (pct_inline) NULL
                    else sprintf("%.0f%% of time in stable region",
                                 pct_interior)) +
    mytheme +
    theme(legend.position = "right")

  if (pct_inline) {
    p <- p + annotate("text",
                      x = xlim[1] + diff(xlim) * 0.03,
                      y = ylim[2] - diff(ylim) * 0.03,
                      label = sprintf("%.0f%% stable", pct_interior),
                      hjust = 0, vjust = 1, size = 5)
  }

  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 7
    h <- if (!is.null(height)) height else 6
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}

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
                                   filename = NULL,
                                   tag_filter = NA, tag_prefix = NULL) {
  if (is.null(es_data)) {
    # load_replicates covers tagged and untagged runs and keeps every replicate
    # (load_sim returned a single run whenever tag_prefix was NULL)
    es_data <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                               conditions = condition, tag_filter = tag_filter,
                               tag_prefix = tag_prefix)
  }
  
  vgrid <- calc_violation_grid(model_name, resolution)
  nash_pt <- find_nash(model_name)
  dom <- TRAIT_DOMAIN[[model_name]]
  
  p <- ggplot(vgrid, aes(v, s)) +
    geom_tile(aes(fill = zone), alpha = 0.7) +
    scale_fill_manual(values = c("compatible" = "#E8E8E8",
                                 "violates (>=1)" = "#FFAAAA",
                                 "violates (<=-1)" = "#AAD4FF"),
                      # plotmath, so the symbols survive the PDF device
                      breaks = c("compatible", "violates (<=-1)", "violates (>=1)"),
                      labels = expression("compatible",
                                          "violates (" <= -1 * ")",
                                          "violates (" >= 1 * ")"),
                      name = "Stability") +
    geom_point(data = nash_pt, aes(x = v, y = s),
               size = 5, color = "black", shape = 16) +
    coord_fixed(xlim = dom, ylim = dom) +
    labs(x = "v", y = "c",
         title = paste0(model_name, " — Nash stability regions")) +
    mytheme
  
  # Overlay simulation trajectory if provided
  if (!is.null(es_data) && nrow(es_data) > 0) {
    thin <- thin_for_plot(es_data)
    # One path per run, so no segments jump from one replicate to the next
    if (!"rep" %in% names(thin)) thin$rep <- 0L
    p <- p + geom_path(data = thin, aes(x = v, y = s, group = rep),
                       color = "grey20", alpha = 0.12, linewidth = 0.2)
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

#' Figure: where each trait spends its time, as one stacked bar per scenario:
#' near its Nash value, at a trait boundary, or elsewhere in the interior.
#' Merges boundary occupancy and time near Nash into one exhaustive partition
#' (a boundary state is counted as boundary even if it is also within
#' nash_radius; with the Nash values used here the two never overlap).
fig_trait_occupancy <- function(model_name = "acute",
                                sigma = 0.1, diploid = NULL,
                                nash_radius = NULL, bound_frac = 0.02,
                                weight = c("time", "substitution"),
                                width = 8, height = 5, filename = NULL,
                                tag_filter = NA, tag_prefix = NULL) {
  weight <- match.arg(weight)

  all_df <- load_all_conditions(model_name, sigma, diploid,
                                tag_filter = tag_filter, tag_prefix = tag_prefix,
                                keep_pre = (weight == "time"))
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }
  # Weighted by evolutionary time: each substitution's pre-state counts for the
  # dwell time spent in it (substitutions per unit time differ ~1e4-fold
  # between scenarios, so a per-substitution fraction is not a time fraction)
  if (weight == "time") {
    all_df <- event_pairs(all_df) %>%
      transmute(scenario, v = v0, s = c0, w = dwell)
  } else {
    all_df <- all_df %>% mutate(w = 1)
  }

  dom  <- TRAIT_DOMAIN[[model_name]]
  span <- dom[2] - dom[1]
  if (is.null(nash_radius)) nash_radius <- 0.1 * span
  thr     <- bound_frac * span
  nash_pt <- find_nash(model_name)
  zones   <- c("near Nash", "interior", "boundary")

  classify <- function(x, x_nash) factor(case_when(
    x < dom[1] + thr | x > dom[2] - thr ~ "boundary",
    abs(x - x_nash) < nash_radius       ~ "near Nash",
    TRUE                                ~ "interior"), levels = zones)

  occ <- bind_rows(
    all_df %>% transmute(scenario, w, trait = "clearance (c)", zone = classify(s, nash_pt$s)),
    all_df %>% transmute(scenario, w, trait = "virulence (v)", zone = classify(v, nash_pt$v))
  ) %>%
    count(scenario, trait, zone, wt = w, .drop = FALSE) %>%
    group_by(scenario, trait) %>%
    mutate(fraction = n / sum(n)) %>%
    ungroup()

  p <- ggplot(occ, aes(x = scenario, y = fraction, fill = zone)) +
    geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
    facet_wrap(~ trait) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = c("near Nash" = "white", "interior" = "grey70",
                                 "boundary" = "grey20"), name = NULL) +
    labs(x = NULL, y = "Time") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  if (!is.null(filename)) {
    safe_ggsave(paste0("figures/", filename, ".pdf"), p, width = width, height = height)
    safe_ggsave(paste0("figures/", filename, ".png"), p, width = width, height = height, dpi = 300)
    cat("Saved:", filename, "\n")
  }
  p
}

# ============================================================================
# §7  SELECTION AND TEMPO -- Fig 4A-B, 4D-E, S7
# ============================================================================


#' Figure: Omega (substitution rate) distributions per condition
#' Violin/boxplots of omega_P and omega_H (log scale) with reference line
#' at omega = 1 separating purifying from positive selection.
#' Panels labelled by trait: Clearance (c) and Virulence (v).
fig_omega <- function(model_name = "acute",
                      sigma = 0.1, diploid = NULL,
                      normalize = TRUE,
                      include_pinned = FALSE,
                      width = NULL, height = NULL,
                      filename = "Figure_omega",
                      tag_filter = NA, tag_prefix = NULL) {

  all_df <- load_all_conditions(model_name, sigma, diploid, include_pinned,
                                tag_filter = tag_filter, tag_prefix = tag_prefix)
  if (nrow(all_df) == 0) {
    warning("No data found"); return(invisible(NULL))
  }

  grp_cols <- "scenario"
  if (is.null(sigma)) grp_cols <- c(grp_cols, "sigma_label")
  if (include_pinned) grp_cols <- c(grp_cols, "run_type")

  normalize <- normalize && all(c("host_pop", "path_pop") %in% names(all_df))
  omega_df <- all_df %>%
    mutate(
      omegaPath = suppressWarnings(as.numeric(omegaPath)),
      omegaHost = suppressWarnings(as.numeric(omegaHost))
    ) %>%
    select(all_of(grp_cols), any_of(c("host_pop", "path_pop")), omegaPath, omegaHost) %>%
    pivot_longer(c(omegaPath, omegaHost),
                 names_to = "player", values_to = "omega") %>%
    filter(!is.na(omega), omega > 0) %>%
    mutate(
      # omega / 2N: each player's ceiling is 2N and N differs 100-fold between
      # players, so raw omega is not comparable across them
      two_n  = if (normalize) 2 * ifelse(player == "omegaPath", path_pop, host_pop) else 1,
      omega  = omega / two_n,
      player = ifelse(player == "omegaPath", "Virulence (v)", "Clearance (c)"))

  # Neutral reference omega = 1, which becomes 1/2N (player-specific) when normalised
  neutral_ref <- omega_df %>% group_by(player) %>%
    summarise(yint = 1 / first(two_n), .groups = "drop")

  if (include_pinned) {
    omega_df <- omega_df %>%
      mutate(x_label = factor(paste0(scenario, "\n", run_type),
                              levels = unique(paste0(scenario, "\n", run_type))))
  } else {
    omega_df <- omega_df %>% mutate(x_label = scenario)
  }

  p <- ggplot(omega_df, aes(x = x_label, y = omega)) +
    geom_violin(fill = "grey80", colour = "grey35", scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white") +
    geom_hline(data = neutral_ref, aes(yintercept = yint),
               linetype = "dashed", color = "gray40") +
    facet_wrap(~ player) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = NULL, y = if (normalize) expression(omega / 2 * N) else expression(omega),
         title = paste0(model_name, " — substitution rate distributions")) +
    mytheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none",
          strip.text = element_text(size = 12))

  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 8
    h <- if (!is.null(height)) height else 5
    ggsave(paste0("figures/", filename, ".pdf"), p, width = w, height = h)
    ggsave(paste0("figures/", filename, ".png"), p, width = w, height = h)
    cat("Saved:", filename, "\n")
  }
  p
}

fig_discriminator <- function(model_name = "minimal",
                              sigma = 0.01, diploid = TRUE,
                              tag_prefix = NULL,
                              filename = "Discriminator") {

  cond_names <- c("EThost_ETpath", "EThost_ERpath",
                  "ERhost_ETpath", "ERhost_ERpath")

  reps <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                          conditions = cond_names, tag_prefix = tag_prefix)
  if (is.null(reps) || nrow(reps) == 0)
    stop("No replicate data loaded for discriminator")

  # Median in log-space is robust to the long upper tail; finite/positive only.
  per_rep <- reps %>%
    mutate(
      omegaHost = suppressWarnings(as.numeric(omegaHost)),
      omegaPath = suppressWarnings(as.numeric(omegaPath))
    ) %>%
    group_by(scenario, rep) %>%
    summarise(
      # omega / 2N, so both players sit on the same 0-1 scale
      omega_H = median(omegaHost[is.finite(omegaHost) & omegaHost > 0],
                       na.rm = TRUE) / (2 * first(host_pop)),
      omega_P = median(omegaPath[is.finite(omegaPath) & omegaPath > 0],
                       na.rm = TRUE) / (2 * first(path_pop)),
      neut_H  = 1 / (2 * first(host_pop)),
      neut_P  = 1 / (2 * first(path_pop)),
      .groups = "drop"
    ) %>%
    filter(is.finite(omega_H), omega_H > 0,
           is.finite(omega_P), omega_P > 0)

  # Greyscale: scenarios separate by shape, so the panel survives B&W printing
  p <- ggplot(per_rep, aes(x = omega_H, y = omega_P, shape = scenario)) +
    # neutral omega = 1, i.e. 1/2N on this scale
    geom_vline(xintercept = unique(per_rep$neut_H), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = unique(per_rep$neut_P), linetype = "dashed", color = "grey60") +
    geom_point(size = 3, colour = "grey15", alpha = 0.9,
               position = position_jitter(width = 0.05, height = 0.05,
                                          seed = 1)) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
    labs(x = expression("median " * omega[H] / 2 * N[H]),
         y = expression("median " * omega[P] / 2 * N[P])) +
    mytheme +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "grey80",
                                       linewidth = 0.3),
      legend.key = element_blank(),
      legend.margin = margin(2, 6, 2, 6),
      legend.text = element_text(size = 11)
    )

  if (!is.null(filename)) {
    safe_ggsave(paste0("figures/", filename, ".pdf"), p,
                width = 7, height = 5)
    safe_ggsave(paste0("figures/", filename, ".png"), p,
                width = 7, height = 5, dpi = 300)
    cat("Saved:", filename, "\n")
  }
  p
}

#' Per-replicate, time-weighted summaries of the tempo runs: each player's mean
#' fitness relative to its Nash fitness, and the fraction of evolutionary time
#' clearance spends at a boundary. Fitness is recomputed from the pre-state
#' traits so each state is weighted by the dwell time spent in it.
tempo_summary <- function(root = "results_gamma", sigma = 0.01) {
  d <- suppressWarnings(load_replicates("minimal", sigma = sigma, diploid = TRUE,
                                        gamma_filter = NULL, results_root = root,
                                        keep_pre = TRUE))
  if (nrow(d) == 0) stop("No tempo runs in ", root)
  d <- drop_unfinished_runs(d)
  mod <- FITNESS_MODELS[["minimal"]]
  np  <- find_nash("minimal")
  wH_nash <- mod$fH(np$v, np$s)
  wP_nash <- mod$fP(np$v, np$s)
  dom <- TRAIT_DOMAIN[["minimal"]]
  thr <- 0.02 * (dom[2] - dom[1])

  event_pairs(d) %>%
    mutate(wH = mod$fH(v0, c0), wP = mod$fP(v0, c0),
           c_bound = c0 < dom[1] + thr | c0 > dom[2] - thr) %>%
    group_by(scenario, gamma, rep, host_pop, path_pop) %>%
    summarise(rel_wH  = weighted.mean(wH, dwell) / wH_nash,
              rel_wP  = weighted.mean(wP, dwell) / wP_nash,
              c_bound = weighted.mean(c_bound, dwell),
              .groups = "drop") %>%
    mutate(R = (1 - gamma) * path_pop / (gamma * host_pop))
}

#' Fig 4D-E: fitness relative to Nash (host, pathogen) and clearance boundary
#' occupancy, both against the tempo ratio R, for the two mixed scenarios and
#' ER/ER. Returns list(fitness, boundary, summary).
tempo_panels <- function(root = "results_gamma", sigma = 0.01) {
  s <- tempo_summary(root, sigma) %>%
    filter(scenario != "ET / ET") %>% droplevels()
  look <- list(
    scale_x_log10(breaks = c(1, 1e2, 1e4),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand = expansion(mult = 0.15)),   # keeps edge labels apart across facets
    scale_shape_manual(values = c("ET host / ER path" = 17, "ER host / ET path" = 15,
                                  "ER / ER" = 16), name = NULL),
    scale_linetype_manual(values = c("ET host / ER path" = "dotted",
                                     "ER host / ET path" = "dashed",
                                     "ER / ER" = "solid"), name = NULL),
    mytheme, theme(legend.position = "top"))
  jit  <- position_jitter(width = 0.06, height = 0, seed = 1)
  line <- stat_summary(aes(group = scenario, linetype = scenario), fun = mean,
                       geom = "line", linewidth = 0.6, colour = "black")

  fit <- s %>%
    pivot_longer(c(rel_wH, rel_wP), names_to = "player", values_to = "rel") %>%
    mutate(player = ifelse(player == "rel_wH", "host", "pathogen"))
  pD <- ggplot(fit, aes(R, rel, shape = scenario)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    line + geom_point(position = jit, size = 2.4, colour = "grey15") +
    facet_wrap(~ player) +
    labs(x = "Tempo ratio R", y = "Fitness / Nash fitness") + look

  pE <- ggplot(s, aes(R, c_bound, shape = scenario)) +
    line + geom_point(position = jit, size = 2.4, colour = "grey15") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Tempo ratio R", y = "Clearance at boundary") + look

  list(fitness = pD, boundary = pE, summary = s)
}

#' S7: tempo symmetry check at R = 1 (gamma = 0.5, N_H = N_P). The minimal model
#' is symmetric under swapping host and pathogen with c <-> v, so the two mixed
#' scenarios must produce the same distributions once relabelled as "reactive"
#' and "fixed" player. Differences point to an asymmetry in the implementation.
#' Distributions, not trajectories: the runs share seeds but not sample paths.
fig_tempo_symmetry <- function(root = "results_gamma", sigma = 0.01, gamma = 0.5) {
  d <- suppressWarnings(load_replicates("minimal", sigma = sigma, diploid = TRUE,
                                        gamma_filter = gamma, results_root = root,
                                        conditions = c("EThost_ERpath", "ERhost_ETpath")))
  if (nrow(d) == 0) stop("No gamma = ", gamma, " runs in ", root)
  d <- drop_unfinished_runs(d)
  runs <- c("ET host / ER path", "ER host / ET path, v and c swapped")
  qty  <- c(reactive_trait = "Reactive player's trait",
            fixed_trait    = "Fixed player's trait",
            reactive_fit   = "Reactive player's fitness",
            fixed_fit      = "Fixed player's fitness")
  mirror <- bind_rows(
    d %>% filter(condition == "EThost_ERpath") %>%
      transmute(run = runs[1], reactive_trait = v, fixed_trait = s,
                reactive_fit = pathFit, fixed_fit = hostFit),
    d %>% filter(condition == "ERhost_ETpath") %>%
      transmute(run = runs[2], reactive_trait = s, fixed_trait = v,
                reactive_fit = hostFit, fixed_fit = pathFit)
  ) %>%
    pivot_longer(-run, names_to = "quantity", values_to = "value") %>%
    mutate(quantity = factor(qty[quantity], levels = qty),
           run = factor(run, levels = runs))
  # KS distance as an effect size only: the samples are autocorrelated, so a
  # p-value would be meaningless
  ks <- mirror %>% group_by(quantity) %>%
    summarise(D = suppressWarnings(ks.test(value[run == runs[1]],
                                           value[run == runs[2]])$statistic),
              .groups = "drop")
  cat("  symmetry check (KS distance, 0 = identical):\n")
  for (i in seq_len(nrow(ks)))
    cat(sprintf("    %-26s D = %.3f\n", ks$quantity[i], ks$D[i]))

  ggplot(mirror, aes(value, linetype = run)) +
    stat_ecdf(linewidth = 0.7, colour = "black") +
    geom_text(data = ks, aes(x = Inf, y = 0.06, label = sprintf("KS D = %.3f", D)),
              hjust = 1.05, size = 5, inherit.aes = FALSE) +
    facet_wrap(~ quantity, scales = "free_x") +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    labs(x = NULL, y = "Cumulative fraction") +
    mytheme +
    theme(legend.position = "top",
          panel.spacing.x = unit(1.5, "lines"))   # free x scales: keep edge labels apart
}

# ============================================================================
# §8  SUPPLEMENTARY DIAGNOSTICS -- S1, S3, S4, S5
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
  
  if (!show_x) p <- strip_x(p) else p <- p + labs(x = "v")
  if (!show_y) p <- strip_y(p) else p <- p + labs(y = "c")
  p
}

#' 2D hex density of (hostFit, pathFit) per scenario.
#' Mirrors hex_landscape_panel for the fitness plane.
#'   - hex density of realised fitness pairs
#'   - black dot at Nash fitness (W_H(v*, s*), W_P(v*, s*))
#'   - gold cross at empirical mean
hex_fitness_panel <- function(data, model_name = "acute",
                              nbins = 100, count_limits = NULL,
                              show_x = TRUE, show_y = TRUE,
                              show_nash = TRUE,
                              x_lims = NULL, y_lims = NULL,
                              show_legend = FALSE) {
  mod     <- FITNESS_MODELS[[model_name]]
  nash_pt <- if (show_nash) find_nash(model_name) else NULL

  # Restrict to finite fitness values for axis limits
  d <- data %>% filter(is.finite(hostFit), is.finite(pathFit))
  if (nrow(d) == 0) return(ggplot() + theme_void())

  if (is.null(x_lims)) {
    q_h   <- quantile(d$hostFit, 0.99, na.rm = TRUE)
    x_max <- if (q_h <= 1.05) 1 else q_h * 1.05
    x_min <- max(0, min(d$hostFit, na.rm = TRUE))
    x_lims <- c(x_min, x_max)
  }
  if (is.null(y_lims)) {
    q_p   <- quantile(d$pathFit, 0.99, na.rm = TRUE)
    y_max <- if (q_p <= 1.05) 1 else q_p * 1.05
    y_min <- max(0, min(d$pathFit, na.rm = TRUE))
    y_lims <- c(y_min, y_max)
  }

  ax_breaks_x <- pretty(x_lims, n = 3)
  ax_breaks_y <- pretty(y_lims, n = 3)

  p <- ggplot() +
    geom_hex(data = d, aes(x = hostFit, y = pathFit),
             bins = nbins, alpha = 0.85)

  fill_args <- list(option = "plasma", name = "Count",
                    trans = "log10", oob = scales::squish)
  if (!is.null(count_limits)) fill_args$limits <- count_limits
  p <- p + do.call(scale_fill_viridis_c, fill_args)

  if (!is.null(nash_pt)) {
    w_h_nash <- mod$fH(nash_pt$v, nash_pt$s)
    w_p_nash <- mod$fP(nash_pt$v, nash_pt$s)
    p <- p + geom_point(data = data.frame(hostFit = w_h_nash,
                                          pathFit = w_p_nash),
                        aes(x = hostFit, y = pathFit),
                        size = 4, color = "black", shape = 16)
  }

  mean_pt <- data.frame(hostFit = mean(d$hostFit, na.rm = TRUE),
                        pathFit = mean(d$pathFit, na.rm = TRUE))
  p <- p + geom_point(data = mean_pt,
                      aes(x = hostFit, y = pathFit),
                      size = 4, color = "#FFD700", shape = 4, stroke = 1.5)

  p <- p +
    scale_x_continuous(breaks = ax_breaks_x, limits = x_lims) +
    scale_y_continuous(breaks = ax_breaks_y, limits = y_lims) +
    coord_fixed() + mytheme + theme(legend.position = "none")

  if (show_legend) {
    p <- p + theme(
      legend.position = "inside",
      legend.position.inside = c(0.97, 0.97),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80",
                                       linewidth = 0.3),
      legend.margin = margin(4, 6, 4, 6),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.35, "cm")
    )
  }

  if (!show_x) p <- strip_x(p) else p <- p + labs(x = expression(W[H]))
  if (!show_y) p <- strip_y(p) else p <- p + labs(y = expression(W[P]))
  p
}

fig_grouped_2_dynamics <- function(model_name = "minimal",
                                   sigma = 0.01, diploid = TRUE,
                                   tag_prefix = NULL,
                                   nbins = 80,
                                   filename = "Fig2_dynamics",
                                   width = 13, height = 8) {

  cond_names <- c("EThost_ETpath", "EThost_ERpath",
                  "ERhost_ETpath", "ERhost_ERpath")
  col_titles <- c("ET / ET", "ET host / ER path",
                  "ER host / ET path", "ER / ER")

  all_rep_data <- load_replicates(
    model_name, sigma = sigma, diploid = diploid,
    conditions = cond_names, tag_prefix = tag_prefix
  )
  if (is.null(all_rep_data) || nrow(all_rep_data) == 0)
    stop("No replicate data loaded for Fig 2 dynamics")

  scen_dfs <- setNames(
    lapply(cond_names, function(cn) {
      all_rep_data %>% filter(condition == cn)
    }),
    cond_names
  )

  trait_panels <- lapply(seq_along(cond_names), function(ci) {
    d <- scen_dfs[[ci]]
    if (is.null(d) || nrow(d) == 0)
      return(ggplot() + theme_void())
    p <- hex_landscape_panel(d, model_name = model_name, nbins = nbins,
                             show_x = TRUE, show_y = (ci == 1))
    p + labs(title = col_titles[ci]) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  })

  # Shared fitness range across scenarios so panels are comparable.
  # Use the 99th percentile of pooled finite fitness as the upper bound,
  # with a small headroom; include the Nash fitness so the reference dot
  # is always visible.
  mod_for_nash <- FITNESS_MODELS[[model_name]]
  nash_pt      <- find_nash(model_name)
  w_h_nash     <- mod_for_nash$fH(nash_pt$v, nash_pt$s)
  w_p_nash     <- mod_for_nash$fP(nash_pt$v, nash_pt$s)
  pooled_fit   <- all_rep_data %>%
    filter(is.finite(hostFit), is.finite(pathFit))
  hi_h <- max(quantile(pooled_fit$hostFit, 0.99, na.rm = TRUE),
              w_h_nash) * 1.10
  hi_p <- max(quantile(pooled_fit$pathFit, 0.99, na.rm = TRUE),
              w_p_nash) * 1.10
  # Round up to the nearest 0.1 so axis breaks land on nice tenths.
  hi_shared <- ceiling(max(hi_h, hi_p) * 10) / 10
  fit_lims  <- c(0, hi_shared)

  fit_panels <- lapply(seq_along(cond_names), function(ci) {
    d <- scen_dfs[[ci]]
    if (is.null(d) || nrow(d) == 0)
      return(ggplot() + theme_void())
    hex_fitness_panel(d, model_name = model_name, nbins = nbins,
                      show_x = TRUE, show_y = (ci == 1),
                      x_lims = fit_lims, y_lims = fit_lims,
                      show_legend = (ci == 1))
  })

  panels <- c(trait_panels, fit_panels)
  out <- wrap_plots(panels, ncol = 4, byrow = TRUE) +
    plot_annotation(tag_levels = "A") &
    .tag_theme &
    theme(plot.margin = margin(14, 14, 6, 18))

  if (!is.null(filename)) {
    safe_ggsave(paste0("figures/", filename, ".pdf"), out,
                width = width, height = height)
    safe_ggsave(paste0("figures/", filename, ".png"), out,
                width = width, height = height, dpi = 300)
    cat("Saved:", filename, "\n")
  }
  out
}

#' S3: per-replicate time-series statistics, trimmed to mean, SD and spectral
#' slope (noise colour) for both traits, plus correlation length for c only:
#' for v it sits at the recording floor in every scenario. SD replaces CV,
#' which mostly tracked small means (e.g. ER host / ET path clearance).
#' Points = replicates, crossbar = mean over replicates.
fig_ts_stats_si <- function(model_name = "minimal", sigma = 0.01, diploid = TRUE,
                            tag_prefix = NULL) {
  d <- suppressWarnings(load_replicates(model_name, sigma = sigma, diploid = diploid,
                                        tag_prefix = tag_prefix))
  dt <- stats::median(diff(sort(unique(d$gen))))   # substitutions per row
  summ <- d %>%
    group_by(scenario, rep) %>%
    summarise(mean_v  = mean(v), mean_c = mean(s),
              sd_v    = sd(v),   sd_c   = sd(s),
              slope_v = calc_spectral_slope(v), slope_c = calc_spectral_slope(s),
              corr_c  = calc_correlation_length(s) * dt,
              .groups = "drop")
  metrics <- c(mean = "Mean", sd = "SD", slope = "Spectral slope",
               corr = "Correlation length")
  long <- summ %>%
    pivot_longer(-c(scenario, rep), names_to = "key", values_to = "value") %>%
    filter(is.finite(value)) %>%
    mutate(trait  = ifelse(grepl("_v$", key), "virulence (v)", "clearance (c)"),
           metric = factor(metrics[sub("_[vc]$", "", key)], levels = metrics))
  floor_note <- tibble(metric = factor("Correlation length", levels = metrics),
                       trait = "virulence (v)", scenario = "ER host / ET path",
                       value = 0, label = "at recording floor")

  ggplot(long, aes(x = scenario, y = value)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.55,
                 linewidth = 0.4, colour = "grey45") +
    geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1),
               size = 2, colour = "grey15") +
    geom_text(data = floor_note, aes(label = label), colour = "grey45", size = 5) +
    facet_grid(metric ~ trait, scales = "free_y", switch = "y") +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.placement = "outside",
          strip.text.y.left = element_text(size = 15, angle = 90))
}

make_snapshot_panel <- function(es_data, gen_num, grid,
                                title_label = "",
                                lag_gens = 100,
                                show_x = TRUE, show_y = TRUE,
                                show_x_label = FALSE, show_y_label = FALSE,
                                model_name = "acute") {
  
  snapshot <- es_data %>% filter(gen == gen_num) %>% slice(1)
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
                          filename = "Figure5_snapshots",
                          show_title = TRUE,
                          tag_filter = NA, tag_prefix = NULL) {

  if (is.null(es_data)) {
    # Snapshots show strategy lines at specific generations — overlaying
    # replicates would be unreadable, so use first replicate only.
    if (!is.null(tag_prefix)) {
      tag_filter <- paste0(tag_prefix, "1")
      cat("  Snapshots: using first replicate (", tag_filter, ")\n")
    }
    es_data <- load_sim(model_name, condition, sigma = sigma,
                        diploid_filter = diploid, tag_filter = tag_filter)
    if (is.null(es_data) || !is.data.frame(es_data) || nrow(es_data) == 0) {
      warning("No data found for ", model_name, "/", condition,
              " (tag_filter=", tag_filter, ")")
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
      title_label = paste0(format(g, big.mark = ","), " substitutions"),
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
      title = if (show_title) paste0("Strategy snapshots (", model_name, ")") else NULL,
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

fig_strategy_evolution <- function(es_data = NULL, model_name = "acute",
                                   condition = "ERhost_ERpath", sigma = 0.1,
                                   diploid = NULL,
                                   filename = "Figure_strategy_params",
                                   width = NULL, height = NULL,
                                   tag_filter = NA, tag_prefix = NULL) {

  has_reps <- !is.null(tag_prefix)

  if (is.null(es_data)) {
    if (has_reps) {
      es_data <- load_replicates(model_name, sigma = sigma, diploid = diploid,
                                 conditions = condition, tag_prefix = tag_prefix)
    } else {
      es_data <- load_sim(model_name, condition, sigma = sigma,
                          diploid_filter = diploid, tag_filter = tag_filter)
    }
    if (is.null(es_data) || nrow(es_data) == 0) {
      warning("No data found for ", model_name, "/", condition)
      return(invisible(NULL))
    }
  }

  # Thin per replicate if needed
  if (has_reps && "rep" %in% names(es_data)) {
    thin <- es_data %>%
      group_by(rep) %>%
      mutate(row_num = row_number()) %>%
      filter(row_num %% 10 == 0) %>%
      ungroup()
  } else {
    thin <- es_data %>%
      mutate(row_num = row_number()) %>%
      filter(row_num %% 10 == 0)
  }

  # Auto-detect appropriate x-axis
  gen_range <- range(thin$gen, na.rm = TRUE)
  use_log <- gen_range[2] > 1e5

  if (use_log) {
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

  # Shared y-axis ranges across rows: intercepts share limits, slopes share limits
  bS_rng <- range(thin$bS, na.rm = TRUE)
  bV_rng <- range(thin$bV, na.rm = TRUE)
  b_lims <- range(c(bS_rng, bV_rng))
  b_pad  <- diff(b_lims) * 0.05
  b_lims <- b_lims + c(-b_pad, b_pad)

  mS_rng <- range(thin$mS, na.rm = TRUE)
  mV_rng <- range(thin$mV, na.rm = TRUE)
  m_lims <- range(c(mS_rng, mV_rng))
  m_pad  <- diff(m_lims) * 0.05
  m_lims <- m_lims + c(-m_pad, m_pad)

  # Theme: no x-axis for top row
  top_theme <- mytheme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (has_reps && "rep" %in% names(thin)) {
    n_reps <- length(unique(thin$rep))
    lw <- if (n_reps <= 3) 0.4 else 0.3
    al <- if (n_reps <= 3) 0.7 else 0.5
    rep_scale <- scale_color_manual(values = REP_COLORS, guide = "none")

    p_bS <- ggplot(thin, aes(gen, bS, color = factor(rep, levels = names(REP_COLORS)), group = rep)) +
      geom_line(linewidth = lw, alpha = al) + rep_scale +
      log_x + ylim(b_lims) + labs(x = NULL, y = "Host intercept") + top_theme

    p_mS <- ggplot(thin, aes(gen, mS, color = factor(rep, levels = names(REP_COLORS)), group = rep)) +
      geom_line(linewidth = lw, alpha = al) + rep_scale +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      log_x + ylim(m_lims) + labs(x = NULL, y = "Host slope") + top_theme

    p_bV <- ggplot(thin, aes(gen, bV, color = factor(rep, levels = names(REP_COLORS)), group = rep)) +
      geom_line(linewidth = lw, alpha = al) + rep_scale +
      log_x + ylim(b_lims) + labs(x = "Substitutions", y = "Pathogen intercept") + mytheme

    p_mV <- ggplot(thin, aes(gen, mV, color = factor(rep, levels = names(REP_COLORS)), group = rep)) +
      geom_line(linewidth = lw, alpha = al) + rep_scale +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      log_x + ylim(m_lims) + labs(x = "Substitutions", y = "Pathogen slope") + mytheme
  } else {
    p_bS <- ggplot(thin, aes(gen, bS)) +
      geom_line(color = "black", linewidth = 0.4) +
      log_x + ylim(b_lims) + labs(x = NULL, y = "Host intercept") + top_theme

    p_mS <- ggplot(thin, aes(gen, mS)) +
      geom_line(color = "black", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      log_x + ylim(m_lims) + labs(x = NULL, y = "Host slope") + top_theme

    p_bV <- ggplot(thin, aes(gen, bV)) +
      geom_line(color = "black", linewidth = 0.4) +
      log_x + ylim(b_lims) + labs(x = "Substitutions", y = "Pathogen intercept") + mytheme

    p_mV <- ggplot(thin, aes(gen, mV)) +
      geom_line(color = "black", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      log_x + ylim(m_lims) + labs(x = "Substitutions", y = "Pathogen slope") + mytheme
  }

  combined <- (p_bS | p_mS) / (p_bV | p_mV) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  
  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 7.5
    h <- if (!is.null(height)) height else 9
    ggsave(paste0("figures/", filename, ".pdf"), combined, width = w, height = h)
  }
  combined
}

# ============================================================================
# §9  TRACKING-STRENGTH MODEL -- S8, S9
# ============================================================================


#' Per-replicate destabilisation summaries for one model and condition.
#' model_name is "minimal", "acute" or tracking_model_name(k).
destab_stats <- function(model_name, condition = "ERhost_ERpath",
                         sigma = 0.01, diploid = TRUE, tag_prefix = NULL) {
  d <- suppressWarnings(load_replicates(
    model_name, sigma = sigma, diploid = diploid, conditions = condition,
    tag_prefix = tag_prefix))
  if (is.null(d) || nrow(d) == 0) return(NULL)
  k_val <- if (grepl("^tracking", model_name)) tracking_k_of(model_name)
           else if (model_name == "minimal") 0 else NA_real_
  d %>%
    group_by(rep) %>%
    summarise(n          = n(),
              v_sd       = sd(v),
              v_boundary = mean(v < 0.02 | v > 0.98),
              r          = suppressWarnings(cor(hostFit, pathFit)),
              w_h        = mean(hostFit),
              .groups = "drop") %>%
    mutate(model = model_name, condition = condition, k = k_val)
}

#' S8: tracking sweep from the tracking runs themselves (k = 0 included, so run
#' length and replicate count match across k).
#'   A  ER/ER virulence SD against k, with ET/ET as the reference
#'   B  corr(W_H, W_P) in ER/ER against k
#'   C  ER-host fitness relative to ET/ET against k: each replicate's mean W_H
#'      divided by the mean W_H of the ET/ET runs at the same k
fig_tracking_sweep_si <- function(ks = TRACKING_KS, sigma = 0.01, diploid = TRUE) {
  register_tracking_k(ks)
  conds  <- c(ERhost_ERpath = "ER / ER", EThost_ETpath = "ET / ET",
              ERhost_ETpath = "ER host / ET path")
  stats <- bind_rows(lapply(ks, function(k)
    bind_rows(lapply(names(conds), function(cnd)
      destab_stats(tracking_model_name(k), cnd, sigma, diploid)))))
  if (nrow(stats) == 0) stop("No tracking runs found")
  stats$scen <- factor(conds[stats$condition], levels = conds)
  nrep <- stats %>% group_by(k, condition) %>% summarise(n = n(), .groups = "drop")
  cat(sprintf("  tracking sweep: %d runs, %d-%d replicates per (k, scenario)\n",
              nrow(stats), min(nrep$n), max(nrep$n)))

  # Integer breaks only: 0.5-spaced labels run into each other at this width
  kx   <- scale_x_continuous(breaks = unique(floor(ks)), minor_breaks = ks)
  jit  <- position_jitter(width = 0.06, height = 0, seed = 1)
  line <- function() stat_summary(aes(group = scen, linetype = scen), fun = mean,
                                  geom = "line", linewidth = 0.6, colour = "black")
  look <- list(kx, mytheme, theme(legend.position = "top"),
               scale_linetype_manual(values = c("ER / ER" = "solid", "ET / ET" = "dashed",
                                                "ER host / ET path" = "dotted"), name = NULL),
               scale_shape_manual(values = c("ER / ER" = 16, "ET / ET" = 1,
                                             "ER host / ET path" = 17), name = NULL))

  pA <- ggplot(filter(stats, condition %in% c("ERhost_ERpath", "EThost_ETpath")) %>% droplevels(),
               aes(k, v_sd, shape = scen)) +
    line() + geom_point(position = jit, size = 2.2, colour = "grey15") +
    labs(x = "Tracking strength k", y = "Virulence SD") + look

  pB <- ggplot(filter(stats, condition == "ERhost_ERpath") %>% droplevels(),
               aes(k, r, shape = scen)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    line() + geom_point(position = jit, size = 2.2, colour = "grey15") +
    labs(x = "Tracking strength k", y = expression(corr(W[H], W[P]))) + look +
    theme(legend.position = "none")

  ref <- stats %>% filter(condition == "EThost_ETpath") %>%
    group_by(k) %>% summarise(w_ref = mean(w_h), .groups = "drop")
  rel <- stats %>% filter(condition %in% c("ERhost_ERpath", "ERhost_ETpath")) %>%
    left_join(ref, by = "k") %>% mutate(rel = w_h / w_ref) %>% droplevels()
  pC <- ggplot(rel, aes(k, rel, shape = scen)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    line() + geom_point(position = jit, size = 2.2, colour = "grey15") +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Tracking strength k", y = "Host fitness / ET-ET") + look

  pA | pB | pC
}

#' Same rows as fig_timeseries (v, c, W_P, W_H, omega_P, omega_H) but the
#' facet columns are the tracking strengths k you ran, for a single condition.
fig_tracking_timeseries_ksweep <- function(condition = "ERhost_ERpath",
                                          ks = c(0.5, 1, 2, 4),
                                          sigma = 0.01, diploid = TRUE,
                                          max_pts = 2000,
                                          replicates = FALSE, reps = NULL,
                                          rep_legend = FALSE,
                                          highlight_rep = NULL, show_omega = TRUE,
                                          show_title = TRUE,
                                          filename = "tracking_timeseries_ksweep",
                                          width = NULL, height = 9) {
  register_tracking_k(ks)

  # Overlay replicates when replicates = TRUE or reps = c(...) is given.
  # Every k, k = 0 included, is read from the tracking-model runs.
  has_reps <- isTRUE(replicates) || !is.null(reps)

  dfs <- setNames(lapply(ks, function(k) {
    if (has_reps) {
      # Every k, k = 0 included, comes from the tracking runs themselves, so
      # run length and replicate count are the same across the sweep
      d <- suppressWarnings(load_replicates(tracking_model_name(k), sigma = sigma,
                                            diploid = diploid, conditions = condition,
                                            reps = reps))
      if (is.null(d) || nrow(d) == 0) return(NULL)
      # Thin per replicate so overlays stay readable
      d %>% group_by(rep) %>%
        group_modify(~thin_for_plot(.x, max_pts = max_pts)) %>% ungroup()
    } else {
      d <- suppressWarnings(load_sim(tracking_model_name(k), condition,
                                     sigma = sigma, diploid_filter = diploid))
      if (is.null(d) || nrow(d) == 0) return(NULL)
      thin_for_plot(d, max_pts = max_pts)
    }
  }), paste0("k", ks))

  avail   <- !vapply(dfs, is.null, logical(1))
  if (!any(avail)) { warning("No tracking data for ", condition); return(NULL) }
  ks_a    <- ks[avail]; dfs_a <- dfs[avail]
  n_cols  <- length(ks_a)
  mn_a    <- tracking_model_name(ks_a[1])   # any tracking model for axis clamps

  all_gens <- unlist(lapply(dfs_a, function(d) d$gen))
  tax <- auto_time_axis(data.frame(gen = all_gens))

  rows <- list(list(var = "v",       ylab = expression(italic(v))),
               list(var = "s",       ylab = expression(italic(c))),
               list(var = "pathFit", ylab = expression(W[P])),
               list(var = "hostFit", ylab = expression(W[H])))
  panels <- list()
  for (ri in seq_along(rows)) {
    row <- rows[[ri]]
    for (ci in seq_along(ks_a)) {
      is_last <- !show_omega && ri == length(rows)
      p <- line_panel(dfs_a[[ci]], row$var, ylab = row$ylab,
                      show_ylab = (ci == 1), show_xlab = is_last,
                      model_name = mn_a,
                      x_lims = tax$lims, x_breaks = tax$breaks, x_labels = tax$labels,
                      has_reps = has_reps,
                      show_legend = (rep_legend && has_reps),
                      highlight_rep = highlight_rep)
      if (is_last && ci == 1) p <- p + labs(x = "Substitutions")
      if (ri == 1)
        p <- p + labs(title = sprintf("k = %g", ks_a[ci])) +
          theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
      panels[[length(panels) + 1]] <- p
    }
  }
  for (who in if (show_omega) c("Path", "Host") else character(0)) {
    ylab <- if (who == "Path") expression(omega[P]) else expression(omega[H])
    is_last <- (who == "Host")
    for (ci in seq_along(ks_a)) {
      p <- omega_panel(dfs_a[[ci]], who, ylab = ylab,
                       show_ylab = (ci == 1), show_xlab = is_last,
                       x_lims = tax$lims, x_breaks = tax$breaks, x_labels = tax$labels,
                       has_reps = has_reps,
                       show_legend = (rep_legend && has_reps))
      if (is_last && ci == 1) p <- p + labs(x = "Substitutions")
      panels[[length(panels) + 1]] <- p
    }
  }

  total <- wrap_plots(panels, ncol = n_cols, byrow = TRUE)
  if (rep_legend && has_reps) total <- total + plot_layout(guides = "collect")
  total <- total +
    plot_annotation(title = if (show_title) sprintf("Tracking model — %s (columns = k)", condition) else NULL,
                    tag_levels = "A") &
    theme(plot.tag.position = "topleft",
          plot.tag = element_text(face = "bold", size = 12))

  if (!is.null(filename)) {
    w <- if (!is.null(width)) width else 2.5 * n_cols + 1
    ggsave(paste0("figures/", filename, ".pdf"), total, width = w, height = height)
    ggsave(paste0("figures/", filename, ".png"), total, width = w, height = height)
    cat("Saved:", filename, "\n")
  }
  total
}
