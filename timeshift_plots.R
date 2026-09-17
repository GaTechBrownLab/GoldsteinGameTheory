# =============================================================================
# timeshift_plots.R — figures for the time-shift (cross-inoculation) assay
#
# Reads what timeshift.py writes to results/timeshift/:
#   timeshift_main.csv          main runs: 4 reps, delta in records
#                               (1 record = 10 substitutions), sym + allo
#   timeshift_zoom.csv          zoom runs: write_every = 1, delta in single
#                               substitutions; sym + allo wherever >1 rep exists
#   timeshift_<which>_pairs.csv one row per cell x (rep_path, rep_host), written
#                               with --pairs; needed for lineage-paired contrasts
#
# Error bars: wherever >= 2 lineages exist, SEs use the lineage as the unit of
# replication (W_P_se_lineage for cell means; within-lineage differences from
# the pairs file for sym - allo contrasts). The block-bootstrap SE over
# focal-time blocks is only a fallback for single-lineage cells, because it
# treats autocorrelated time points from the same run as the replication unit.
#
# Visual grammar follows Gaba & Ebert (2009) TREE 24:226-232 Fig 3:
#   filled points  = host from the PAST     (delta < 0)
#   large point    = CONTEMPORARY host      (delta = 0)
#   open points    = host from the FUTURE   (delta > 0)
#
# Definitions only. Fig 5 is built by manuscript_figures.R (fig5); standalone:
#   source("timeshift_plots.R"); fig_timeshift("minimal")
# =============================================================================

library(tidyverse)
library(patchwork)

if (!exists("mytheme")) {
  mytheme <- theme_bw() +
    theme(
      axis.ticks.length = unit(0.2, "cm"),
      legend.text       = element_text(size = 12),
      axis.text         = element_text(size = 12, color = "black"),
      axis.title        = element_text(size = 13),
      panel.border      = element_rect(fill = NA, colour = "black", linewidth = 1),
      strip.text.x      = element_text(size = 12),
      strip.background  = element_blank(),
      legend.title      = element_blank()
    )
}

dir.create("figures", showWarnings = FALSE)

TS_DIR <- "results/timeshift"

COND_LABELS <- c(
  EThost_ETpath = "ET / ET",
  EThost_ERpath = "ET host / ER path",
  ERhost_ETpath = "ER host / ET path",
  ERhost_ERpath = "ER / ER"
)
COND_ORDER <- names(COND_LABELS)

label_conditions <- function(df) {
  df %>% mutate(condition = factor(condition, levels = COND_ORDER,
                                   labels = COND_LABELS[COND_ORDER]))
}

load_timeshift <- function(which = c("main", "zoom")) {
  which <- match.arg(which)
  f <- file.path(TS_DIR, paste0("timeshift_", which, ".csv"))
  if (!file.exists(f))
    stop("Missing ", f, ". Generate it with timeshift.py first.")
  read_csv(f, show_col_types = FALSE) %>% label_conditions()
}

# Returns NULL when the pairs file was not generated (run without --pairs).
load_pairs <- function(which = c("main", "zoom")) {
  which <- match.arg(which)
  f <- file.path(TS_DIR, paste0("timeshift_", which, "_pairs.csv"))
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE) %>% label_conditions()
}

# Lineage-level SE where the cell has >= 2 lineages; block-bootstrap SE otherwise.
pick_se <- function(d, response = "W_P") {
  blk <- d[[paste0(response, "_se")]]
  lin <- paste0(response, "_se_lineage")
  if (!lin %in% names(d)) return(blk)
  coalesce(suppressWarnings(as.numeric(d[[lin]])), blk)
}

# Sympatric - allopatric contrast, differenced WITHIN each pathogen lineage and
# summarised across lineages. A lineage's sympatric and allopatric means are
# positively correlated, so this SE is smaller -- and correct -- compared with
# adding the two cell SEs in quadrature.
lineage_contrast <- function(which, model, protocol,
                             settle_mode = "solver-path",
                             response = "W_P", k = NULL) {
  p <- load_pairs(which)
  if (is.null(p)) return(NULL)
  mean_col <- paste0(response, "_mean")
  p %>%
    filter(fitness_model == model, protocol == !!protocol) %>%
    { if (protocol == "rule") filter(., settle == settle_mode) else . } %>%
    { if (!is.null(k)) filter(., tracking_k == k) else . } %>%
    group_by(tracking_k, condition, delta, rep_path) %>%
    summarise(sym  = mean(.data[[mean_col]][sympatric == 1]),
              allo = mean(.data[[mean_col]][sympatric == 0]),
              .groups = "drop") %>%
    filter(!is.na(sym), !is.na(allo)) %>%
    mutate(diff = sym - allo) %>%
    group_by(tracking_k, condition, delta) %>%
    summarise(est = mean(diff), se = sd(diff) / sqrt(n()), n_lineages = n(),
              .groups = "drop")
}

# --- Gaba & Ebert point styling -------------------------------------------
# Encode past / contemporary / future in fill and size rather than colour, so
# the panel still reads in greyscale and matches the source figure.
ts_point_layers <- function(size_small = 2.2, size_big = 5) {
  list(
    geom_line(linewidth = 0.4),
    geom_point(aes(fill = tense, size = tense), shape = 21, stroke = 0.7),
    scale_fill_manual(values = c(past = "black", contemporary = "black",
                                 future = "white"), guide = "none"),
    scale_size_manual(values = c(past = size_small, contemporary = size_big,
                                 future = size_small), guide = "none")
  )
}

# --- Guard against magnifying numerical noise ------------------------------
# ET/ET is degenerate: both traits sit at the Nash point and W_P is constant to
# ~1e-5. With scales = "free_y" that residual float noise fills the panel and
# reads as structure. Force each facet to span at least `frac` of its own mean
# so genuinely flat panels look flat.
min_span_blank <- function(df, group_col = "condition", frac = 0.06) {
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(mid = mean(y, na.rm = TRUE),
              half = max(diff(range(c(y - se, y + se), na.rm = TRUE)) / 2,
                         abs(mean(y, na.rm = TRUE)) * frac / 2),
              .groups = "drop") %>%
    transmute(!!group_col := .data[[group_col]],
              delta = 0, ymin = mid - half, ymax = mid + half) %>%
    pivot_longer(c(ymin, ymax), values_to = "y") %>%
    select(-name)
}

add_tense <- function(df) {
  df %>% mutate(tense = factor(
    ifelse(delta < 0, "past", ifelse(delta > 0, "future", "contemporary")),
    levels = c("past", "contemporary", "future")))
}

# =============================================================================
# Figure: time-shift profile (the Gaba & Ebert Fig 3 analogue)
# =============================================================================
fig_timeshift_profile <- function(which = "zoom",
                                  model = "minimal",
                                  protocol = "rule",
                                  settle_mode = "solver-path",
                                  response = c("W_P", "W_H"),
                                  sympatric_only = TRUE,
                                  delta_range = NULL,
                                  filename = NULL,
                                  width = 11, height = 3.2) {
  response <- match.arg(response)
  d <- load_timeshift(which) %>%
    filter(fitness_model == model, protocol == !!protocol) %>%
    { if (protocol == "rule") filter(., settle == settle_mode) else . } %>%
    { if (sympatric_only) filter(., sympatric == 1) else . }
  if (!is.null(delta_range))
    d <- d %>% filter(delta >= delta_range[1], delta <= delta_range[2])
  if (nrow(d) == 0) stop("No rows matched.")

  d$y  <- d[[paste0(response, "_mean")]]
  d$se <- pick_se(d, response)
  d <- add_tense(d)

  x_lab <- if (which == "zoom") "Time shift (substitutions)"
           else "Time shift (100s of substitutions)"
  y_lab <- if (response == "W_P") "Average pathogen fitness"
           else "Average host fitness"

  ggplot(d, aes(delta, y)) +
    geom_blank(data = min_span_blank(d), aes(delta, y)) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
    geom_ribbon(aes(ymin = y - se, ymax = y + se), fill = "grey70",
                alpha = 0.35, colour = NA) +
    ts_point_layers() +
    facet_wrap(~ condition, nrow = 1, scales = "free_y") +
    labs(x = x_lab, y = y_lab) +
    mytheme
}

# =============================================================================
# Figure: sympatric vs allopatric across the delta grid
#
# The allopatric series is the control that separates lineage-specific matching
# from shared directional trends in trait level. Without it a peak at delta = 0
# cannot be attributed to coevolutionary history.
# `which = "zoom"` gives single-substitution resolution; it needs >1 zoom rep.
# =============================================================================
fig_timeshift_sympatry <- function(model = "minimal",
                                   protocol = "rule",
                                   settle_mode = "solver-path",
                                   filename = NULL,
                                   width = 11, height = 3.4,
                                   which = "main") {
  d <- load_timeshift(which) %>%
    filter(fitness_model == model, protocol == !!protocol) %>%
    { if (protocol == "rule") filter(., settle == settle_mode) else . } %>%
    mutate(Pairing = ifelse(sympatric == 1, "Sympatric (same lineage)",
                            "Allopatric (different lineage)"))
  d$se <- pick_se(d, "W_P")

  span <- d %>% mutate(y = W_P_mean) %>% min_span_blank()

  x_lab <- if (which == "zoom") "Time shift (substitutions)"
           else "Time shift (100s of substitutions)"

  ggplot(d, aes(delta, W_P_mean, colour = Pairing, fill = Pairing)) +
    geom_blank(data = span, aes(delta, y), inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
    geom_ribbon(aes(ymin = W_P_mean - se, ymax = W_P_mean + se),
                alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.6) +
    geom_point(size = if (which == "zoom") 0.9 else 1.6) +
    scale_colour_manual(values = c("Sympatric (same lineage)" = "#D95F02",
                                   "Allopatric (different lineage)" = "#7570B3")) +
    scale_fill_manual(values = c("Sympatric (same lineage)" = "#D95F02",
                                 "Allopatric (different lineage)" = "#7570B3")) +
    facet_wrap(~ condition, nrow = 1, scales = "free_y") +
    labs(x = x_lab, y = "Average pathogen fitness") +
    mytheme + theme(legend.position = "bottom")
}

# =============================================================================
# Figure: what the sympatric advantage is actually made of
#
# For the minimal model W_P = v(1-v) * (1-c) factorises as f(v) * g(c), so
#     E[W_P] = E[f] E[g] + Cov(f, g)
# exactly. Only Cov(f, g) reflects the two players being matched to each other;
# the other two terms are where each player's trait happens to sit. Splitting
# the sympatric-minus-allopatric contrast this way is what distinguishes
# "reciprocal adaptation in trait levels" from genotype-specific local
# adaptation. Defined for model = "minimal" only: timeshift.py leaves
# E_f/E_g/cov_fg blank for every other model, so any other model gets an
# explanatory placeholder instead of silently showing minimal-model bars.
# =============================================================================
fig_timeshift_decomposition <- function(protocol = "rule",
                                        settle_mode = "solver-path",
                                        delta_at = 0,
                                        filename = NULL,
                                        width = 7.5, height = 4,
                                        model = "minimal") {
  if (model != "minimal") {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, size = 4.2,
                 label = paste0("Not defined for the '", model, "' model:\n",
                                "the decomposition needs W_P = f(v) g(c),\n",
                                "which only the minimal model satisfies.")) +
        theme_void()
    )
  }

  d <- load_timeshift("main") %>%
    filter(fitness_model == "minimal", protocol == !!protocol,
           delta == delta_at) %>%
    { if (protocol == "rule") filter(., settle == settle_mode) else . } %>%
    select(condition, sympatric, E_f, E_g, cov_fg, W_P_mean) %>%
    pivot_wider(names_from = sympatric,
                values_from = c(E_f, E_g, cov_fg, W_P_mean),
                names_prefix = "s")

  parts <- d %>% transmute(
    condition,
    `Pathogen trait position` = (E_f_s1 - E_f_s0) * E_g_s0,
    `Host trait level`        = E_f_s1 * (E_g_s1 - E_g_s0),
    `Matching (Cov f,g)`      = cov_fg_s1 - cov_fg_s0,
    total                     = W_P_mean_s1 - W_P_mean_s0
  ) %>%
    pivot_longer(-c(condition, total), names_to = "component",
                 values_to = "value") %>%
    mutate(component = factor(component,
      levels = c("Pathogen trait position", "Host trait level",
                 "Matching (Cov f,g)")))

  ggplot(parts, aes(condition, value, fill = component)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "black", linewidth = 0.3) +
    geom_point(data = distinct(parts, condition, total),
               aes(condition, total), inherit.aes = FALSE,
               shape = 18, size = 3.5) +
    scale_fill_manual(values = c("Pathogen trait position" = "#E8A33D",
                                 "Host trait level"        = "#7FB3D5",
                                 "Matching (Cov f,g)"      = "#C0392B")) +
    labs(x = NULL, y = expression("Contribution to sym" - "allo " * W[P])) +
         #caption = paste("Diamond = total contrast. The Cov term is NOT an",
         #                "identified matching estimate under this protocol\n(see",
         #                "fig_timeshift_matching): re-settling makes v* and c*",
         #                "jointly determined, so Cov > 0 even between",
         #                "independent lineages.")) +
    mytheme +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 20, hjust = 1))
}

# =============================================================================
# Figure: the identified matching (genotype-specificity) estimate
#
# PLOTTED QUANTITY: sympatric minus allopatric mean W_P under the PHENOTYPE
# protocol. This is a sym - allo CONTRIBUTION, exactly the same kind of quantity
# as the three bars in fig_timeshift_decomposition, so D and E are directly
# comparable. It is NOT a raw covariance.
#
# Why it is the identified matching term: under the phenotype protocol each
# player's realised trait comes from its own trajectory regardless of partner,
# so E[f] and E[g] are bit-identical between the sympatric and allopatric sets
# (verified: difference is exactly 0). Hence
#     sym - allo  =  Cov_sym(f,g) - Cov_allo(f,g)
# to machine precision, and since allopatric pairs are independent lineages
# Cov_allo is ~1e-5, so the contrast isolates the matching component with no
# trait-level confound. Use this, not the Cov bar in the rule-protocol
# decomposition, which is contaminated by joint determination of the equilibrium.
#
# MODEL DEPENDENCE: the f/g factorisation exists only for the minimal model.
# For acute, chronic or tracking, W_P does not factor into f(v)*g(c), so the
# contrast is a general interaction (dependence) term rather than a covariance
# of those two functions. The computation is identical and valid either way --
# only the interpretation and the axis label change, which is handled below.
#
# ERROR BARS: differenced within lineage when the pairs file exists. Without it
# the fallback adds the two cell SEs in quadrature, which overstates the SE.
# =============================================================================
fig_timeshift_matching <- function(model = "minimal", delta_at = 0,
                                   filename = NULL, width = 6, height = 4) {
  lc <- lineage_contrast("main", model, "phenotype")
  d <- if (!is.null(lc) && any(lc$delta == delta_at)) {
    lc %>% filter(delta == delta_at) %>% select(condition, est, se)
  } else {
    load_timeshift("main") %>%
      filter(fitness_model == model, protocol == "phenotype",
             delta == delta_at) %>%
      select(condition, sympatric, W_P_mean, W_P_se) %>%
      pivot_wider(names_from = sympatric, values_from = c(W_P_mean, W_P_se),
                  names_prefix = "s") %>%
      transmute(condition,
                est = W_P_mean_s1 - W_P_mean_s0,
                se  = sqrt(W_P_se_s1^2 + W_P_se_s0^2))
  }

  # The plotted quantity is sym - allo, the same footing as panel D -- not a raw
  # covariance. It EQUALS Cov_sym(f,g) - Cov_allo(f,g) only for the minimal
  # model, the one case where W_P factorises into f(v)*g(c). For any other model
  # it is a general interaction term and must not be labelled a covariance.
  ylab <- if (model == "minimal") {
    expression(Delta * "Cov(" * f(v) * "," ~ g(c) * "):  sym" - "allo")
  } else {
    expression("Matching contribution to sym" - "allo " * W[P])
  }

  ggplot(d, aes(condition, est)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(fill = "#C0392B", colour = "black", linewidth = 0.3, width = 0.6) +
    geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.18) +
    labs(x = NULL, y = ylab)+
         #caption = "Phenotype protocol. Error bars are +/- 1 block-bootstrap SE.") +
    mytheme + theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# =============================================================================
# Combined figure
#
# Panel C uses the zoom runs (single-substitution resolution) whenever that
# model has cross-lineage zoom pairings; otherwise it falls back to the main
# runs. Panel titles report the lineage count so a single-lineage panel is
# never mistaken for a replicated one.
# =============================================================================
fig_timeshift <- function(model = "minimal",
                          filename = paste0("Figure_timeshift_", model),
                          width = 13, height = 14) {
  zm <- load_timeshift("zoom") %>% filter(fitness_model == model)
  n_lin <- if (nrow(zm)) max(zm$n_pairings[zm$sympatric == 1]) else 0L
  lin_txt <- sprintf("%d lineage%s", n_lin, if (n_lin == 1) "" else "s")
  sym_from <- if (any(zm$sympatric == 0)) "zoom" else "main"

  pA <- fig_timeshift_profile("zoom", model, "rule", "solver-path") +
    labs(title = paste0("A  Rule shift (live cross-infection), zoom runs, ", lin_txt))
  pB <- fig_timeshift_profile("zoom", model, "phenotype") +
    labs(title = paste0("B  Phenotype shift (induction blocked), zoom runs, ", lin_txt))
  pC <- fig_timeshift_sympatry(model, which = sym_from) +
    labs(title = paste0("C  Sympatric vs allopatric, ", sym_from, " runs"))
  pD <- fig_timeshift_decomposition(model = model) +
    labs(title = "D  Composition of the sympatric advantage (rule protocol)")
  pE <- fig_timeshift_matching(model) +
    labs(title = "E  Identified matching (phenotype protocol)")

  out <- pA / pB / pC / (pD | pE) +
    plot_layout(heights = c(1, 1, 1.25, 1.5))
  if (!is.null(filename)) {
    ggsave(paste0("figures/", filename, ".pdf"), out, width = width, height = height)
    ggsave(paste0("figures/", filename, ".png"), out, width = width, height = height,
           dpi = 200)
  }
  out
}
