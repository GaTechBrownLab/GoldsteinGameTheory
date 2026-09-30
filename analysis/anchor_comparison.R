# ============================================================================
# anchor_comparison.R — ER/ER under the two equilibrium-selection rules
#
#   "mutant"  the equilibrium best for the player that just substituted
#   "anchor"  the one the behavioural dynamics reach from the phenotype the
#             players were already expressing
#
# Matched seeds and settings; the runs come from
#   run_experiments.py ... --eq-selection {mutant,anchor} -o results_anchor
#
#   Rscript anchor_comparison.R            # figure + summary table
# ============================================================================

source("analysis/Plots.R")

# Default: the short matched runs, both rules in one tree. Otherwise pass the
# trees to compare, e.g. the production runs against the existing controls:
#   Rscript anchor_comparison.R data/main/minimal data/rule_mutant/minimal
args <- commandArgs(trailingOnly = TRUE)
ROOTS <- if (length(args)) args else "data/main/minimal"
OUT   <- if (length(args)) "output/figures/anchor_comparison_100k" else "output/figures/anchor_comparison"

runs <- unlist(lapply(ROOTS, list.dirs, recursive = FALSE))
runs <- runs[grepl("ERhost_ERpath", basename(runs)) &
             file.exists(file.path(runs, "simulation.csv"))]
if (!length(runs)) stop("no ER/ER runs in ", paste(ROOTS, collapse = ", "))

load_run <- function(d) {
  f <- file.path(d, "simulation.csv")
  # A run still in burn-in has an empty file
  if (file.size(f) < 100) { cat("  skipping (no rows yet):", basename(d), "\n"); return(NULL) }
  df <- read.csv(f)
  if (!nrow(df)) return(NULL)
  df$rule <- if (grepl("_anchor", basename(d))) "anchor" else "mutant"
  df$rep  <- as.integer(sub(".*rep([0-9]+).*", "\\1", basename(d)))
  df
}
dat <- bind_rows(lapply(runs, load_run))
dat$rule <- factor(dat$rule, levels = c("mutant", "anchor"),
                   labels = c("Best for the mutant", "Anchored on the prior state"))
post <- dat %>% filter(event == "post")

cat("\nrows per rule and replicate:\n")
print(with(post, table(rule, rep)))

# --- summary the reviewer would want ---------------------------------------
BND <- 1e-9
summ <- post %>%
  group_by(rule, rep) %>%
  summarise(
    n            = n(),
    boundary     = mean(v <= BND | v >= 1 - BND | s <= BND | s >= 1 - BND),
    prod_gt1     = mean(mS * mV > 1),
    prod_ltm1    = mean(mS * mV < -1),
    W_H          = mean(hostFit),
    W_P          = mean(pathFit),
    sd_v         = sd(v),
    sd_c         = sd(s),
    omega_H      = median(suppressWarnings(as.numeric(omegaHost)), na.rm = TRUE),
    omega_P      = median(suppressWarnings(as.numeric(omegaPath)), na.rm = TRUE),
    .groups = "drop")
cat("\nper replicate:\n"); print(as.data.frame(summ), digits = 3)
cat("\nmean over replicates:\n")
print(as.data.frame(summ %>% group_by(rule) %>%
        summarise(across(-rep, mean), .groups = "drop")), digits = 3)

# --- time series ------------------------------------------------------------
long <- post %>%
  group_by(rule, rep) %>%
  group_modify(~ thin_for_plot(.x, max_pts = 600)) %>%
  ungroup() %>%
  select(rule, rep, gen, v, s, pathFit, hostFit) %>%
  pivot_longer(c(v, s, pathFit, hostFit), names_to = "trait", values_to = "y") %>%
  mutate(trait = factor(trait, levels = c("v", "s", "pathFit", "hostFit"),
                        labels = c("Virulence~italic(v)", "Clearance~italic(c)",
                                   "italic(W)[P]", "italic(W)[H]")))

p <- ggplot(long, aes(gen, y, group = rep)) +
  geom_line(data = ~ subset(.x, rep != 1), colour = "grey75", linewidth = 0.25) +
  geom_line(data = ~ subset(.x, rep == 1), colour = "black", linewidth = 0.3) +
  facet_grid(trait ~ rule, scales = "free_y", switch = "y",
             labeller = labeller(trait = label_parsed)) +
  labs(x = "Substitutions", y = NULL) +
  mytheme +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
        panel.spacing.x = unit(0.9, "lines"))

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
safe_ggsave(paste0(OUT, ".pdf"), p, width = 11, height = 9)
safe_ggsave(paste0(OUT, ".png"), p, width = 11, height = 9, dpi = 200)
cat("\nSaved:", OUT, "\n")
