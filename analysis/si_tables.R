# ============================================================================
# si_tables.R — two supplementary tables
#   Table S1  every run set: design, length, recording interval, figures
#   Table S2  per-replicate statistics for the minimal model
# Writes CSV into output/tables/.
# ============================================================================
suppressMessages({library(tidyverse); library(jsonlite)})

# ---- Table S1: run inventory ------------------------------------------------
cfgs <- list.files(c("data/main", "data/zoom", "data/tempo", "data/rule_mutant",
                     "data/random"),
                   pattern = "^config\\.json$", recursive = TRUE, full.names = TRUE)
inv <- map_dfr(cfgs, function(f) {
  c1 <- fromJSON(f)
  tibble(tree = strsplit(f, "/")[[1]][1],
         model = c1$fitness_model,
         condition = c1$condition,
         rep = c1$rep %||% 0,
         gens = c1$parameters$max_gens,
         burn = c1$parameters$burn_in_gens,
         write_every = c1$parameters$write_every,
         gamma = c1$prob_host_mutate,
         k = c1$TRACKING_K %||% NA_real_,
         NH = c1$HOST_POP_N, NP = c1$PATH_POP_N,
         rule = c1$EQ_SELECTION %||% "mutant",
         proposals = c1$PROPOSALS %||% "grid")
})
s1 <- inv %>%
  group_by(tree, model, rule, proposals, gens, burn, write_every) %>%
  summarise(scenarios = n_distinct(condition), reps = n_distinct(rep),
            sweep = case_when(n_distinct(k[!is.na(k)]) > 1 ~ paste0("k = ", paste(sort(unique(k)), collapse = ", ")),
                              n_distinct(gamma) > 1 ~ paste0("gamma = ", paste(sort(unique(gamma)), collapse = ", ")),
                              TRUE ~ "-"),
            runs = n(), .groups = "drop") %>%
  arrange(tree, model)
write_csv(s1, "output/tables/TableS1_run_inventory.csv")
cat("Table S1: run inventory\n\n")
print(as.data.frame(s1), row.names = FALSE)

# ---- Table S2: per-replicate statistics -------------------------------------
B <- 1e-9
one <- function(d) {
  raw <- read.csv(file.path(d, "simulation.csv"))
  df <- raw %>% filter(event == "post")
  cfg <- fromJSON(file.path(d, "config.json"))
  # Time-weighted occupancy uses each recorded generation's PRE state weighted by
  # the dwell on its POST row, which is the time actually spent in that state.
  # (Weighting the post state by the next recorded dwell is wrong at
  # write_every = 10: that dwell belongs to a state ten substitutions later.)
  pre <- raw %>% filter(event == "pre")
  stopifnot(nrow(pre) == nrow(df))
  t_c_bound <- weighted.mean(pre$s <= 0.02 | pre$s >= 0.98, df$dwell) * 100
  tibble(scenario = cfg$condition, rep = cfg$rep,
         escalating = 100*mean(df$mS*df$mV > 1),
         either_trait_at_bound_pct_states =
           100*mean(df$v <= B | df$v >= 1-B | df$s <= B | df$s >= 1-B),
         clearance_at_bound_pct_time = t_c_bound,
         W_H = mean(df$hostFit), W_P = mean(df$pathFit),
         SD_v = sd(df$v), SD_c = sd(df$s),
         omega_H = median(suppressWarnings(as.numeric(df$omegaHost)), na.rm = TRUE),
         omega_P = median(suppressWarnings(as.numeric(df$omegaPath)), na.rm = TRUE))
}
dirs <- list.dirs("data/main/minimal", recursive = FALSE)
s2 <- map_dfr(dirs, one) %>%
  mutate(scenario = factor(scenario, c("EThost_ETpath","EThost_ERpath","ERhost_ETpath","ERhost_ERpath"),
                           c("ET/ET","ET host / ER path","ER host / ET path","ER/ER"))) %>%
  arrange(scenario, rep)
write_csv(s2, "output/tables/TableS2_per_replicate.csv")
cat("\n\nTable S2: per-replicate statistics, minimal model, 100K substitutions each\n")
cat("  either_trait_at_bound_pct_states: % of recorded substitutions with v or c exactly at 0 or 1\n")
cat("  clearance_at_bound_pct_time:      % of evolutionary time with c within 0.02 of 0 or 1",
    "(the Fig 3F measure)\n\n")
print(as.data.frame(s2 %>% mutate(across(where(is.numeric), ~signif(.x, 3)))), row.names = FALSE)
cat("\nranges across replicates:\n")
print(as.data.frame(s2 %>% group_by(scenario) %>%
  summarise(across(c(escalating, either_trait_at_bound_pct_states,
                     clearance_at_bound_pct_time, W_H, W_P),
                   ~sprintf("%.3g - %.3g", min(.x), max(.x))), .groups = "drop")), row.names = FALSE)
