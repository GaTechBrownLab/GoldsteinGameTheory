# ============================================================================
# si_extra_figures.R — two supplementary figures that have no figure behind
# them in the current draft:
#
#   FigS12_escalation  where the punctuated jumps come from (slope-product
#                      class, episode phase, episode length, jump size)
#   FigS13_sawtooth    the host-push / pathogen-recover cycle in EThost/ERpath,
#                      the mechanism the mechanism-tree box describes
#
# Data are precomputed into output/tables/ (see the block in the session notes);
# rerun that step after new simulations.
#
#   Rscript si_extra_figures.R
# ============================================================================

source("analysis/Plots.R")
D <- "output/tables"
CLS <- c("m_c m_v < -1", "|m_c m_v| < 1", "m_c m_v > +1")
# note the quotes around the comparisons: plotmath reads a bare "<-" as an arrow
lab_cls <- c("m_c m_v < -1"  = 'italic(m[c])~italic(m[v])~"< -1"',
             "|m_c m_v| < 1" = 'group("|",italic(m[c])~italic(m[v]),"|")~"< 1"',
             "m_c m_v > +1"  = 'italic(m[c])~italic(m[v])~"> +1"')

# ---------------------------------------------------------------------------
# S12: escalation and jumps
# ---------------------------------------------------------------------------
cl <- read.csv(file.path(D, "jump_classes.csv")) %>%
  mutate(class = factor(class, CLS)) %>%
  pivot_longer(c(pct_states, pct_jumps), names_to = "what", values_to = "pct") %>%
  mutate(what = factor(what, c("pct_states", "pct_jumps"),
                       c("all substitutions", "jumps")))
enr <- read.csv(file.path(D, "jump_classes.csv")) %>% mutate(class = factor(class, CLS))

pA <- ggplot(cl, aes(class, pct, fill = what)) +
  geom_col(position = position_dodge(0.75), width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(data = enr, aes(class, pmax(pct_states, pct_jumps) + 6,
                            label = sprintf("%.1f x", enrichment)),
            inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("all substitutions" = "grey80", "jumps" = "grey25")) +
  scale_x_discrete(labels = function(x) parse(text = lab_cls[x])) +
  labs(x = NULL, y = "Share (%)",
       subtitle = "Jumps concentrate in the escalating class") +
  mytheme + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
                  plot.subtitle = element_text(size = 15))

ph <- read.csv(file.path(D, "jump_phases.csv")) %>%
  mutate(phase = factor(phase, c("entry", "inside", "exit", "no crossing"),
                        c("entry", "inside", "exit", "no crossing")))
pB <- ggplot(ph, aes(phase, pct)) +
  geom_col(width = 0.7, fill = "grey45", colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 5) +
  expand_limits(y = max(ph$pct) * 1.15) +
  labs(x = "Episode phase of the substitution", y = "Share of jumps (%)",
       subtitle = "Jumps happen on the way in, not on the way out") +
  mytheme + theme(plot.subtitle = element_text(size = 15))

ep <- read.csv(file.path(D, "episode_lengths.csv"))
pC <- ggplot(ep, aes(length)) +
  geom_histogram(binwidth = 1, fill = "grey45", colour = "black", linewidth = 0.3) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks = c(0, 10, 100, 1000)) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30)) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 5,
           label = sprintf("median %d, mean %.1f\n90th percentile %d, n = %d",
                           median(ep$length), mean(ep$length),
                           quantile(ep$length, 0.9), nrow(ep))) +
  labs(x = expression("Episode length (substitutions with"~italic(m[c])~italic(m[v])~"> 1)"),
       y = "Count",
       subtitle = "Escalation episodes are brief") +
  mytheme + theme(plot.subtitle = element_text(size = 15))

sz <- read.csv(file.path(D, "jump_sizes.csv")) %>% mutate(class = factor(class, CLS))
pD <- ggplot(sz, aes(size, colour = class, linetype = class)) +
  stat_ecdf(linewidth = 0.9) +
  scale_colour_manual(values = c("grey60", "grey35", "black"),
                      labels = function(x) parse(text = lab_cls[x])) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"),
                        labels = function(x) parse(text = lab_cls[x])) +
  labs(x = "Jump size (largest trait change in one substitution)",
       y = "Cumulative fraction of jumps",
       subtitle = "The largest jumps are escalation jumps") +
  mytheme + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
                  legend.title = element_blank(), plot.subtitle = element_text(size = 15))

p <- ((pA | pB) / (pC | pD)) + plot_annotation(tag_levels = "A") &
  .tag_theme & theme(plot.margin = margin(16, 10, 6, 16))
safe_ggsave("output/figures/FigS12_escalation.pdf", p, width = 15, height = 11)
safe_ggsave("output/figures/FigS12_escalation.png", p, width = 15, height = 11, dpi = 200)
cat("Saved: FigS12_escalation\n")

# ---------------------------------------------------------------------------
# S13: the sawtooth in EThost/ERpath
# ---------------------------------------------------------------------------
st <- read.csv(file.path(D, "sawtooth.csv")) %>%
  mutate(who = factor(mutator, c("host", "path"),
                      c("host substitution", "pathogen substitution")),
         step = gen - min(gen))

pE <- ggplot(st, aes(step, v)) +
  geom_step(colour = "grey55", linewidth = 0.4, direction = "hv") +
  geom_point(aes(shape = who, fill = who), size = 2.6, colour = "black", stroke = 0.4) +
  scale_shape_manual(values = c("host substitution" = 25, "pathogen substitution" = 21)) +
  scale_fill_manual(values = c("host substitution" = "black", "pathogen substitution" = "white")) +
  labs(x = "Substitutions", y = expression("Realised virulence"~italic(v)),
       subtitle = "One host substitution knocks virulence down; the pathogen climbs back\n(x axis counts substitutions, not evolutionary time)") +
  mytheme + theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
                  legend.title = element_blank(), plot.subtitle = element_text(size = 15))

# Denominators matter here: 10.8% of pathogen substitutions leave v exactly
# unchanged (they only reshape the rule), so the share that raises v is 74.1%
# of all pathogen substitutions but 83.1% of those that move it. Values from the
# three long runs; the high-resolution runs give 95.4 / 74.2 / 82.8.
shares <- data.frame(
  who = factor(c("host\nsubstitutions\nthat lower v",
                 "pathogen\nsubstitutions\nthat raise v\n(of those that move it)"),
               levels = c("host\nsubstitutions\nthat lower v",
                          "pathogen\nsubstitutions\nthat raise v\n(of those that move it)")),
  pct = c(96.3, 83.1))
pF <- ggplot(shares, aes(who, pct)) +
  geom_col(width = 0.6, fill = c("black", "grey75"), colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 5.5) +
  expand_limits(y = 108) +
  annotate("text", x = 1.5, y = 20, size = 4.6, colour = "grey25",
           label = "10.8% of pathogen substitutions\nleave v unchanged: they only\nreshape the rule") +
  labs(x = NULL, y = "Share of substitutions (%)",
       subtitle = "Direction is near-deterministic") +
  mytheme + theme(plot.subtitle = element_text(size = 15))

q <- (pE | pF) + plot_layout(widths = c(2.6, 1)) + plot_annotation(tag_levels = "A") &
  .tag_theme & theme(plot.margin = margin(16, 10, 6, 16))
safe_ggsave("output/figures/FigS13_sawtooth.pdf", q, width = 16, height = 6)
safe_ggsave("output/figures/FigS13_sawtooth.png", q, width = 16, height = 6, dpi = 200)
cat("Saved: FigS13_sawtooth\n")
