# Goldstein et al. Game Theory 
# Python implementation figure functions 
# Canan Karakoc 
# October 2025 

#======================================== Set up ===================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(pracma)
library(scales)

mytheme <- theme_bw() +
  theme(axis.ticks.length = unit(.2, "cm")) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 16, color = "black"), axis.title = element_text(size = 18)) +
  theme(plot.title = element_text(size = 18))+
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    size = 1
  )) +
  theme(strip.text.x = element_text(size = 16), strip.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    linewidth = 1
  )) +
  theme(
    axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_blank()
  ) +
  theme(
    axis.title.x = element_text(margin = margin(16, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 16, 0, 0)),
    axis.text.x = element_text(margin = margin(16, 0, 0, 0)),
    axis.text.y = element_text(margin = margin(0, 16, 0, 0))
  )

theme_set(mytheme)


#================= Theory plots - Fitness landscapes - Optimum strategies ===================

# Core Parameters
beta <- 1
d_0 <- 0.1
m_c <- 0.1
m_v <- 1
epsilon <- 1e-4

c_vals <- seq(0, 1, length.out = 300)
v_vals <- seq(0, 1, length.out = 300)
grid <- expand.grid(c = c_vals, v = v_vals)


# Acute Functions (original paper)

d_acute <- function(c, v) {
  d_0 + m_c * (((1+epsilon) * c) / ((1+epsilon) - c)) + m_v * (((1+epsilon) * v )/ ((1+epsilon) - v))
}

r_acute <- function(v) v^beta

f_H_acute <- function(c, v) {
  d <- d_acute(c, v)
  c / (d + c)
}


f_P_acute <- function(c, v) {
  d <- d_acute(c, v)
  r <- r_acute(v)
  r / (d + c)
}

# Evaluate fitness functions
grid <- grid %>%
  mutate(
    f_H_acute = mapply(f_H_acute, c, v),
    f_P_acute = mapply(f_P_acute, c, v),
    joint_fitness_acute = f_H_acute * f_P_acute,
  )

fmax_H_acute <- max(grid$f_H_acute)
fmax_P_acute <- max(grid$f_P_acute)

grid <- grid %>%
  mutate(
    f_H_acute = f_H_acute / fmax_H_acute,
    f_P_acute = f_P_acute / fmax_P_acute,
    joint_fitness_acute = f_H_acute * f_P_acute,
  )

# Optimal strategies
# Best Response Curves
v_vals <- seq(0.01, 0.99, length.out = 300)
c_vals <- seq(0.01, 0.99, length.out = 300)

opt_c_ac <- data.frame(
  v = v_vals,
  c = sapply(v_vals, function(v) optimize(function(c) -f_H_acute(c, v), c(1e-3, 0.99))$minimum)
)

# Evaluate over grid instead of optimize (more robust for flat/edge behavior)
opt_v_ac <- data.frame(
  c = c_vals,
  v = sapply(c_vals, function(c) {
    v_seq <- seq(0.01, 0.99, length.out = 200)
    f_vals <- sapply(v_seq, function(v) f_P_acute(c, v))
    v_seq[which.max(f_vals)]
  })
)

# Nash Point (best overlap of host/pathogen response)
find_nash_brute_force <- function(host_df, path_df) {
  expand.grid(i = 1:nrow(host_df), j = 1:nrow(path_df)) %>%
    mutate(
      c_host = host_df$c[i],
      v_host = host_df$v[i],
      c_path = path_df$c[j],
      v_path = path_df$v[j],
      dist = sqrt((c_host - c_path)^2 + (v_host - v_path)^2)
    ) %>%
    arrange(dist) %>%
    slice(1) %>%
    transmute(v = v_host, c = c_host)
}

# Acute
opt_joint_acute <- find_nash_brute_force(opt_c_ac, opt_v_ac)


# PLOTS 
A <- ggplot(grid, aes(x=v, y=c, z=f_H_acute)) +
  geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +  # Creates 4 bins
  geom_line(data = opt_c_ac, aes(x=v, y=c), color="steelblue", size=2, inherit.aes=FALSE) +
  labs(x="v (virulence)", y="c (clearance)") +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  scale_fill_viridis_d(option = "viridis", 
                       labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", 
                                  "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))

B <- ggplot(grid, aes(x=v, y=c, z=f_P_acute)) +
  geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
  geom_line(data = opt_v_ac, aes(x=v, y=c), color="lightcoral", size=2, inherit.aes=FALSE) +
  labs(x="v (virulence)", y="c (clearance)") +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  scale_fill_viridis_d(option = "viridis", 
                       labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", 
                                  "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))

C <- ggplot(grid, aes(x = v, y = c, z = joint_fitness_acute)) +
  geom_contour(aes(z = f_H_acute), bins = 10, color = "steelblue", alpha = 0.5) +
  geom_contour(aes(z = f_P_acute), bins = 10, color = "lightcoral", alpha = 0.5) +
  geom_line(data = opt_c_ac, aes(x = v, y = c), color = "darkblue", inherit.aes = FALSE, size = 1) +
  geom_line(data = opt_v_ac, aes(x = v, y = c), color = "firebrick", inherit.aes = FALSE, size = 1) +
  geom_point(data = opt_joint_acute, aes(x = v, y = c), color = "black", size = 3, inherit.aes = FALSE) +
  labs(x = "v (virulence)", y = "c (clearance)") +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))

# Plot as facets
strip_y <- function(p) {
  p + labs(y = NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

A1 <- A + labs(x = NULL, y = "c (clearance)")
B1 <- strip_y(B) + labs(x = "v (virulence)")
C1 <- strip_y(C) + labs(x = NULL)

p <- (A1 | B1 | C1) +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.position = "right",
    text = element_text(size = 14)
  )

final_fig <- p + plot_annotation(tag_levels = "A")

ggsave("figures/theory_figure.png", final_fig, width = 10, height = 4, units = "in")
ggsave("figures/theory_figure.pdf", final_fig, width = 10, height = 4, units = "in")

C2 <- ggplot(grid, aes(x = v, y = c, z = joint_fitness_acute)) +
  geom_contour(aes(z = f_H_acute), bins = 10, color = "steelblue", alpha = 0.5) +
  geom_contour(aes(z = f_P_acute), bins = 10, color = "lightcoral", alpha = 0.5) +
  geom_vline(xintercept = 0.5278595, size = 1, color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = 0.7773371, size = 1, color = "firebrick", linetype = "dashed") +
  geom_line(data = opt_c_ac, aes(x = v, y = c), color = "darkblue", inherit.aes = FALSE, size = 1) +
  geom_line(data = opt_v_ac, aes(x = v, y = c), color = "firebrick", inherit.aes = FALSE, size = 1) +
  geom_point(data = opt_joint_acute, aes(x = v, y = c), color = "black", size = 3, inherit.aes = FALSE) +
  labs(x = "v (virulence)", y = "c (clearance)") +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  geom_point(data = opt_joint_acute, aes(x = v, y = c), color = "black", size = 3, inherit.aes = FALSE)

ggsave("figures/optimum.png", C2, width = 4.2, height = 3.6, units = "in")
ggsave("figures/optimum.pdf", C2, width = 4.2, height = 3.6, units = "in")


#============= Multiple Nash equlibria and example for tree intersections ==============

# Core parameters
beta <- 1; d0 <- 0.1; mc <- 0.1; mv <- 1.0; epsilon <- 1e-4

clamp01 <- function(x) pmin(1, pmax(0, x))

mortality <- function(s, v) {
  onep <- 1 + epsilon
  d0 + (mc * onep * s)/(onep - s) + (mv * onep * v)/(onep - v)
}

hostFit <- function(v, s) {
  m <- mortality(s, v)
  s / (s + m)
}

pathFit <- function(v, s) {
  m <- mortality(s, v)
  (v^beta) / (s + m)
}

# Gradient for weak Nash calculation
grad <- function(f, v, s, h = 1e-5) {
  v1 <- clamp01(v - h); v2 <- clamp01(v + h)
  s1 <- clamp01(s - h); s2 <- clamp01(s + h)
  dv <- (f(v2, s) - f(v1, s)) / max(v2 - v1, 1e-12)
  ds <- (f(v, s2) - f(v, s1)) / max(s2 - s1, 1e-12)
  list(dv = dv, ds = ds)
}

# Build contour grid
v_vals <- seq(0, 1, length.out = 300)
s_vals <- seq(0, 1, length.out = 300)
grid <- expand.grid(v = v_vals, s = s_vals)
grid$fH <- mapply(hostFit, grid$v, grid$s)
grid$fP <- mapply(pathFit, grid$v, grid$s)

v_plot <- seq(0, 1, length.out = 900)
s_plot <- seq(0, 1, length.out = 900)

# ==== PANEL B: Multiple Nash Equilibria ====

# 1. First weak Nash (filled circle)
v_weak1 <- 0.33
s_weak1 <- 0.62

gH1 <- grad(hostFit, v_weak1, s_weak1)
gP1 <- grad(pathFit, v_weak1, s_weak1)

mS1 <- -gP1$dv / gP1$ds
mV1 <- -gH1$ds / gH1$dv

bS1 <- s_weak1 - mS1 * v_weak1
bV1 <- v_weak1 - mV1 * s_weak1

# Strategy lines for first weak Nash (dashed)
B_host_line1 <- data.frame(v = v_plot, s = clamp01(bS1 + mS1 * v_plot))
B_path_line1 <- data.frame(v = clamp01(bV1 + mV1 * s_plot), s = s_plot)

# 2. Second weak Nash (gray/open circle) - different location
# Try a point with different strategy slopes
v_weak2 <- 0.42
s_weak2 <- 0.70

gH2 <- grad(hostFit, v_weak2, s_weak2)
gP2 <- grad(pathFit, v_weak2, s_weak2)

mS2 <- -gH2$dv / gP2$ds
mV2 <- -gH2$ds / gH2$dv

bS2 <- s_weak2 - mS2 * v_weak2
bV2 <- v_weak2 - mV2 * s_weak2

# Strategy lines for second weak Nash (dotted)
B_host_line2 <- data.frame(v = v_plot, s = clamp01(bS2 + mS2 * v_plot))
B_path_line2 <- data.frame(v = clamp01(bV2 + mV2 * s_plot), s = s_plot)

# 3. Strong Nash (star)
best_response_host <- function(v) {
  optimize(function(s) -hostFit(v, s), c(1e-6, 1-1e-6))$minimum
}

best_response_path <- function(s) {
  optimize(function(v) -pathFit(v, s), c(1e-6, 1-1e-6))$minimum
}

v_grid <- seq(0.02, 0.98, length.out = 200)
s_grid <- seq(0.02, 0.98, length.out = 200)

BR_host <- data.frame(v = v_grid, s = sapply(v_grid, best_response_host))
BR_path <- data.frame(s = s_grid, v = sapply(s_grid, best_response_path))

distances <- expand.grid(i = 1:nrow(BR_host), j = 1:nrow(BR_path))
distances$dist <- sqrt((BR_host$v[distances$i] - BR_path$v[distances$j])^2 + 
                       (BR_host$s[distances$i] - BR_path$s[distances$j])^2)
closest <- distances[which.min(distances$dist), ]
strong_nash <- data.frame(v = BR_host$v[closest$i], s = BR_host$s[closest$i])

# Create points dataframe
nash_points <- data.frame(
  v = c(v_weak1, v_weak2, strong_nash$v),
  s = c(s_weak1, s_weak2, strong_nash$s),
  type = c("weak1", "weak2", "strong")
)

# Common axis theme for opposite-side ticks
opposite_ticks <- theme(
  axis.ticks.length = unit(3, "pt"),
  axis.ticks.y.right = element_line(),
  axis.ticks.x.top = element_line(),
  axis.text.y.right = element_text(margin = margin(l = 4)),
  axis.text.x.top = element_text(margin = margin(b = 4))
)


pB <- ggplot(grid, aes(v, s)) +
  geom_contour(aes(z = fP), color = "lightcoral", bins = 10, linewidth = 0.3, alpha = 0.85) +
  geom_contour(aes(z = fH), color = "steelblue", bins = 10, linewidth = 0.3, alpha = 0.85) +
  # First weak Nash strategy lines (dashed)
  geom_line(data = B_host_line1, aes(v, s), linetype = "dashed", linewidth = 1.0, color = "darkblue") +
  geom_line(data = B_path_line1, aes(v, s), linetype = "dashed", linewidth = 1.0, color = "firebrick") +
  # Second weak Nash strategy lines (dotted)
  geom_line(data = B_host_line2, aes(v, s), linetype = "dotted", linewidth = 1.0, color = "darkblue") +
  geom_line(data = B_path_line2, aes(v, s), linetype = "dotted", linewidth = 1.0, color = "firebrick") +
  # Equilibrium points
  geom_point(data = subset(nash_points, type == "weak1"), aes(x = v, y = s), 
             size = 5, colour = "black") +  # Weak Nash 1 (filled)
  geom_point(data = subset(nash_points, type == "weak2"), aes(x = v, y = s), 
             size = 4, shape = 21, fill = "gray70", colour = "black", stroke = 1) +  # Weak Nash 2 (gray)
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(x = "v (virulence)", y = "c (clearance)") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  mytheme+
    scale_y_continuous(
    breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
    sec.axis = dup_axis(name = NULL)  # mirror y ticks to right
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
    sec.axis = dup_axis(name = NULL)  # mirror x ticks to top
  ) +
  labs(x = NULL, y = "c (clearance)") +  # no x title here
  opposite_ticks

# ==== PANEL C: Three intersections ====

v_star <- 0.38
s_star <- 0.60
mS_C <- 2.6
mV_C <- 2.2

bS_C <- s_star - mS_C * v_star
bV_C <- v_star - mV_C * s_star

margin <- 1e-3
if (bS_C > -margin) bS_C <- -margin
if (bS_C + mS_C < 1 + margin) bS_C <- (1 + margin) - mS_C
if (bV_C > -margin) bV_C <- -margin
if (bV_C + mV_C < 1 + margin) bV_C <- (1 + margin) - mV_C

C_host_line <- data.frame(v = v_plot, s = clamp01(bS_C + mS_C * v_plot))
C_path_line <- data.frame(v = clamp01(bV_C + mV_C * s_plot), s = s_plot)

C_interior <- data.frame(v = v_star, s = s_star)
C_exterior <- data.frame(
  v = c(0, 1),
  s = c(clamp01(bS_C + mS_C * 0), clamp01(bS_C + mS_C * 1))
)

pC <- ggplot(grid, aes(v, s)) +
  geom_contour(aes(z = fP), color = "lightcoral", bins = 10, linewidth = 0.3, alpha = 0.85) +
  geom_contour(aes(z = fH), color = "steelblue", bins = 10, linewidth = 0.3, alpha = 0.85) +
  geom_line(data = C_host_line, aes(v, s), linetype = "dashed", linewidth = 1.0, color = "darkblue") +
  geom_line(data = C_path_line, aes(v, s), linetype = "dashed", linewidth = 1.0, color = "firebrick") +
  geom_point(data = C_interior, aes(v, s), size = 4, shape = 21, fill = "gray70", colour = "black", stroke = 1) +  # Unstable (filled)
  geom_point(data = C_exterior, aes(v, s), size = 5, fill = "black", colour = "black", stroke = 1) +  # Stable (open)
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(x = "v (virulence)", y = "c (clearance)") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  mytheme


strip_y <- function(p) {
  p + labs(y = NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

A2 <- pB + labs(x = "v (virulence)", y = "c (clearance)")
B2 <- strip_y(pC) + labs(x = "v (virulence)")

p2 <- (A2 | B2) +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.position = "right",
    text = element_text(size = 14)
  )

final_fig2 <- p2 + plot_annotation(tag_levels = "A")  

ggsave("figures/nash_figure5.png", final_fig2, width = 9, height = 4, units = "in")
ggsave("figures/nash_figure5.pdf", final_fig2, width = 9, height = 4, units = "in")

#================================== Simulation results ==========================================
## Trait Evolution, fitness over time, relative substitution rates

# Working directory
setwd("~/Documents/GitHub/GoldsteinGameTheory")

library(zoo)

# ---------- Omega spike panel using actual CSV omega columns ----------
omega_spike_panel <- function(df, who = c("Path", "Host"),
                              k_per_decade = 40,
                              cap = 1e6,  # Cap extreme outliers for visualization
                              show_xlab = FALSE, 
                              show_ylab = TRUE, 
                              ylab = NULL) {
  who <- match.arg(who)
  
  # Select correct omega column from CSV
  omega_col <- if (who == "Path") "omegaPath" else "omegaHost"
  
  # KEY FIX: Read omega directly from CSV, don't calculate from dwell
dat <- df %>%
  filter(event == "post") %>%
  transmute(
    x = gen,
    y = pmin(cap, pmax(1e-12, as.numeric(.data[[omega_col]]))),  # Use CSV omega
    lbin = floor(log10(x) * k_per_decade)
  ) %>%
  group_by(lbin) %>%
  slice_max(order_by = y, n = 1, with_ties = FALSE) %>%  # Keep tallest spike per bin
  ungroup()

# Spike plot (vertical lines)
p <- ggplot(dat) +
  geom_segment(aes(x = x, xend = x, y = 1e-2, yend = pmax(1e-2, y)),
               linewidth = 0.35, alpha = 0.9) +
  #scale_x_log10(limits = X_LIMS, breaks = X_BREAKS, labels = X_LABS) +
  scale_x_log10(limits = X_LIMS, breaks = X_BREAKS, labels = X_LABS)+
  scale_y_log10(limits = W_LIMS, breaks = W_BREAKS, labels = W_LABS,
                minor_breaks = NULL) +
  coord_cartesian(xlim = X_LIMS, ylim = W_LIMS) +
  mytheme

if (show_ylab) {
  p <- p + labs(y = ylab)
} else {
  p <- p + labs(y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

if (show_xlab) {
  p <- p + labs(x = "Evolutionary time")
} else {
  p <- p + labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

p
}

# ---------- Line panels for v, c, F (keep existing) ----------
make_line <- function(df, y, ylab = NULL, show_xlab = FALSE, show_ylab = TRUE) {
  p <- ggplot(df, aes(x = gen, y = {{y}})) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    scale_x_log10(limits = X_LIMS, breaks = X_BREAKS, labels = X_LABS) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    coord_cartesian(xlim = X_LIMS) +
    mytheme
  
  if (show_ylab) p <- p + labs(y = ylab) 
  else p <- p + labs(y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  if (show_xlab) p <- p + labs(x = "Evolutionary time") 
  else p <- p + labs(x = NULL) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  p
}

# ---------- Thin data for line plots (optional, speeds up plotting) ----------
thin_for_plot <- function(df_post, every = 100) {
  df_post %>%
    filter(event == "post") %>%
    mutate(.row = row_number()) %>%
    filter(.row %% every == 0)
}

# ---------- Build Figure 1 with CORRECTED omega panels ----------
# Load your data
ei_acute <- read.csv("results/ei_simulation.csv")
es_acute <- read.csv("results/es_simulation.csv")

ei_post2 <- ei_acute %>% filter(event == "post")
es_post2 <- es_acute %>% filter(event == "post")

ei_thin <- thin_for_plot(ei_post2, every = 100)
es_thin <- thin_for_plot(es_post2, every = 100)


# ---------- LOG VERSION (original manuscript style) ----------
X_LIMS_LOG   <- c(1e4, 1e6)
X_BREAKS_LOG <- c(1e4, 1e5, 1e6)
X_LABS_LOG   <- trans_format("log10", math_format(10^.x))

W_LIMS_LOG   <- c(1e-2, 1e6)
W_BREAKS_LOG <- c(1e-2, 1e2, 1e6)  # Sparser breaks
W_LABS_LOG   <- trans_format("log10", math_format(10^.x))

# Log-scale line panels
make_line_log <- function(df, y, ylab = NULL, show_xlab = FALSE, show_ylab = TRUE) {
  p <- ggplot(df, aes(x = gen, y = {{y}})) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    scale_x_log10(limits = X_LIMS_LOG, breaks = X_BREAKS_LOG, labels = X_LABS_LOG) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    coord_cartesian(xlim = X_LIMS_LOG) +
    mytheme
  
  if (show_ylab) p <- p + labs(y = ylab) 
  else p <- p + labs(y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  if (show_xlab) p <- p + labs(x = "Evolutionary time") 
  else p <- p + labs(x = NULL) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  p
}

# Log-scale omega panels
omega_panel_log <- function(df, who = c("Path", "Host"),
                            k_per_decade = 40,
                            cap = 1e6,
                            show_xlab = FALSE, 
                            show_ylab = TRUE, 
                            ylab = NULL) {
  who <- match.arg(who)
  omega_col <- if (who == "Path") "omegaPath" else "omegaHost"
  
  dat <- df %>%
    filter(event == "post") %>%
    transmute(
      x = gen,
      y = pmin(cap, pmax(1e-12, as.numeric(.data[[omega_col]]))),
      lbin = floor(log10(x) * k_per_decade)
    ) %>%
    group_by(lbin) %>%
    slice_max(order_by = y, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  p <- ggplot(dat) +
    geom_segment(aes(x = x, xend = x, y = 1e-2, yend = pmax(1e-2, y)),
                 linewidth = 0.35, alpha = 0.9) +
    scale_x_log10(limits = X_LIMS_LOG, breaks = X_BREAKS_LOG, labels = X_LABS_LOG) +
    scale_y_log10(limits = W_LIMS_LOG, breaks = W_BREAKS_LOG, labels = W_LABS_LOG,
                  minor_breaks = NULL) +
    coord_cartesian(xlim = X_LIMS_LOG, ylim = W_LIMS_LOG) +
    mytheme
  
  if (show_ylab) {
    p <- p + labs(y = ylab)
  } else {
    p <- p + labs(y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  if (show_xlab) {
    p <- p + labs(x = "Evolutionary time")
  } else {
    p <- p + labs(x = NULL) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  p
}

# Build log-scale figure
pA_log <- make_line_log(ei_post2, y = v,       ylab = expression(italic(v)))
pC_log <- make_line_log(ei_post2, y = s,       ylab = expression(italic(c)))
pE_log <- make_line_log(ei_post2, y = pathFit, ylab = expression(italic(F)[P]))
pG_log <- make_line_log(ei_post2, y = hostFit, ylab = expression(italic(F)[H]))
pI_log <- omega_panel_log(ei_post2, "Path", ylab = expression(omega[P]))
pK_log <- omega_panel_log(ei_post2, "Host", ylab = expression(omega[H]), show_xlab = TRUE)

pB_log <- make_line_log(es_post2, y = v,       show_ylab = FALSE)
pD_log <- make_line_log(es_post2, y = s,       show_ylab = FALSE)
pF_log <- make_line_log(es_post2, y = pathFit, show_ylab = FALSE)
pH_log <- make_line_log(es_post2, y = hostFit, show_ylab = FALSE)
pJ_log <- omega_panel_log(es_post2, "Path", show_ylab = FALSE)
pL_log <- omega_panel_log(es_post2, "Host", show_ylab = FALSE, show_xlab = TRUE)

final_plot_log <- (
  (pA_log | pB_log) /
    (pC_log | pD_log) /
    (pE_log | pF_log) /
    (pG_log | pH_log) /
    (pI_log | pJ_log) /
    (pK_log | pL_log)
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = "topleft",
    plot.tag = element_text(face = "bold", size = 12)
  )

print(final_plot_log)


# ---------- Summary statistics (for paper) ----------
cat("\n=== SUMMARY STATISTICS ===\n\n")

cat("EI Mode:\n")
ei_post2 %>%
  summarise(
    omega_P_median = median(omegaPath),
    omega_H_median = median(omegaHost), 
    v_median = median(v), 
    c_median = median(s), 
    fh_median = median(hostFit), 
    fp_median = median(pathFit)
  ) %>% print()


cat("\nES Mode:\n")
es_post2 %>%
  summarise(
    omega_P_median = median(omegaPath),
    omega_H_median = median(omegaHost), 
    v_median = median(v), 
    c_median = median(s), 
    fh_median = median(hostFit), 
    fp_median = median(pathFit)
  ) %>% print()

ggsave("figures/simulation.png", final_plot_log, width = 9, height = 9, units = "in")
ggsave("figures/simulation.pdf", final_plot_log, width = 9, height = 9, units = "in")

# Supplementary figures

# ---------- LINEAR X settings ----------
library(scales)

X_LIMS_LIN   <- c(0, 1e6)
X_BREAKS_LIN <- c(0, 2e5, 4e5, 6e5, 8e5, 1e6)
X_LABS_LIN   <- label_number(accuracy = 1, big.mark = ",")

# Zoom window: 500k - 1M
X_LIMS_ZOOM   <- c(5e5, 1e6)
X_BREAKS_ZOOM <- c(5e5, 6e5, 7e5, 8e5, 9e5, 1e6)
X_LABS_ZOOM   <- label_number(accuracy = 1, big.mark = ",")

# Reuse omega Y settings from your LOG version
W_LIMS_LOG   <- c(1e-2, 1e6)
W_BREAKS_LOG <- c(1e-2, 1e2, 1e6)
W_LABS_LOG   <- trans_format("log10", math_format(10^.x))

# ---------- Line panels: linear x ----------
make_line_linear <- function(df, y, x_lims, x_breaks, x_labels,
                             ylab = NULL, show_xlab = FALSE, show_ylab = TRUE) {
  p <- ggplot(df, aes(x = gen, y = {{y}})) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    scale_x_continuous(limits = x_lims, breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    coord_cartesian(xlim = x_lims) +
    mytheme
  
  if (show_ylab) p <- p + labs(y = ylab) else
    p <- p + labs(y = NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  if (show_xlab) p <- p + labs(x = "Evolutionary time") else
    p <- p + labs(x = NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

# ---------- Omega spike panels: linear x (ω still log y) ----------
# bins controls the number of linear x-bins; keep fairly high so spikes look similar to the log-binned version
omega_panel_linear <- function(df, who = c("Path", "Host"),
                               cap = 1e6,
                               x_lims, x_breaks, x_labels,
                               bins = 2000,
                               show_xlab = FALSE, show_ylab = TRUE, ylab = NULL) {
  who <- match.arg(who)
  omega_col <- if (who == "Path") "omegaPath" else "omegaHost"
  
  # Thin by tallest spike within each linear x-bin across the chosen window
  rng <- range(x_lims)
  edges <- seq(rng[1], rng[2], length.out = bins + 1)
  
  dat <- df %>%
    dplyr::filter(event == "post", gen >= rng[1], gen <= rng[2]) %>%
    dplyr::transmute(
      x = gen,
      y = pmin(cap, pmax(1e-12, as.numeric(.data[[omega_col]]))),
      lbin = cut(gen, breaks = edges, include.lowest = TRUE, right = TRUE)
    ) %>%
    dplyr::group_by(lbin) %>%
    dplyr::slice_max(order_by = y, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  p <- ggplot(dat) +
    geom_segment(aes(x = x, xend = x, y = 1e-2, yend = pmax(1e-2, y)),
                 linewidth = 0.35, alpha = 0.9) +
    scale_x_continuous(limits = x_lims, breaks = x_breaks, labels = x_labels) +
    scale_y_log10(limits = W_LIMS_LOG, breaks = W_BREAKS_LOG, labels = W_LABS_LOG,
                  minor_breaks = NULL) +
    coord_cartesian(xlim = x_lims, ylim = W_LIMS_LOG) +
    mytheme
  
  if (show_ylab) {
    p <- p + labs(y = ylab)
  } else {
    p <- p + labs(y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  if (show_xlab) {
    p <- p + labs(x = "Evolutionary time")
  } else {
    p <- p + labs(x = NULL) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  p
}

# ---------- Assemble: linear FULL (0–1e6) ----------
pA_lin <- make_line_linear(ei_post2, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(v)))
pC_lin <- make_line_linear(ei_post2, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(c)))
pE_lin <- make_line_linear(ei_post2, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(F)[P]))
pG_lin <- make_line_linear(ei_post2, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(F)[H]))
pI_lin <- omega_panel_linear(ei_post2, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             ylab = expression(omega[P]))
pK_lin <- omega_panel_linear(ei_post2, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             ylab = expression(omega[H]), show_xlab = TRUE)

pB_lin <- make_line_linear(es_post2, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pD_lin <- make_line_linear(es_post2, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pF_lin <- make_line_linear(es_post2, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pH_lin <- make_line_linear(es_post2, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pJ_lin <- omega_panel_linear(es_post2, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE)
pL_lin <- omega_panel_linear(es_post2, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE, show_xlab = TRUE)

final_plot_linear_full <- (
  (pA_lin | pB_lin) /
    (pC_lin | pD_lin) /
    (pE_lin | pF_lin) /
    (pG_lin | pH_lin) /
    (pI_lin | pJ_lin) /
    (pK_lin | pL_lin)
) + plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(face = "bold", size = 12))

print(final_plot_linear_full)

# ---------- Assemble: linear ZOOM (500k–1M) ----------
pA_z <- make_line_linear(ei_post2, y = v,       X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, ylab = expression(italic(v)))
pC_z <- make_line_linear(ei_post2, y = s,       X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, ylab = expression(italic(c)))
pE_z <- make_line_linear(ei_post2, y = pathFit, X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, ylab = expression(italic(F)[P]))
pG_z <- make_line_linear(ei_post2, y = hostFit, X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, ylab = expression(italic(F)[H]))
pI_z <- omega_panel_linear(ei_post2, "Path", x_lims = X_LIMS_ZOOM, x_breaks = X_BREAKS_ZOOM, x_labels = X_LABS_ZOOM,
                           ylab = expression(omega[P]))
pK_z <- omega_panel_linear(ei_post2, "Host", x_lims = X_LIMS_ZOOM, x_breaks = X_BREAKS_ZOOM, x_labels = X_LABS_ZOOM,
                           ylab = expression(omega[H]), show_xlab = TRUE)

pB_z <- make_line_linear(es_post2, y = v,       X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, show_ylab = FALSE)
pD_z <- make_line_linear(es_post2, y = s,       X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, show_ylab = FALSE)
pF_z <- make_line_linear(es_post2, y = pathFit, X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, show_ylab = FALSE)
pH_z <- make_line_linear(es_post2, y = hostFit, X_LIMS_ZOOM, X_BREAKS_ZOOM, X_LABS_ZOOM, show_ylab = FALSE)
pJ_z <- omega_panel_linear(es_post2, "Path", x_lims = X_LIMS_ZOOM, x_breaks = X_BREAKS_ZOOM, x_labels = X_LABS_ZOOM,
                           show_ylab = FALSE)
pL_z <- omega_panel_linear(es_post2, "Host", x_lims = X_LIMS_ZOOM, x_breaks = X_BREAKS_ZOOM, x_labels = X_LABS_ZOOM,
                           show_ylab = FALSE, show_xlab = TRUE)

final_plot_linear_zoom <- (
  (pA_z | pB_z) /
    (pC_z | pD_z) /
    (pE_z | pF_z) /
    (pG_z | pH_z) /
    (pI_z | pJ_z) /
    (pK_z | pL_z)
) + plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(face = "bold", size = 12))

print(final_plot_linear_zoom)

#============================== Distribution plots ====================================

# Find the range of densities across BOTH datasets
ei_max <- max(table(cut(ei_post2$v, breaks=100), cut(ei_post2$s, breaks=100)))
es_max <- max(table(cut(es_post2$v, breaks=100), cut(es_post2$s, breaks=100)))
density_max <- max(ei_max, es_max)

# Set shared limits (you can adjust these)
density_limits <- c(1, density_max)  # Start at 1 for log scale

density_limits <- c(1, density_max)  # e.g., 1..147

# nice log breaks within that range
log_brks <- scales::log_breaks(n = 5)(density_limits)

common_fill <- scale_fill_viridis_c(
  trans  = "log10",
  limits = density_limits,
  breaks = log_brks,
  labels = trans_format("log10", math_format(10^.x)),
  na.value = "white",
  oob = scales::squish,
  guide = guide_colorbar(ticks = TRUE, title.position = "top")
)

#  make the fill use the bin count (floored at 1 to avoid log(0))
a <- ggplot(ei_post2, aes(v, s)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(pmax(count, 1)))) +
  common_fill +
  labs(x = "v (virulence)", y = "c (clearance)")+
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  annotate("segment", x = 0.7, xend = 0.54, y = 0.5, yend = 0.76,
           color = "firebrick", linewidth = 1,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"))

a_inset <- ggplot(ei_post2, aes(v, s)) +
  stat_bin2d(bins = 50, aes(fill = after_stat(pmax(count, 1)))) +
  common_fill + guides(fill = "none") +
  labs(x = NULL, y = NULL)+
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.text = element_text(size = 7)
  )

b <- ggplot(es_post2, aes(v, s)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(pmax(count, 1)))) +
  common_fill +
  labs(x = "v (virulence)", y = "c (clearance)")+
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))


a_with_inset <- a +
  inset_element(a_inset, left = 0.05, bottom = 0.05, right = 0.75, top = 0.65)

# drop y stuff on the right panel
strip_y <- function(p) p + theme(axis.title.y = element_blank(),
                                 axis.text.y  = element_blank(),
                                 axis.ticks.y = element_blank())

b_plot <- strip_y(b) + labs(x = "v (virulence)")

# Compose with one merged legend
p3 <- (a_with_inset | b_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        plot.tag = element_text(face = "bold"),
        text = element_text(size = 14))

distribution_fig <- p3 + plot_annotation(tag_levels = "A")


ggsave("figures/distribution.png", distribution_fig, width = 7, height = 3.5, units = "in")
ggsave("figures/distribution.pdf", distribution_fig, width = 7, height = 3.5, units = "in")


#===================================== Dwell times =====================================

# Calculate CCDF as COUNTS not proportions
plot_data <- bind_rows(
  ei_post2 %>% filter(dwell > 0, mutator == "Path") %>% mutate(scenario = "EI Path"),
  ei_post2 %>% filter(dwell > 0, mutator == "Host") %>% mutate(scenario = "EI Host"),
  es_post2 %>% filter(dwell > 0, mutator == "Path") %>% mutate(scenario = "ES Path"),
  es_post2 %>% filter(dwell > 0, mutator == "Host") %>% mutate(scenario = "ES Host")
) %>%
  group_by(scenario) %>%
  arrange(dwell) %>%
  mutate(
    rank = row_number(),
    ccdf_count = n() - rank + 1,  # ← Actual count, not proportion
    ccdf_prop = ccdf_count / n()
  ) %>%
  ungroup()

# Plot using counts
p_es <- plot_data %>%
  filter(grepl("ES", scenario)) %>%
  ggplot(aes(x = dwell, y = ccdf_count, color = mutator)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(limits = c(0, 0.2), breaks = c(0, 0.1, 0.2)) +
  scale_color_manual(values = c("Path" = "darkblue", "Host" = "firebrick")) +
  labs(x = "Dwell time", y = "Number of events ≥ x") +
  mytheme+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.8))

p_ei <- plot_data %>%
  filter(grepl("EI", scenario)) %>%
  ggplot(aes(x = dwell, y = ccdf_count, color = mutator)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(limits = c(0, 5), breaks = c(0, 2, 4)) +
  scale_color_manual(values = c("Path" = "darkblue", "Host" = "firebrick")) +
  labs(x = "Dwell time", y = NULL) +
  mytheme+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.8))


strip_y <- function(p) {
  p + labs(y = NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

A3 <- p_ei + labs(x = "v (virulence)", y = "c (clearance)")
B3 <- strip_y(p_es) + labs(x = "v (virulence)")

p3 <- (A3 | B3) +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.position = "right",
    text = element_text(size = 14)
  )

dwell_fig <- p3 + plot_annotation(tag_levels = "A")  

ggsave("figures/dwell.png", dwell_fig, width = 7, height = 3.5, units = "in")
ggsave("figures/dwell.pdf", dwell_fig, width = 7, height = 3.5, units = "in")


#========================= Snapshots of strategies =================================
### CURRENT VS. 100 GENS AGO

# Load the time series data
fig6_data <- read.csv("postproc/fig6_timeseries.csv")

# Load data
es_full <- read.csv("results/es_simulation.csv")
es_post <- es_full %>% filter(event == "post")

# Fitness functions
beta <- 1; d0 <- 0.1; mc <- 0.1; mv <- 1.0; epsilon <- 1e-4

mortality <- function(s, v) {
  onep <- 1 + epsilon
  d0 + (mc * onep * s)/(onep - s) + (mv * onep * v)/(onep - v)
}

hostFit <- function(v, s) {
  m <- mortality(s, v)
  s / (s + m)
}

pathFit <- function(v, s) {
  m <- mortality(s, v)
  (v^beta) / (s + m)
}

# Build fitness grid
grid <- expand.grid(
  v = seq(0, 1, length.out = 100),
  s = seq(0, 1, length.out = 100)
)
grid$fH <- mapply(hostFit, grid$v, grid$s)
grid$fP <- mapply(pathFit, grid$v, grid$s)

# helper: get previous snapshot k generations back (fallback to nearest earlier)
get_prev_by_lag <- function(df, gen_now, lag_gens = 100) {
  tgt <- gen_now - lag_gens
  cand <- df %>% filter(gen <= tgt) %>% arrange(desc(gen)) %>% slice(1)
  if (nrow(cand) == 0) cand <- df %>% filter(gen < gen_now) %>%
    arrange(desc(gen)) %>% slice(1)
  cand
}

# intersection of host s = bS + mS*v  and pathogen v = bV + mV*s
line_intersection <- function(bS, mS, bV, mV) {
  den <- (1 - mV*mS)
  if (abs(den) < 1e-9) return(c(v = NA_real_, s = NA_real_))
  v <- (bV + mV*bS) / den
  s <- bS + mS*v
  c(v = v, s = s)
}

get_prev_adaptive <- function(df, gen_now,
                              step = 100,  # try 100-gens back each step
                              min_delta = 0.02,  # require ≥ 0.02 movement
                              max_steps = 20) {
  cur <- df %>% filter(gen == gen_now)
  if (nrow(cur) == 0) return(NULL)
  
  cur_x <- line_intersection(cur$bS, cur$mS, cur$bV, cur$mV)
  
  for (i in 1:max_steps) {
    tgt <- gen_now - i*step
    prev <- df %>% filter(gen <= tgt) %>% arrange(desc(gen)) %>% slice(1)
    if (nrow(prev) == 0) break
    prev_x <- line_intersection(prev$bS, prev$mS, prev$bV, prev$mV)
    if (anyNA(prev_x) || anyNA(cur_x)) next
    if (max(abs(prev_x - cur_x)) >= min_delta) return(prev)
  }
  
  # fallback: nearest earlier
  df %>% filter(gen < gen_now) %>% arrange(desc(gen)) %>% slice(1)
}

# --- helper to clip a line to [0,1]x[0,1] square (keeps only visible segment) ---
clip_unit_square <- function(df) {
  df <- df %>% filter(v >= -0.5, v <= 1.5, s >= -0.5, s <= 1.5)
  df$v <- pmin(1, pmax(0, df$v))
  df$s <- pmin(1, pmax(0, df$s))
  df
}

# --- panel function ---
make_snapshot_panel <- function(gen_num,
                                title_label = "",
                                show_x = TRUE,
                                show_y = TRUE,
                                prev_mode = c("adaptive","lag"),
                                lag_gens = 100,
                                step = 100,
                                min_delta = 0.02) {
  
  prev_mode <- match.arg(prev_mode)
  
  # current snapshot
  snapshot_current <- es_post %>% filter(gen == gen_num)
  if (nrow(snapshot_current) == 0) return(NULL)
  
  # choose previous snapshot
  prev_gen <- switch(
    prev_mode,
    adaptive = get_prev_adaptive(es_post, gen_now = gen_num,
                                 step = step, min_delta = min_delta),
    lag      = get_prev_by_lag(es_post, gen_now = gen_num, lag_gens = lag_gens)
  )
  if (is.null(prev_gen) || nrow(prev_gen) == 0) prev_gen <- snapshot_current
  
  # strategy lines (extend a bit then clip to unit square)
  v_line <- seq(-0.1, 1.1, length.out = 500)
  s_line <- seq(-0.1, 1.1, length.out = 500)
  
  host_line_current <- tibble(v = v_line,
                              s = snapshot_current$bS + snapshot_current$mS * v_line) |>
    clip_unit_square()
  path_line_current <- tibble(s = s_line,
                              v = snapshot_current$bV + snapshot_current$mV * s_line) |>
    clip_unit_square()
  host_line_prev    <- tibble(v = v_line,
                              s = prev_gen$bS + prev_gen$mS * v_line) |>
    clip_unit_square()
  path_line_prev    <- tibble(s = s_line,
                              v = prev_gen$bV + prev_gen$mV * s_line) |>
    clip_unit_square()
  
  # plot
  p <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
              fill = "white", color = "black", linewidth = 0.5) +
    geom_contour(data = grid, aes(v, s, z = fP),
                 color = "lightcoral", bins = 10, linewidth = 0.3, alpha = 0.7) +
    geom_contour(data = grid, aes(v, s, z = fH),
                 color = "steelblue",  bins = 10, linewidth = 0.3, alpha = 0.7) +
    # previous (dotted) + current (longdash)
    geom_line(data = host_line_prev, aes(v, s),
              color = "darkblue", linetype = "dotted", linewidth = 1.1, alpha = 0.85) +
    geom_line(data = path_line_prev, aes(v, s),
              color = "firebrick", linetype = "dotted", linewidth = 1.1, alpha = 0.85) +
    geom_line(data = host_line_current, aes(v, s),
              color = "darkblue", linetype = "longdash", linewidth = 1.2) +
    geom_line(data = path_line_current, aes(v, s),
              color = "firebrick", linetype = "longdash", linewidth = 1.2) +
    geom_point(aes(x = snapshot_current$v, y = snapshot_current$s),
               size = 3, color = "black") +
    coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE, clip = "on") +
    scale_x_continuous(breaks = c(0,0.5,1), labels = c( "0", "0.5", "1")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c( "0", "0.5", "1")) +
    labs(title = title_label) +
    mytheme +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white"),
      panel.border     = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.title       = element_text(hjust = 0.02, vjust = -1, size = 14),
      panel.grid       = element_blank()
    )
  
  if (!show_x) p <- p + labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  else p <- p + labs(x = "v (virulence)")
  
  if (!show_y) p <- p + labs(y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  else p <- p + labs(y = "c (clearance)")
  
  p
}

# Choose snapshots
snapshots <- c(100, 1000, 10000, 100000, 500000, 990000)
labels <- c("100 gen.", "1K gen.", "10K gen.", "100K gen.", "500K gen.", "990K gen.")

# Create panels with conditional axis labels
panels <- list()

for (i in 1:6) {
  row <- ceiling(i / 3)
  col <- ((i - 1) %% 3) + 1
  
  show_x <- (row == 2)
  show_y <- (col == 1)
  
  cat(sprintf("Creating panel %s (gen %d) - row %d, col %d\n", 
              labels[i], snapshots[i], row, col))
  
  p <- make_snapshot_panel(snapshots[i], labels[i], show_x = show_x, show_y = show_y, prev_mode = "lag", lag_gens = 100)
  if (!is.null(p)) {
    panels[[i]] <- p
  }
}

# Combine
fig6 <- (panels[[1]] | panels[[2]] | panels[[3]]) /
  (panels[[4]] + xlab(NULL) + ylab(NULL) | panels[[5]] | panels[[6]] + xlab(NULL))

print(fig6)

ggsave("figures/fig6.png", fig6, width = 8, height = 6, units = "in")
ggsave("figures/fig6.pdf", fig6, width = 8, height = 6, units = "in")
```

#=========================== Nash Regions and Occupancy ==================================

# ---------- Load Data ----------
nash_region <- read.csv("postproc/fig7_region_acute.csv")
occupancy <- read.csv("postproc/fig7_occupancy.csv")

# ES persistent states
dwell_threshold <- 0.01
es_persistent <- es_post2 %>%
  filter(event == "post", dwell > dwell_threshold) %>%
  mutate(v_bin = cut(v, breaks = 50, labels = FALSE),
         s_bin = cut(s, breaks = 50, labels = FALSE)) %>%
  group_by(v_bin, s_bin) %>%
  summarise(total_dwell = sum(dwell), .groups = "drop") %>%
  mutate(v_mid = (v_bin - 0.5) / 50,
         s_mid = (s_bin - 0.5) / 50)

# ---------- Fitness Grid for Contours ----------
nS <- 0.1; nV <- 1.0; d0 <- 0.1; beta <- 1.0
ONE_PLUS_EPS <- 1 + 1/999

mortality <- function(v, s) {
  s_term <- (nS * ONE_PLUS_EPS * s) / (ONE_PLUS_EPS - s)
  v_term <- (nV * ONE_PLUS_EPS * v) / (ONE_PLUS_EPS - v)
  d0 + s_term + v_term
}

host_fitness <- function(v, s) {
  m <- mortality(v, s)
  s / (s + m)
}

path_fitness <- function(v, s) {
  m <- mortality(v, s)
  (v^beta) / (s + m)
}

fitness_grid <- expand.grid(
  v = seq(0.01, 0.99, length.out = 150),
  s = seq(0.01, 0.99, length.out = 150)
) %>%
  mutate(
    f_H = mapply(host_fitness, v, s),
    f_P = mapply(path_fitness, v, s)
  )


# Finite difference step
h <- 1e-4

# Gradient helper functions
grad_host_v <- function(v, s) {
  v1 <- max(0.01, v - h); v2 <- min(0.99, v + h)
  (host_fitness(v2, s) - host_fitness(v1, s)) / (v2 - v1)
}

grad_host_s <- function(v, s) {
  s1 <- max(0.01, s - h); s2 <- min(0.99, s + h)
  (host_fitness(v, s2) - host_fitness(v, s1)) / (s2 - s1)
}

grad_path_v <- function(v, s) {
  v1 <- max(0.01, v - h); v2 <- min(0.99, v + h)
  (path_fitness(v2, s) - path_fitness(v1, s)) / (v2 - v1)
}

grad_path_s <- function(v, s) {
  s1 <- max(0.01, s - h); s2 <- min(0.99, s + h)
  (path_fitness(v, s2) - path_fitness(v, s1)) / (s2 - s1)
}

# ---------- Calculate Violation Grid ----------

# Create grid (smaller for speed)
violation_grid <- expand.grid(
  v = seq(0.02, 0.98, length.out = 80),  # Smaller grid = faster
  s = seq(0.02, 0.98, length.out = 80)
)

# Calculate gradients and slopes
violation_grid <- violation_grid %>%
  rowwise() %>%
  mutate(
    # Host gradient
    FH_v = grad_host_v(v, s),
    FH_s = grad_host_s(v, s),
    
    # Pathogen gradient
    FP_v = grad_path_v(v, s),
    FP_s = grad_path_s(v, s),
    
    # Tangent slopes (perpendicular to gradients)
    # Strategy slope = -gradient_perpendicular / gradient_parallel
    m_c = -FP_v / (FP_s + 1e-12),  # Host strategy slope
    m_v = -FH_s / (FH_v + 1e-12),  # Pathogen strategy slope
    
    # Stability product
    prod_mv = m_c * m_v,
    
    # Classify violation
    zone = case_when(
      prod_mv >= 1  ~ "violates (≥1)",   # RED
      prod_mv <= -1 ~ "violates (≤-1)",  # BLUE
      TRUE          ~ "compatible"        # WHITE
    )
  ) %>%
  ungroup()

# ---------- Panel A: Violation Map ----------

pA_violations <- ggplot() +
  # Violation zones
  geom_tile(data = violation_grid, 
            aes(x = v, y = s, fill = zone)) +
  scale_fill_manual(
    values = c("violates (≥1)" = "lightpink",    # Red
               "compatible" = "gray80",          # White
               "violates (≤-1)" = "lightblue"),   # Blue
    name = "Stability"
  ) +
  geom_contour(data = fitness_grid, 
               aes(x = v, y = s, z = f_H), 
               color = "steelblue", alpha = 0.5, bins = 10, linewidth = 0.3) +
  geom_contour(data = fitness_grid, 
               aes(x = v, y = s, z = f_P), 
               color = "lightcoral", alpha = 0.5, bins = 10, linewidth = 0.3) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "virulence (v)", y = "clearance (c)") +
  mytheme+
  theme(panel.grid = element_blank(), axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))

# Create a binary field for contour drawing
contour_data <- violation_grid %>%
  mutate(compatible_binary = ifelse(zone == "compatible", 1, 0))

# ---------- Panel B with Smooth Red Contour ----------

pB <- ggplot() +
  # Occupancy heatmap
  geom_tile(data = occupancy %>% filter(occupancy > 0), 
            aes(x = v_bin_mid, y = s_bin_mid, fill = occupancy)) +
  scale_fill_viridis_c(
    trans = "log10",
    limits = c(1e-4, 1e2),  # Adjusted limits
    breaks = c(1e-4, 1e-2, 1e0, 1e2),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    name = "Time",
    option = "magma"
  )  +
  # Smooth red contour around compatible region
  geom_contour(data = contour_data,
               aes(x = v, y = s, z = compatible_binary),
               breaks = 0.5,  # Boundary between 0 and 1
               color = "gold", linewidth = 1) +
  # Fitness contours
  geom_contour(data = fitness_grid, 
               aes(x = v, y = s, z = f_H), 
               color = "steelblue", alpha = 0.3, bins = 8, linewidth = 0.3) +
  geom_contour(data = fitness_grid, 
               aes(x = v, y = s, z = f_P), 
               color = "light coral", alpha = 0.3, bins = 8, linewidth = 0.3) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "virulence (v)", y = NULL) +
  mytheme +
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))

# ---------- Panel C with Smooth Red Contour ----------

pC <- ggplot() +
  # Persistent occupancy heatmap
  geom_tile(data = es_persistent %>% filter(total_dwell > 0), 
            aes(x = v_mid, y = s_mid, fill = total_dwell)) +
  scale_fill_viridis_c(
    trans = "log10",
    limits = c(1e-4, 1e2),  # Adjusted limits
    breaks = c(1e-4, 1e-2, 1e0, 1e2),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    name = "Time",
    option = "magma"
  ) +
  # Smooth red contour
  geom_contour(data = contour_data,
               aes(x = v, y = s, z = compatible_binary),
               breaks = 0.5,
               color = "gold", linewidth = 1) +
  # Fitness contours
  geom_contour(data = fitness_grid, 
               aes(x = v, y = s, z = f_H), 
               color = "steelblue", alpha = 0.5, bins = 10, linewidth = 0.3) +
  geom_contour(data = fitness_grid, 
               aes(x = v, y = s, z = f_P), 
               color = "lightcoral", alpha = 0.5, bins = 10, linewidth = 0.3) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "virulence (v)", y = "clearance (c)") +
  mytheme +
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())+
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))

# Combine all three
fig7_full <- pA_violations + pB + pC + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.position = "right",
    text = element_text(size = 14)
  )

print(fig7_full)

ggsave("figures/fig7.png", fig7_full, width = 8, height = 4, units = "in")
ggsave("figures/fig7.pdf", fig6, width = 8, height = 4, units = "in")

