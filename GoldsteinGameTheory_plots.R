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
library(zoo)


# Working directory
setwd("~/Documents/GitHub/GoldsteinGameTheory")

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

set_theme(mytheme)
#===================================== FIGURE 1 =============================================
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

#-----------------------------------------
# FOR MINIMAL FITNESS PLOTS USE THESE!!! 
#-----------------------------------------

#f_H_minimal <- function(c, v) {
#  c * (1 - c) * (1 - v)
#}

#f_P_minimal <- function(c, v) {
#  v * (1 - v) * (1 - c)
#}

#-----------------------------------------

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
  geom_line(data = opt_c_ac, aes(x=v, y=c), color="steelblue", size=1.5, inherit.aes=FALSE, linetype = "dashed") +
  labs(x="v (virulence)", y="c (clearance)") +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  scale_fill_viridis_d(option = "viridis", 
                       labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", 
                                  "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))

B <- ggplot(grid, aes(x=v, y=c, z=f_P_acute)) +
  geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
  geom_line(data = opt_v_ac, aes(x=v, y=c), color="lightcoral", size=1.5, inherit.aes=FALSE, linetype = "dashed") +
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
  geom_vline(xintercept = 0.5278595, size = 1.5, color = "firebrick", linetype = "solid") +
  geom_hline(yintercept = 0.7773371, size = 1.5, color = "darkblue", linetype = "solid") +
  geom_line(data = opt_c_ac, aes(x = v, y = c), color = "steelblue", inherit.aes = FALSE, size = 1.5, linetype = "dashed") +
  geom_line(data = opt_v_ac, aes(x = v, y = c), color = "lightcoral", inherit.aes = FALSE, size = 1.5, linetype = "dashed") +
  geom_point(data = opt_joint_acute, aes(x = v, y = c), color = "grey20", size = 5, inherit.aes = FALSE) +
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

ggsave("figures/Figure1.png", final_fig, width = 10, height = 4, units = "in")
ggsave("figures/Figure1.pdf", final_fig, width = 10, height = 4, units = "in")


#===============================================================================
# MINIMAL FITNESS
#===============================================================================

# ==========================================
# MINIMAL FITNESS MODEL - THREE PANEL PLOT
# ==========================================

# Define minimal fitness functions
f_H_minimal <- function(c, v) {
  c * (1 - c) * (1 - v)
}

f_P_minimal <- function(c, v) {
  v * (1 - v) * (1 - c)
}

# Create grid
c_vals <- seq(0, 1, length.out = 300)
v_vals <- seq(0, 1, length.out = 300)
grid_minimal <- expand.grid(c = c_vals, v = v_vals)

# Evaluate fitness functions
grid_minimal <- grid_minimal %>%
  mutate(
    f_H_minimal = mapply(f_H_minimal, c, v),
    f_P_minimal = mapply(f_P_minimal, c, v),
    joint_fitness_minimal = f_H_minimal * f_P_minimal
  )

# Find max fitness for normalization
fmax_H_minimal <- max(grid_minimal$f_H_minimal)
fmax_P_minimal <- max(grid_minimal$f_P_minimal)

# Normalize
grid_minimal <- grid_minimal %>%
  mutate(
    f_H_minimal = f_H_minimal / fmax_H_minimal,
    f_P_minimal = f_P_minimal / fmax_P_minimal,
    joint_fitness_minimal = f_H_minimal * f_P_minimal
  )

# Calculate best response curves
v_vals_br <- seq(0.01, 0.99, length.out = 300)
c_vals_br <- seq(0.01, 0.99, length.out = 300)

# Host best response to pathogen virulence
opt_c_minimal <- data.frame(
  v = v_vals_br,
  c = sapply(v_vals_br, function(v) {
    optimize(function(c) -f_H_minimal(c, v), c(1e-3, 0.99))$minimum
  })
)

# Pathogen best response to host clearance
opt_v_minimal <- data.frame(
  c = c_vals_br,
  v = sapply(c_vals_br, function(c) {
    v_seq <- seq(0.01, 0.99, length.out = 200)
    f_vals <- sapply(v_seq, function(v) f_P_minimal(c, v))
    v_seq[which.max(f_vals)]
  })
)

# Find Nash equilibrium
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

opt_joint_minimal <- find_nash_brute_force(opt_c_minimal, opt_v_minimal)

cat("Nash equilibrium for minimal fitness:\n")
cat("v =", opt_joint_minimal$v, ", c =", opt_joint_minimal$c, "\n")

# ==========================================
# CREATE THREE PANELS
# ==========================================

# Panel A: Host fitness landscape
A_minimal <- ggplot(grid_minimal, aes(x = v, y = c, z = f_H_minimal)) +
  geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
  geom_line(data = opt_c_minimal, aes(x = v, y = c), 
            color = "steelblue", linewidth = 1.5, inherit.aes = FALSE, linetype = "dashed") +
  labs(x = "v (virulence)", y = "c (clearance)") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_fill_viridis_d(option = "viridis", 
                       name = "Host\nFitness",
                       labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", 
                                  "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1")) +
  mytheme

# Panel B: Pathogen fitness landscape
B_minimal <- ggplot(grid_minimal, aes(x = v, y = c, z = f_P_minimal)) +
  geom_contour_filled(breaks = seq(0, 1, length.out = 10)) +
  geom_line(data = opt_v_minimal, aes(x = v, y = c), 
            color = "lightcoral", linewidth = 1.5, inherit.aes = FALSE, linetype = "dashed") +
  labs(x = "v (virulence)", y = "c (clearance)") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_fill_viridis_d(option = "viridis", 
                       name = "Pathogen\nFitness",
                       labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", 
                                  "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1")) +
  mytheme

# Panel C: Joint fitness with best responses and Nash
C_minimal <- ggplot(grid_minimal, aes(x = v, y = c, z = joint_fitness_minimal)) +
  geom_contour(aes(z = f_H_minimal), bins = 10, color = "steelblue", alpha = 0.5) +
  geom_contour(aes(z = f_P_minimal), bins = 10, color = "lightcoral", alpha = 0.5) +
  geom_vline(xintercept = opt_joint_minimal$v, linewidth = 1.5, 
             color = "firebrick", linetype = "solid") +
  geom_hline(yintercept = opt_joint_minimal$c, linewidth = 1.5, 
             color = "darkblue", linetype = "solid") +
  geom_line(data = opt_c_minimal, aes(x = v, y = c), 
            color = "steelblue", inherit.aes = FALSE, linewidth = 1.5, linetype = "dashed") +
  geom_line(data = opt_v_minimal, aes(x = v, y = c), 
            color = "lightcoral", inherit.aes = FALSE, linewidth = 1.5, linetype = "dashed") +
  geom_point(data = opt_joint_minimal, aes(x = v, y = c), 
             color = "grey20", size = 5, inherit.aes = FALSE) +
  labs(x = "v (virulence)", y = "c (clearance)") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  mytheme

# ==========================================
# ARRANGE AS FACETS (same style as original)
# ==========================================

strip_y <- function(p) {
  p + labs(y = NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

A1_minimal <- A_minimal + labs(x = NULL, y = "c (clearance)")
B1_minimal <- strip_y(B_minimal) + labs(x = "v (virulence)")
C1_minimal <- strip_y(C_minimal) + labs(x = NULL)

p_minimal <- (A1_minimal | B1_minimal | C1_minimal) +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.position = "right",
    text = element_text(size = 14)
  )

final_fig_minimal <- p_minimal + plot_annotation(
  tag_levels = "A",
)

print(final_fig_minimal)
ggsave("fitness_landscapes_minimal.pdf", final_fig_minimal, width = 14, height = 5)


#================================ FIGURE 2 =============================================
#============= Multiple Nash equlibria and example for tree intersections ==============


# Core parameters and functions
# =============================

beta <- 1
d0 <- 0.1
mc <- 0.1
mv <- 1.0
epsilon <- 1e-4

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

# Gradient calculation for Nash equilibria
grad <- function(f, v, s, h = 1e-5) {
  v1 <- clamp01(v - h); v2 <- clamp01(v + h)
  s1 <- clamp01(s - h); s2 <- clamp01(s + h)
  dv <- (f(v2, s) - f(v1, s)) / max(v2 - v1, 1e-12)
  ds <- (f(v, s2) - f(v, s1)) / max(s2 - s1, 1e-12)
  list(dv = dv, ds = ds)
}

# Strategy line helpers
host_line <- function(v, bS, mS) clamp01(bS + mS * v)
path_line <- function(s, bV, mV) clamp01(bV + mV * s)

# Build fitness landscape grid
# =============================

v_vals <- seq(0, 1, length.out = 300)
s_vals <- seq(0, 1, length.out = 300)
grid <- expand.grid(v = v_vals, s = s_vals)
grid$fH <- mapply(hostFit, grid$v, grid$s)
grid$fP <- mapply(pathFit, grid$v, grid$s)

# High resolution for smooth strategy lines
v_plot <- seq(0, 1, length.out = 900)
s_plot <- seq(0, 1, length.out = 900)

# Find ET (Evolved Trait) equilibrium - strong Nash
# =================================================

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

# Find intersection of best response curves
distances <- expand.grid(i = 1:nrow(BR_host), j = 1:nrow(BR_path))
distances$dist <- sqrt((BR_host$v[distances$i] - BR_path$v[distances$j])^2 + 
                         (BR_host$s[distances$i] - BR_path$s[distances$j])^2)
closest <- distances[which.min(distances$dist), ]
v_star <- BR_host$v[closest$i]
s_star <- BR_host$s[closest$i]

# Colors
col_host <- "darkblue"  # Blue
col_path <- "firebrick"  # Red/orange

# Panel plotting function
# ========================

plot_panel <- function(bV, mV, bS, mS, v_int, s_int, 
                       show_stable = TRUE, 
                       boundary_points = NULL,
                       show_yaxis = TRUE,
                       show_best_response = TRUE) {
  
  # Create strategy line data
  host_data <- data.frame(v = v_plot, s = host_line(v_plot, bS, mS))
  path_data <- data.frame(v = path_line(s_plot, bV, mV), s = s_plot)
  
  # Base plot with fitness contours
  p <- ggplot(grid, aes(v, s)) +
    geom_contour(aes(z = fP), color = "lightcoral", bins = 10, 
                 linewidth = 0.3, alpha = 0.7) +
    geom_contour(aes(z = fH), color = "steelblue", bins = 10, 
                 linewidth = 0.3, alpha = 0.7)
  
  # Add best response curves (optimal strategies)
  if (show_best_response) {
    p <- p +
      geom_line(data = BR_host, aes(v, s), 
                color = "lightcoral", linewidth = 2, linetype = "dashed") +
      geom_line(data = BR_path, aes(v, s), 
                color = "steelblue", linewidth = 2, linetype = "dashed")
  }
  
  # Add strategy lines (dashed)
  p <- p +
    geom_line(data = host_data, aes(v, s), 
              linetype = "solid", linewidth = 1.5, color = col_host) +
    geom_line(data = path_data, aes(v, s), 
              linetype = "solid", linewidth = 1.5, color = col_path) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    labs(x = "v (virulence)") +
    mytheme
  
  # Y-axis handling
  if (show_yaxis) {
    p <- p + labs(y = "c (clearance)")
  } else {
    p <- p + labs(y = NULL) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  # Add equilibrium points
  if (show_stable) {
    # Stable equilibrium - filled black circle
    p <- p + geom_point(aes(x = v_int, y = s_int), 
                        size = 5, colour = "black")
  } else {
    # Unstable equilibrium - gray filled circle
    p <- p + geom_point(aes(x = v_int, y = s_int), 
                        size = 5, shape = 21, 
                        fill = "gray70", colour = "black", stroke = 1)
  }
  
  # Add boundary equilibria if provided
  if (!is.null(boundary_points)) {
    p <- p + geom_point(data = boundary_points, aes(x = v, y = s),
                        size = 5, colour = "black")
  }
  
  return(p)
}

# PANEL A: ET equilibrium (mS = mV = 0)
# =====================================

cat("\nGenerating Panel A: ET equilibrium...\n")
bS_A <- s_star  # Fixed at s*
mS_A <- 0.0
bV_A <- v_star  # Fixed at v*
mV_A <- 0.0

pA <- plot_panel(bV_A, mV_A, bS_A, mS_A, v_star, s_star,
                 show_stable = TRUE, show_yaxis = TRUE)


# PANEL B: Host slope varies, point conserved
# ============================================

# Host develops positive responsiveness
mS_B <- 0.5
bS_B <- s_star - mS_B * v_star  # Adjust intercept to maintain point
bV_B <- v_star  # Pathogen still fixed
mV_B <- 0.0

pB <- plot_panel(bV_B, mV_B, bS_B, mS_B, v_star, s_star,
                 show_stable = TRUE, show_yaxis = FALSE)


# PANEL C: Pathogen slope varies, point conserved
# ===============================================

bS_C <- s_star  # Host still fixed
mS_C <- 0.0
mV_C <- 0.8  # Pathogen develops positive responsiveness
bV_C <- v_star - mV_C * s_star  # Adjust intercept to maintain point

pC <- plot_panel(bV_C, mV_C, bS_C, mS_C, v_star, s_star,
                 show_stable = TRUE, show_yaxis = FALSE)

# PANEL D: Large slopes, destabilization (|mS * mV| > 1)
# ======================================================

# Both slopes large, creating instability
mS_D <- 2.5
mV_D <- 2.5
bS_D <- s_star - mS_D * v_star
bV_D <- v_star - mV_D * s_star

# Adjust intercepts to force boundary behavior
margin <- 1e-3
if (bS_D > -margin) bS_D <- -margin
if (bS_D + mS_D < 1 + margin) bS_D <- (1 + margin) - mS_D
if (bV_D > -margin) bV_D <- -margin
if (bV_D + mV_D < 1 + margin) bV_D <- (1 + margin) - mV_D

# Find interior intersection (unstable)
v_int_D <- v_star
s_int_D <- s_star

# Boundary equilibria (stable attractors)
boundary_D <- data.frame(
  v = c(clamp01(bV_D), clamp01(bV_D + mV_D)),
  s = c(clamp01(bS_D), clamp01(bS_D + mS_D))
)

pD <- plot_panel(bV_D, mV_D, bS_D, mS_D, v_int_D, s_int_D,
                 show_stable = FALSE, 
                 boundary_points = boundary_D,
                 show_yaxis = FALSE)

# Combine panels and save
# =======================

final_fig_2 <- (pA | pB | pC | pD) +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(widths = c(1, 1, 1, 1)) &  # Note the & instead of +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(-0.01, 0.76)  # (x, y) on 0-1 scale
  )

ggsave("figures/Figure2.png", final_fig_2, width = 9, height = 6, units = "in")
ggsave("figures/Figure2.pdf", final_fig_2, width = 9, height = 6, units = "in")

# Print stability conditions
cat("\n" %+% strrep("=", 70) %+% "\n")
cat("Stability conditions:\n")
cat(sprintf("Panel A: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (stable)\n", 
            mS_A, mV_A, abs(mS_A * mV_A)))
cat(sprintf("Panel B: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (stable)\n", 
            mS_B, mV_B, abs(mS_B * mV_B)))
cat(sprintf("Panel C: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (stable)\n", 
            mS_C, mV_C, abs(mS_C * mV_C)))
cat(sprintf("Panel D: mS=%.2f, mV=%.2f, |mS*mV|=%.2f (UNSTABLE, >1)\n", 
            mS_D, mV_D, abs(mS_D * mV_D)))
cat(strrep("=", 70) %+% "\n")

# Create standalone legend
# =========================

legend_data <- data.frame(
  x = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
  y = c(0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25),
  label = c(
    "Host response rule: c = bc + mc v",
    "Pathogen response rule: v = bv + mv c",
    "Host best response (optimal c given v)",
    "Pathogen best response (optimal v given c)",
    "Host fitness contours",
    "Pathogen fitness contours",
    "Stable equilibrium",
    "Unstable equilibrium"
  )
)

p_legend <- ggplot(legend_data, aes(x, y)) +
  # Strategy lines (dashed)
  annotate("segment", x = 0.05, xend = 0.15, y = 0.95, yend = 0.95,
           linetype = "dashed", linewidth = 1.2, color = col_host) +
  annotate("segment", x = 0.05, xend = 0.15, y = 0.85, yend = 0.85,
           linetype = "dashed", linewidth = 1.2, color = col_path) +
  # Best response curves (solid)
  annotate("segment", x = 0.05, xend = 0.15, y = 0.75, yend = 0.75,
           linetype = "solid", linewidth = 0.8, color = col_host, alpha = 0.9) +
  annotate("segment", x = 0.05, xend = 0.15, y = 0.65, yend = 0.65,
           linetype = "solid", linewidth = 0.8, color = col_path, alpha = 0.9) +
  # Contours
  annotate("segment", x = 0.05, xend = 0.15, y = 0.55, yend = 0.55,
           linewidth = 0.5, color = col_host, alpha = 0.7) +
  annotate("segment", x = 0.05, xend = 0.15, y = 0.45, yend = 0.45,
           linewidth = 0.5, color = col_path, alpha = 0.7) +
  # Points
  annotate("point", x = 0.1, y = 0.35, size = 3, color = "black") +
  annotate("point", x = 0.1, y = 0.25, size = 3, shape = 21, 
           fill = "gray70", color = "black", stroke = 1) +
  # Labels
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.03, size = 3.5) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

#=============================== FIGURE 3 ======================================================
#============================Simulation results ================================================
## Trait Evolution, fitness over time, relative substitution rates

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
  
  # Read omega directly from CSV, don't calculate from dwell
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

# ---------- Line panels for v, c, W ----------
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

# ---------- Build Figure 3 with omega panels ----------
# Original simulations
# Acute fitness functions, all simulations
ei_acute                <- read.csv("organized_test_simulations/Acute_Fitness_Model/20251028_Host&Path_Reactive_m&a_0.01/results/ei_simulation.csv")%>% filter(event == "post", gen > 10000)
es_acute                <- read.csv("organized_test_simulations/Acute_Fitness_Model/20251028_Host&Path_Reactive_m&a_0.01/results/es_simulation.csv")%>% filter(event == "post", gen > 10000)
es_acute_fixed_host     <- read.csv("organized_test_simulations/Acute_Fitness_Model/20251104_Fixed_Host_m&a_0.01/results/es_simulation.csv")%>% filter(event == "post", gen > 10000)
es_acute_fixed_pathogen <- read.csv("organized_test_simulations/Acute_Fitness_Model/20251104_Fixed_Path_m&a_0.01/results/es_simulation.csv")%>% filter(event == "post", gen > 10000)

# Minimal fitness, all simulations 
ei_acute_both_reactive_minimal_fit  <- read.csv("organized_test_simulations/Minimal_Fitness_Model/20251030_Host&Path_Reactive_m&a_0.01/results/ei_simulation.csv")%>% filter(event == "post", gen > 10000)
es_acute_both_reactive_minimal_fit  <- read.csv("organized_test_simulations/Minimal_Fitness_Model/20251030_Host&Path_Reactive_m&a_0.01/results/es_simulation.csv")%>% filter(event == "post", gen > 10000)
es_acute_fixed_host_minimal_fit     <- read.csv("organized_test_simulations/Minimal_Fitness_Model/20251030_Fixed_Host_m&a_0.01/results/es_simulation.csv")%>% filter(event == "post", gen > 10000)
es_acute_fixed_pathogen_minimal_fit <- read.csv("organized_test_simulations/Minimal_Fitness_Model/20251030_Fixed_Path_m&a_0.01/results/es_simulation.csv")%>% filter(event == "post", gen > 10000)

ei_acute_thin                <- thin_for_plot(ei_acute, every = 100)
es_acute_thin                <- thin_for_plot(es_acute, every = 100)
es_acute_fixed_host_thin     <- thin_for_plot(es_acute_fixed_host, every = 100)
es_acute_fixed_pathogen_thin <- thin_for_plot(es_acute_fixed_pathogen, every = 100)

ei_acute_both_reactive_minimal_fit_thin  <- thin_for_plot(ei_acute_both_reactive_minimal_fit, every = 100)
es_acute_both_reactive_minimal_fit_thin  <- thin_for_plot(es_acute_both_reactive_minimal_fit, every = 100)
es_acute_fixed_host_minimal_fit_thin     <- thin_for_plot(es_acute_fixed_host_minimal_fit, every = 100)
es_acute_fixed_pathogen_minimal_fit_thin <- thin_for_plot(es_acute_fixed_pathogen_minimal_fit, every = 100)

# ---------- LINEAR X settings ----------
library(scales)

X_LIMS_LIN   <- c(1e4, 1e6)
X_BREAKS_LIN <- c(10000, 505000, 1e6)
X_LABS_LIN   <- c("10K", "505K", "1M")

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
pA_lin <- make_line_linear(ei_acute_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(v)))
pE_lin <- make_line_linear(ei_acute_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(c)))
pI_lin <- make_line_linear(ei_acute_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(W[P]))
pM_lin <- make_line_linear(ei_acute_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(W[H]))
pQ_lin <- omega_panel_linear(ei_acute_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             ylab = expression(omega[P]))
pU_lin <- omega_panel_linear(ei_acute_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             ylab = expression(omega[H]), show_xlab = TRUE)

pB_lin <- make_line_linear(es_acute_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pF_lin <- make_line_linear(es_acute_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pJ_lin <- make_line_linear(es_acute_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pN_lin <- make_line_linear(es_acute_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pR_lin <- omega_panel_linear(es_acute_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE)
pV_lin <- omega_panel_linear(es_acute_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE, show_xlab = TRUE)

pC_lin <- make_line_linear(es_acute_fixed_host_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pG_lin <- make_line_linear(es_acute_fixed_host_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pK_lin <- make_line_linear(es_acute_fixed_host_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pO_lin <- make_line_linear(es_acute_fixed_host_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pS_lin <- omega_panel_linear(es_acute_fixed_host_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE)
pW_lin <- omega_panel_linear(es_acute_fixed_host_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE, show_xlab = TRUE)

pD_lin <- make_line_linear(es_acute_fixed_pathogen_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pH_lin <- make_line_linear(es_acute_fixed_pathogen_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pL_lin <- make_line_linear(es_acute_fixed_pathogen_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pP_lin <- make_line_linear(es_acute_fixed_pathogen_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pT_lin <- omega_panel_linear(es_acute_fixed_pathogen_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE)
pX_lin <- omega_panel_linear(es_acute_fixed_host_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE, show_xlab = TRUE)

total_thin <- (
  (pA_lin | pC_lin | pD_lin| pB_lin ) /
    (pE_lin  | pG_lin | pH_lin| pF_lin) /
    (pI_lin  | pK_lin | pL_lin| pJ_lin) /
    (pM_lin  | pO_lin | pP_lin| pN_lin) /
    (pQ_lin  | pS_lin | pT_lin| pR_lin) /
    (pU_lin  | pW_lin | pX_lin| pV_lin)
) + plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(face = "bold", size = 12))

ggsave("figures/Figure3_lin.png", total_thin, width = 9, height = 9, units = "in")
ggsave("figures/Figure3_lin.pdf", total_thin, width = 9, height = 9, units = "in")


# ---------- Assemble: linear FULL (0–1e6) - MINIMAL FITNESS ----------
pA_lin_m <- make_line_linear(ei_acute_both_reactive_minimal_fit_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(v)))
pE_lin_m <- make_line_linear(ei_acute_both_reactive_minimal_fit_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(italic(c)))
pI_lin_m <- make_line_linear(ei_acute_both_reactive_minimal_fit_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(W[P]))
pM_lin_m <- make_line_linear(ei_acute_both_reactive_minimal_fit_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, ylab = expression(W[H]))
pQ_lin_m <- omega_panel_linear(ei_acute_both_reactive_minimal_fit_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             ylab = expression(omega[P]))
pU_lin_m <- omega_panel_linear(ei_acute_both_reactive_minimal_fit_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             ylab = expression(omega[H]), show_xlab = TRUE)

pB_lin_m <- make_line_linear(es_acute_both_reactive_minimal_fit_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pF_lin_m <- make_line_linear(es_acute_both_reactive_minimal_fit_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pJ_lin_m <- make_line_linear(es_acute_both_reactive_minimal_fit_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pN_lin_m <- make_line_linear(es_acute_both_reactive_minimal_fit_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pR_lin_m <- omega_panel_linear(es_acute_both_reactive_minimal_fit_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE)
pV_lin_m <- omega_panel_linear(es_acute_both_reactive_minimal_fit_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE, show_xlab = TRUE)

pC_lin_m <- make_line_linear(es_acute_fixed_host_minimal_fit_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pG_lin_m <- make_line_linear(es_acute_fixed_host_minimal_fit_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pK_lin_m <- make_line_linear(es_acute_fixed_host_minimal_fit_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pO_lin_m <- make_line_linear(es_acute_fixed_host_minimal_fit_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pS_lin_m <- omega_panel_linear(es_acute_fixed_host_minimal_fit_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE)
pW_lin_m <- omega_panel_linear(es_acute_fixed_host_minimal_fit_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE, show_xlab = TRUE)

pD_lin_m <- make_line_linear(es_acute_fixed_pathogen_minimal_fit_thin, y = v,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pH_lin_m <- make_line_linear(es_acute_fixed_pathogen_minimal_fit_thin, y = s,       X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pL_lin_m <- make_line_linear(es_acute_fixed_pathogen_minimal_fit_thin, y = pathFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pP_lin_m <- make_line_linear(es_acute_fixed_pathogen_minimal_fit_thin, y = hostFit, X_LIMS_LIN, X_BREAKS_LIN, X_LABS_LIN, show_ylab = FALSE)
pT_lin_m <- omega_panel_linear(es_acute_fixed_pathogen_minimal_fit_thin, "Path", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                             show_ylab = FALSE)
pX_lin_m   <- omega_panel_linear(es_acute_fixed_pathogen_minimal_fit_thin, "Host", x_lims = X_LIMS_LIN, x_breaks = X_BREAKS_LIN, x_labels = X_LABS_LIN,
                              show_ylab = FALSE, show_xlab = TRUE)

total_thin_minimal <- (
  (pA_lin_m   | pC_lin_m | pD_lin_m| pB_lin_m ) /
    (pE_lin_m | pG_lin_m | pH_lin_m| pF_lin_m) /
    (pI_lin_m | pK_lin_m | pL_lin_m| pJ_lin_m ) /
    (pM_lin_m | pO_lin_m | pP_lin_m| pN_lin_m ) /
    (pQ_lin_m | pS_lin_m | pT_lin_m| pR_lin_m) /
    (pU_lin_m | pW_lin_m | pX_lin_m| pV_lin_m)
) + plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(face = "bold", size = 12))

ggsave("figures/Figure3_lin_minimal.png", total_thin_minimal, width = 9, height = 9, units = "in")
ggsave("figures/Figure3_lin_minimal.pdf", total_thin_minimal, width = 9, height = 9, units = "in")

#============================== Distribution plots ====================================

# Find the range of densities across BOTH datasets
ei_max <- max(table(cut(ei_acute$v, breaks=100), cut(ei_acute$s, breaks=100)))
es_max <- max(table(cut(es_acute$v, breaks=100), cut(es_acute$s, breaks=100)))
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
a <- ggplot(ei_acute, aes(v, s)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(pmax(count, 1)))) +
  common_fill +
  labs(x = "v (virulence)", y = "c (clearance)")+
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  annotate("segment", x = 0.7, xend = 0.54, y = 0.5, yend = 0.76,
           color = "firebrick", linewidth = 1,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"))

a_inset <- ggplot(ei_acute, aes(v, s)) +
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

b <- ggplot(es_acute, aes(v, s)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(pmax(count, 1)))) +
  common_fill +
  labs(x = "v (virulence)", y = "c (clearance)")+
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))


c <- ggplot(es_acute_fixed_host, aes(v, s)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(pmax(count, 1)))) +
  common_fill +
  labs(x = "v (virulence)", y = "c (clearance)")+
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))


d <- ggplot(es_acute_fixed_pathogen, aes(v, s)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(pmax(count, 1)))) +
  common_fill +
  labs(x = "v (virulence)", y = "c (clearance)")+
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))

# drop y stuff on the right panel
strip_y <- function(p) p + theme(axis.title.y = element_blank(),
                                 axis.text.y  = element_blank(),
                                 axis.ticks.y = element_blank())
strip_x <- function(p) p + theme(axis.title.x = element_blank(),
                                 axis.text.x  = element_blank(),
                                 axis.ticks.x = element_blank())

a_plot <- strip_x(a) 

a_with_inset <- a_plot +
  inset_element(a_inset, left = 0.05, bottom = 0.05, right = 0.75, top = 0.65)

b_plot <- strip_y(b) 
b_plot <- strip_x(b_plot)
c_plot <- c
d_plot <- strip_y(d) + labs(x = "v (virulence)")

# Compose with one merged legend
p3 <- ((a_with_inset | b_plot) / (c_plot | d_plot)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        plot.tag = element_text(face = "bold"),
        text = element_text(size = 14))

distribution_fig <- p3 + plot_annotation(tag_levels = "A")

ggsave("figures/distribution.png", distribution_fig, width = 7.5, height = 6, units = "in")
ggsave("figures/distribution.pdf", distribution_fig, width = 7.5, height = 6, units = "in")

#===================================== Dwell times =====================================

# Calculate CCDF as COUNTS not proportions
plot_data <- bind_rows(
  ei_acute %>% filter(dwell > 0, mutator == "Path") %>% mutate(scenario = "EI Path"),
  ei_acute %>% filter(dwell > 0, mutator == "Host") %>% mutate(scenario = "EI Host"),
  es_acute %>% filter(dwell > 0, mutator == "Path") %>% mutate(scenario = "ES Path"),
  es_acute %>% filter(dwell > 0, mutator == "Host") %>% mutate(scenario = "ES Host")
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
  scale_x_continuous(limits = c(0, 0.1), breaks = c(0, 0.1)) +
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
  scale_x_continuous(limits = c(0, 3), breaks = c(0, 3)) +
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

pdwell <- (A3 | B3) +
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
fig6_data <- read.csv("~/Documents/GitHub/GoldsteinGameTheory/organized_test_simulations/Acute_Fitness_Model/20251028_Host&Path_Reactive_m&a_0.01/postproc/fig6_timeseries.csv")

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
  snapshot_current <- es_acute %>% filter(gen == gen_num)
  if (nrow(snapshot_current) == 0) return(NULL)
  
  # choose previous snapshot
  prev_gen <- switch(
    prev_mode,
    adaptive = get_prev_adaptive(es_acute, gen_now = gen_num,
                                 step = step, min_delta = min_delta),
    lag      = get_prev_by_lag(es_acute, gen_now = gen_num, lag_gens = lag_gens)
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
              color = "darkblue", linetype = "dotted", linewidth = 1.5, alpha = 0.85) +
    geom_line(data = path_line_prev, aes(v, s),
              color = "firebrick", linetype = "dotted", linewidth = 1.5, alpha = 0.85) +
    geom_line(data = host_line_current, aes(v, s),
              color = "darkblue", linetype = "solid", linewidth = 1.5) +
    geom_line(data = path_line_current, aes(v, s),
              color = "firebrick", linetype = "solid", linewidth = 1.5) +
    geom_point(aes(x = snapshot_current$v, y = snapshot_current$s),
               size = 5, color = "black") +
    coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE, clip = "on") +
    scale_x_continuous(breaks = c(0,0.5,1), labels = c( "0", "0.5", "1")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c( "0", "0.5", "1")) +
    labs(title = title_label) +
    mytheme +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white"),
      panel.border     = element_rect(fill = NA, color = "grey20", linewidth = 1),
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
snapshots <- c(20000, 50000, 100000, 200000, 500000, 999000)
labels <- c("20K gen.", "50K gen.", "100K gen.", "200K gen.", "500K gen.", "1M gen.")

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
snapshot_figure <- (panels[[1]] | panels[[2]] | panels[[3]]) /
  (panels[[4]] + xlab(NULL) + ylab(NULL) | panels[[5]] | panels[[6]] + xlab(NULL))

ggsave("figures/Snapshots.png", snapshot_figure, width = 7.5, height = 5.5, units = "in")
ggsave("figures/Snapshots.pdf", snapshot_figure, width = 7.5, height = 5.5, units = "in")


#=========================== Nash Regions and Occupancy ==================================

# ---------- Load Data ----------
nash_region <- read.csv("/Users/canankarakoc/Documents/GitHub/GoldsteinGameTheory/organized_test_simulations/Acute_Fitness_Model/20251028_Host&Path_Reactive_m&a_0.01/postproc/fig7_region_acute.csv")
occupancy <- read.csv("/Users/canankarakoc/Documents/GitHub/GoldsteinGameTheory/organized_test_simulations/Acute_Fitness_Model/20251028_Host&Path_Reactive_m&a_0.01/postproc/fig7_occupancy.csv")

# ES persistent states
dwell_threshold <- 0.01
es_persistent <- es_acute %>%
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
nash_ocupancy <- pA_violations + pB + pC + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.position = "right",
    text = element_text(size = 14)
  )

ggsave("figures/nash_ocupancy.png", nash_ocupancy, width = 8, height = 4, units = "in")
ggsave("figures/nash_ocupancy..pdf", nash_ocupancy, width = 8, height = 4, units = "in")

#!!!! this figure is better when the simulation is recorded every generation 

#==================== Stability analysis ==============================================

# --- Plot 1: Host Strategy Parameters ---
# Thinned data!!!

# Host intercept (bS)
p_bS <- ggplot(es_acute, aes(x = gen, y = bS)) +
  geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Host strategy intercept"
  ) +
  mytheme

# Host slope (mS)
p_mS <- ggplot(es_acute, aes(x = gen, y = mS)) +
  geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Host strategy slope"
  ) +
  mytheme

# --- Plot 2: Pathogen Strategy Parameters ---

# Pathogen intercept (bV)
p_bV <- ggplot(es_acute, aes(x = gen, y = bV)) +
  geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Pathogen strategy intercept"
  ) +
  mytheme

# Pathogen slope (mV)
p_mV <- ggplot(es_acute, aes(x = gen, y = mV)) +
  geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Pathogen strategy slope"
  ) +
  mytheme

# --- Combine plots ---
combined_plot <- (p_bS | p_mS) / (p_bV | p_mV) +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(widths = c(1, 1, 1, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 16)
  )

# Save plot
ggsave("strategy_parameters_evolution.pdf", combined_plot, width = 7.5, height = 9)

# --- Additional: Plot product mS * mV (stability criterion) ---

p_product <- es_acute %>%
  mutate(product = mS * mV) %>%
  ggplot(aes(x = gen, y = product)) +
  geom_line(alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "firebrick") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  annotate("rect", xmin = 1e4, xmax = 1e6, ymin = -1, ymax = 1,
           alpha = 0.1, fill = "green") +
  annotate("text", x = 1e5, y = 0.5, 
           label = "Stable region\n(|mS·mV| < 1)", 
           size = 6, color = "darkgreen") +
  labs(
    title = "Stability Criterion: Product of Strategy Slopes",
    x = "Generation",
    y = expression(m[S] %.% m[V])
  ) +
  mytheme

ggsave("strategy_stability_product.pdf", p_product, width = 6.5, height = 5.5)

# --- Phase plot: mS vs mV ---
p_phase <- ggplot(es_acute, aes(x = mV, y = mS)) +
  geom_density_2d_filled(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "white") +
  # Add stability boundaries (hyperbolas: mS*mV = ±1)
  stat_function(fun = function(x) 1/x, 
                xlim = c(0.01, 10), 
                color = "red", linetype = "dashed", linewidth = 0.8) +
  stat_function(fun = function(x) -1/x, 
                xlim = c(0.01, 10), 
                color = "red", linetype = "dashed", linewidth = 0.8) +
  stat_function(fun = function(x) 1/x, 
                xlim = c(-10, -0.01), 
                color = "red", linetype = "dashed", linewidth = 0.8) +
  stat_function(fun = function(x) -1/x, 
                xlim = c(-10, -0.01), 
                color = "red", linetype = "dashed", linewidth = 0.8) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  labs(
    title = "Phase Space: Strategy Slopes",
    subtitle = "Red lines show stability boundaries (mS·mV = ±1)",
    x = expression(m[V]~"(Pathogen slope)"),
    y = expression(m[S]~"(Host slope)")
  ) +
  mytheme +
  theme(legend.position = "none")

# Filter to interior equilibria (away from boundaries)
# ============================================================================
# When equilibria are at boundaries (v=0,1 or s=0,1), intercepts/slopes
# become ill-defined (many combinations pass through boundary points)

es_interior <- es_acute %>%
  filter(v > 0.05, v < 0.95, s > 0.05, s < 0.95)

# Cap extreme slopes at biologically meaningful limits
# ============================================================================
# Very large slopes (>10 or <-10) all function similarly due to clamping
# The distinction between mS=-100 and mS=-1000 doesn't matter biologically

es_capped <- es_acute %>%
  mutate(
    mS_capped = pmin(pmax(mS, -10), 10),
    mV_capped = pmin(pmax(mV, -10), 10),
    bS_capped = pmin(pmax(bS, -5), 5),
    bV_capped = pmin(pmax(bV, -5), 5)
  )

# Thin data for plotting
es_interior_thin <- es_interior %>%
  mutate(row_num = row_number()) %>%
  filter(row_num %% 10 == 0)

es_capped_thin <- es_capped %>%
  mutate(row_num = row_number()) %>%
  filter(row_num %% 10 == 0)


# Interior equilibria only 
# ============================================================================

p1_bS <- ggplot(es_interior_thin, aes(x = gen, y = bS)) +
  geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Host intercept"
  ) +
  mytheme

p1_mS <- ggplot(es_interior_thin, aes(x = gen, y = mS)) +
  geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Host slope"
  ) +
  mytheme

p1_bV <- ggplot(es_interior_thin, aes(x = gen, y = bV)) +
  geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Pathogen intercept"
  ) +
  mytheme

p1_mV <- ggplot(es_interior_thin, aes(x = gen, y = mV)) +
  geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Pathogen slope"
  ) +
  mytheme

plot1 <- (p1_bS | p1_mS) / (p1_bV | p1_mV) +
  plot_annotation(
    title = "Strategy Parameters: Interior Equilibria Only",
    subtitle = "Filtered to 0.05 < v,s < 0.95 to exclude boundary artifacts",
    tag_levels = "A"
  )+
  plot_layout(widths = c(1, 1, 1, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 16)
  )

# Save plot
ggsave("strategy_parameters_evolution.pdf", plot1, width = 7.5, height = 9)

# ============================================================================
# PLOT SET 2: All data with capped values (shows boundary behavior)
# ============================================================================

p2_bS <- ggplot(es_capped_thin, aes(x = gen, y = bS_capped)) +
  geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Host Intercept (capped ±5)"
  ) +
  mytheme

p2_mS <- ggplot(es_capped_thin, aes(x = gen, y = mS_capped)) +
  geom_line(color = "steelblue", alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = c(-10, 10), linetype = "dotted", color = "firebrick", alpha = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Host Slope (capped ±10)"
  ) +
  mytheme

p2_bV <- ggplot(es_capped_thin, aes(x = gen, y = bV_capped)) +
  geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Pathogen Intercept (capped ±5)"
  ) +
  mytheme

p2_mV <- ggplot(es_capped_thin, aes(x = gen, y = mV_capped)) +
  geom_line(color = "lightcoral", alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = c(-10, 10), linetype = "dotted", color = "firebrick", alpha = 0.5) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  labs(
    x = "Generation",
    y = "Pathogen Slope (capped ±10)"
  ) +
  mytheme

plot2 <- (p2_bS | p2_mS) / (p2_bV | p2_mV) +
  plot_annotation(
    tag_levels = "A"
  )+
  plot_layout(widths = c(1, 1, 1, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 16)
  )

ggsave("figures/strategy_params_capped.pdf", plot2, width = 12, height = 8)

# ============================================================================
# PLOT 3: Stability product (mS × mV) - most important!
# ============================================================================

es_stability <- es_acute %>%
  mutate(
    product = mS * mV,
    product_capped = pmin(pmax(product, -5), 5),
    stable = abs(product) < 1,
    at_boundary = (v < 0.05 | v > 0.95 | s < 0.05 | s > 0.95)
  )

es_stability_thin <- es_stability %>%
  mutate(row_num = row_number()) %>%
  filter(row_num %% 10 == 0)

p_stab <- ggplot(es_stability, aes(x = gen, y = product_capped)) +
  geom_line(aes(color = stable), alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  annotate("rect", xmin = 1e4, xmax = 1e6, ymin = -1, ymax = 1,
           alpha = 0.05, fill = "green") +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6)) +
  scale_color_manual(
    values = c("TRUE" = "darkgreen", "FALSE" = "darkred"),
    labels = c("TRUE" = "Stable (|mS·mV| < 1)", 
               "FALSE" = "Unstable (|mS·mV| ≥ 1)")
  ) +
  labs(
    title = "Stability Criterion: mS × mV Product",
    subtitle = "Values capped at ±5 for visualization",
    x = "Generation",
    y = expression(m[S] %.% m[V]~"(capped)"),
    color = "Stability"
  ) +
  mytheme +
  theme(legend.position = c(0.25, 0.85),
        legend.background = element_rect(fill = "white", color = "black"))

ggsave("figures/strategy_stability.pdf", p_stab, width = 7, height = 6)

# ============================================================================
# PLOT 4: Phase space (mS vs mV)
# ============================================================================
#Phase Space: Strategy Slopes (Interior Equilibria)",
#"Red hyperbolas: stability boundaries (mS·mV = ±1)

p_phase_interior <- ggplot(es_interior, aes(x = mV, y = mS)) +
  geom_density_2d_filled(alpha = 0.7, bins = 100) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "white", linewidth = 0.8) +
  # Stability boundaries: mS*mV = ±1
  geom_function(fun = function(x) 1/x, xlim = c(0.1, 10), 
                color = "firebrick", linetype = "dashed", linewidth = 1) +
  geom_function(fun = function(x) -1/x, xlim = c(0.1, 10), 
                color = "firebrick", linetype = "dashed", linewidth = 1) +
  geom_function(fun = function(x) 1/x, xlim = c(-10, -0.1), 
                color = "firebrick", linetype = "dashed", linewidth = 1) +
  geom_function(fun = function(x) -1/x, xlim = c(-10, -0.1), 
                color = "firebrick", linetype = "dashed", linewidth = 1) +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(
    x = expression(m[V]~"(Pathogen slope)"),
    y = expression(m[S]~"(Host slope)")
  ) +
  mytheme +
  theme(legend.position = "none")

ggsave("figures/strategy_phase_interior.pdf", p_phase_interior, width = 5, height = 5.5)

# ============================================================================
# PLOT 5: Distinguish boundary vs interior over time
# ============================================================================

boundary_summary <- es_stability %>%
  mutate(time_bin = cut(gen, breaks = seq(0, 1e6, by = 1e4))) %>%
  group_by(time_bin) %>%
  summarise(
    gen_mid = mean(gen),
    frac_boundary = mean(at_boundary),
    frac_stable = mean(stable),
    mean_v = mean(v),
    mean_s = mean(s),
    .groups = "drop"
  ) %>%
  filter(!is.na(time_bin))

p_boundary <- ggplot(boundary_summary, aes(x = gen_mid)) +
  geom_line(aes(y = frac_boundary, color = "At boundary"), linewidth = 1) +
  geom_line(aes(y = frac_stable, color = "Stable"), linewidth = 1) +
  scale_x_log10(breaks = c(1e4, 1e5, 1e6)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(
    values = c("At boundary" = "purple", "Stable" = "darkgreen")
  ) +
  labs(
    title = "Fraction of Time at Boundaries vs in Stable Region",
    x = "Generation",
    y = "Fraction",
    color = "State"
  ) +
  mytheme +
  theme(legend.position = c(0.2, 0.8))

print(p_boundary)
ggsave("figures/boundary_vs_stable.pdf", p_boundary, width = 10, height = 6)


# SUMMARY STATISTICS
# ============================================================================

cat("  Fraction at boundaries:", 
    round(mean(es_stability$at_boundary), 3), "\n")
cat("  Fraction stable:", 
    round(mean(es_stability$stable), 3), "\n")
cat("  Fraction stable & interior:", 
    round(mean(es_stability$stable & !es_stability$at_boundary), 3), "\n\n")

cat("Interior equilibria statistics:\n")
es_interior %>%
  summarise(
    mean_bS = mean(bS),
    sd_bS = sd(bS),
    mean_mS = mean(mS),
    sd_mS = sd(mS),
    mean_bV = mean(bV),
    sd_bV = sd(bV),
    mean_mV = mean(mV),
    sd_mV = sd(mV),
    mean_product = mean(mS * mV),
    frac_stable = mean(abs(mS * mV) < 1)
  ) %>%
  print()


#===============================================================================================
# TIME SERIES ANALYSIS
#===============================================================================================
library(pracma)  # for spectral analysis
library(forecast)  # for time series analysis

# ==========================================
# 1. ORGANIZE DATA INTO TIDY FORMAT
# ==========================================

# Function to add metadata to each dataset
add_metadata <- function(df, fitness, evolution) {
  df %>%
    mutate(
      fitness = fitness,         # "regular" or "minimal"
      evolution = evolution,     # "ET-ET", "ER-ER", "ERpath-EThost", "ETpath-ERhost"
    )
}

# Combine all datasets
all_data <- bind_rows(
  # Regular acute fitness
  add_metadata(ei_acute, "regular", "ET-ET"),
  add_metadata(es_acute, "regular", "ER-ER"),
  add_metadata(es_acute_fixed_host, "regular", "ERpath-EThost"),
  add_metadata(es_acute_fixed_pathogen, "regular", "ETpath-ERhost"),
  
  # Minimal fitness
  add_metadata(ei_acute_both_reactive_minimal_fit, "minimal", "ET-ET"),
  add_metadata(es_acute_both_reactive_minimal_fit, "minimal", "ER-ER"),
  add_metadata(es_acute_fixed_host_minimal_fit, "minimal", "ERpath-EThost"),
  add_metadata(es_acute_fixed_pathogen_minimal_fit, "minimal", "ETpath-ERhost")
) %>%
  mutate(
    condition = paste(fitness, evolution, sep = "_"),
    # Convert omega columns to numeric, handling "NA" strings
    omegaPath = as.numeric(ifelse(omegaPath == "NA", NA, omegaPath)),
    omegaHost = as.numeric(ifelse(omegaHost == "NA", NA, omegaHost))
  )

# ==========================================
# 2. TIME SERIES STATISTICS FUNCTIONS
# ==========================================

# Coefficient of Variation
calc_cv <- function(x) {
  x_clean <- na.omit(x)
  if(length(x_clean) < 2) return(NA)
  sd(x_clean) / mean(x_clean)
}

# Spectral slope (noise color): 1/f^β where β indicates noise color
# β ≈ 0: white noise, β ≈ 1: pink noise, β ≈ 2: brown noise
calc_spectral_slope <- function(x, time = NULL) {
  x_clean <- na.omit(x)
  if(length(x_clean) < 50) return(NA)
  
  tryCatch({
    # Detrend
    x_detrend <- residuals(lm(x_clean ~ seq_along(x_clean)))
    
    # Compute power spectral density
    spec <- spectrum(x_detrend, plot = FALSE)
    
    # Fit linear model in log-log space (exclude DC component)
    freq <- spec$freq[-1]
    power <- spec$spec[-1]
    
    # Filter out zero/negative values
    valid <- freq > 0 & power > 0
    if(sum(valid) < 10) return(NA)
    
    fit <- lm(log10(power[valid]) ~ log10(freq[valid]))
    -coef(fit)[2]  # Return negative slope (β)
  }, error = function(e) NA)
}

# Autocorrelation at lag 1
calc_acf1 <- function(x) {
  x_clean <- na.omit(x)
  if(length(x_clean) < 10) return(NA)
  tryCatch({
    acf(x_clean, lag.max = 1, plot = FALSE)$acf[2]
  }, error = function(e) NA)
}

# Variance
calc_variance <- function(x) {
  x_clean <- na.omit(x)
  if(length(x_clean) < 2) return(NA)
  var(x_clean)
}

# Mean
calc_mean <- function(x) {
  mean(x, na.rm = TRUE)
}


# Function to calculate ACF at multiple specific lags
calc_acf_multilags <- function(x, lags = c(1, 10, 100, 1000)) {
  x_clean <- na.omit(x)
  if(length(x_clean) < max(lags) + 10) return(rep(NA, length(lags)))
  
  tryCatch({
    acf_result <- acf(x_clean, lag.max = max(lags), plot = FALSE)
    # Return ACF values at specified lags (+1 because acf includes lag 0)
    sapply(lags, function(lag) {
      if(lag + 1 <= length(acf_result$acf)) {
        acf_result$acf[lag + 1]
      } else {
        NA
      }
    })
  }, error = function(e) rep(NA, length(lags)))
}

# Function to calculate "ACF decay rate" - how fast does correlation drop?
calc_acf_decay <- function(x, max_lag = 1000) {
  x_clean <- na.omit(x)
  if(length(x_clean) < max_lag + 10) return(NA)
  
  tryCatch({
    acf_result <- acf(x_clean, lag.max = max_lag, plot = FALSE)
    acf_vals <- as.numeric(acf_result$acf[-1])  # Remove lag 0
    lags <- seq_along(acf_vals)
    
    # Fit exponential decay: ACF(lag) ~ exp(-decay_rate * lag)
    # Only use positive ACF values for fitting
    valid <- acf_vals > 0 & is.finite(acf_vals)
    if(sum(valid) < 10) return(NA)
    
    fit <- lm(log(acf_vals[valid]) ~ lags[valid])
    decay_rate <- -coef(fit)[2]  # Negative slope = decay rate
    
    return(decay_rate)
  }, error = function(e) NA)
}

# Function to find "effective correlation length" - lag where ACF drops below threshold
calc_correlation_length <- function(x, threshold = 0.1, max_lag = 5000) {
  x_clean <- na.omit(x)
  if(length(x_clean) < max_lag + 10) return(NA)
  
  tryCatch({
    acf_result <- acf(x_clean, lag.max = max_lag, plot = FALSE)
    acf_vals <- as.numeric(acf_result$acf[-1])  # Remove lag 0
    
    # Find first lag where ACF drops below threshold
    below_threshold <- which(abs(acf_vals) < threshold)
    if(length(below_threshold) > 0) {
      return(below_threshold[1])
    } else {
      return(max_lag)  # Correlation persists beyond max_lag
    }
  }, error = function(e) NA)
}

# ==========================================
# UPDATED STATISTICS WITH MULTI-LAG ACF
# ==========================================

time_series_stats_enhanced <- all_data %>%
  group_by(fitness, evolution, condition) %>%
  summarise(
    # Virulence (v)
    v_mean = calc_mean(v),
    v_cv = calc_cv(v),
    v_var = calc_variance(v),
    v_spectral_slope = calc_spectral_slope(v),
    v_acf1 = calc_acf1(v),
    v_acf10 = calc_acf_multilags(v, lags = 10)[1],
    v_acf100 = calc_acf_multilags(v, lags = 100)[1],
    v_acf1000 = calc_acf_multilags(v, lags = 1000)[1],
    v_acf_decay = calc_acf_decay(v),
    v_corr_length = calc_correlation_length(v)*100,
    
    # Clearance (c)
    c_mean = calc_mean(s),
    c_cv = calc_cv(s),
    c_var = calc_variance(s),
    c_spectral_slope = calc_spectral_slope(s),
    c_acf1 = calc_acf1(s),
    c_acf10 = calc_acf_multilags(s, lags = 10)[1],
    c_acf100 = calc_acf_multilags(s, lags = 100)[1],
    c_acf1000 = calc_acf_multilags(s, lags = 1000)[1],
    c_acf_decay = calc_acf_decay(s),
    c_corr_length = calc_correlation_length(s)*100,
    
    # Pathogen fitness
    pathFit_mean = calc_mean(pathFit),
    pathFit_cv = calc_cv(pathFit),
    pathFit_var = calc_variance(pathFit),
    pathFit_spectral_slope = calc_spectral_slope(pathFit),
    pathFit_acf1 = calc_acf1(pathFit),
    pathFit_acf10 = calc_acf_multilags(pathFit, lags = 10)[1],
    pathFit_acf100 = calc_acf_multilags(pathFit, lags = 100)[1],
    pathFit_acf1000 = calc_acf_multilags(pathFit, lags = 1000)[1],
    pathFit_acf_decay = calc_acf_decay(pathFit),
    pathFit_corr_length = calc_correlation_length(pathFit)*100,
    
    # Host fitness
    hostFit_mean = calc_mean(hostFit),
    hostFit_cv = calc_cv(hostFit),
    hostFit_var = calc_variance(hostFit),
    hostFit_spectral_slope = calc_spectral_slope(hostFit),
    hostFit_acf1 = calc_acf1(hostFit),
    hostFit_acf10 = calc_acf_multilags(hostFit, lags = 10)[1],
    hostFit_acf100 = calc_acf_multilags(hostFit, lags = 100)[1],
    hostFit_acf1000 = calc_acf_multilags(hostFit, lags = 1000)[1],
    hostFit_acf_decay = calc_acf_decay(hostFit),
    hostFit_corr_length = calc_correlation_length(hostFit)*100,
    
    # Omega pathogen
    omegaPath_mean = calc_mean(omegaPath),
    omegaPath_cv = calc_cv(omegaPath),
    omegaPath_var = calc_variance(omegaPath),
    omegaPath_spectral_slope = calc_spectral_slope(omegaPath),
    omegaPath_acf1 = calc_acf1(omegaPath),
    omegaPath_acf10 = calc_acf_multilags(omegaPath, lags = 10)[1],
    omegaPath_acf100 = calc_acf_multilags(omegaPath, lags = 100)[1],
    omegaPath_acf1000 = calc_acf_multilags(omegaPath, lags = 1000)[1],
    omegaPath_acf_decay = calc_acf_decay(omegaPath),
    omegaPath_corr_length = calc_correlation_length(omegaPath)*100,
    
    # Omega host
    omegaHost_mean = calc_mean(omegaHost),
    omegaHost_cv = calc_cv(omegaHost),
    omegaHost_var = calc_variance(omegaHost),
    omegaHost_spectral_slope = calc_spectral_slope(omegaHost),
    omegaHost_acf1 = calc_acf1(omegaHost),
    omegaHost_acf10 = calc_acf_multilags(omegaHost, lags = 10)[1],
    omegaHost_acf100 = calc_acf_multilags(omegaHost, lags = 100)[1],
    omegaHost_acf1000 = calc_acf_multilags(omegaHost, lags = 1000)[1],
    omegaHost_acf_decay = calc_acf_decay(omegaHost),
    omegaHost_corr_length = calc_correlation_length(omegaHost)*100,
    
    n_obs = n(),
    .groups = "drop"
  )

# Save enhanced statistics
write.csv(time_series_stats_enhanced, 
          "time_series_statistics_enhanced.csv", 
          row.names = FALSE)

# ==========================================
# 4. VISUALIZATION - COMPARISON PLOTS
# ==========================================

# Define labels once
variable_labels <- c(
  "v" = "italic(v)",
  "c" = "italic(c)",
  "pathFit" = "W[P]",
  "hostFit" = "W[H]",
  "omegaPath" = "omega[P]",
  "omegaHost" = "omega[H]"
)

# Updated variable labels (no tilde needed, just space)
variable_labels_full <- c(
  "v" = "italic(v)~'(virulence)'",
  "c" = "italic(c)~'(clearance)'"
)

# Reshape for plotting
stats_long <- time_series_stats_enhanced %>%
  pivot_longer(
    cols = -c(fitness, evolution, condition, n_obs),
    names_to = "metric",
    values_to = "value"
  ) %>%
  extract(
    metric,
    into = c("variable", "stat_type"),
    regex = "^(v|c|pathFit|hostFit|omegaPath|omegaHost)_(.+)$"
  ) %>%
  mutate(
    variable = factor(variable, 
                      levels = c("v", "c", "pathFit", "hostFit", "omegaPath", "omegaHost")),
    evolution = factor(evolution, 
                       levels = c("ET-ET", "ERpath-EThost", "ETpath-ERhost", "ER-ER"))
  )

# Select key metrics and filter to traits only
key_metrics <- c("mean", "cv", "spectral_slope", "acf1", "corr_length")

stats_traits <- stats_long %>%
  filter(
    variable %in% c("v", "c"),
    stat_type %in% key_metrics
  ) %>%
  mutate(
    stat_type = factor(
      stat_type,
      levels = key_metrics,
      labels = c(
        "Mean",
        "CV",
        "Noise color",
        "ACF100",
        "Memory"
      )
    )
  )

# Check what you have
table(stats_traits$stat_type)

# ==========================================
# FIGURE WITH PROPER REFERENCE LINES
# ==========================================

p_traits_summary <- stats_traits %>%
  ggplot(aes(x = evolution, y = value, fill = fitness)) +
  geom_col(position = "dodge", width = 0.7) +
  
  # Reference lines for spectral slope (noise color)
  geom_hline(
    data = data.frame(
      stat_type = factor("Noise color", 
                         levels = levels(stats_traits$stat_type)),
      yint = c(0, 1, 2)
    ),
    aes(yintercept = yint),
    linetype = "dotted",
    color = "gray50",
    alpha = 0.5
  ) +
  
  # Reference line for ACF (zero correlation)
  geom_hline(
    data = data.frame(
      stat_type = factor("ACF100", 
                         levels = levels(stats_traits$stat_type)),
      yint = 0
    ),
    aes(yintercept = yint),
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5
  ) +
  
  # Faceting
  facet_grid(
    stat_type ~ variable,
    scales = "free_y",
    labeller = labeller(
      variable = as_labeller(variable_labels_full, default = label_parsed)
    )
  ) +
  
  # Styling
  scale_fill_manual(
    values = c("minimal" = "#E8E8E8", "regular" = "#404040"),
    name = "Fitness landscape"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  
  labs(
    x = "Evolution mode",
    y = "Value"
  ) +
  
  mytheme +
  theme(
    strip.text.y = element_text(size = 13, margin = margin(r = 5, l = 5)),
    strip.text.x = element_text(size = 13),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    panel.spacing = unit(0.6, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )

p_traits_summary

# Save with adequate dimensions
ggsave(
  "traits_timeseries_summary.pdf",
  p_traits_summary,
  width = 5.5,
  height = 9)  # Taller to accommodate 5 rows

# ==========================================
# PLOT REALIZED STATES ON FITNESS LANDSCAPE
# ==========================================
# ==========================================
# CORE FITNESS FUNCTIONS
# ==========================================

# Regular (acute) fitness
d_acute <- function(c, v) {
  d_0 + m_c * (((1+epsilon) * c) / ((1+epsilon) - c)) + 
    m_v * (((1+epsilon) * v )/ ((1+epsilon) - v))
}

r_acute <- function(v) v^beta

f_H_acute <- function(c, v) {
  d <- d_acute(c, v)
  c / (c + d)
}

f_P_acute <- function(c, v) {
  d <- d_acute(c, v)
  r <- r_acute(v)
  r / (c + d)
}

# Minimal fitness (simple polynomial)
f_H_minimal <- function(c, v) {
  c * (1 - c) * (1 - v)
}

f_P_minimal <- function(c, v) {
  v * (1 - v) * (1 - c)
}

# ==========================================
# HELPER: Calculate hex counts for global limits
# ==========================================

calculate_hex_counts <- function(data, thin_factor = 1, nbins = 100) {
  data_thin <- data %>%
    filter(row_number() %% thin_factor == 0) %>%
    filter(!is.na(v), !is.na(s), v >= 0, v <= 1, s >= 0, s <= 1)
  
  if (nrow(data_thin) < 10) return(NULL)
  
  hex_data <- hexbin::hexbin(data_thin$v, data_thin$s, xbins = nbins)
  hex_data@count
}

# ==========================================
# CALCULATE GLOBAL LIMITS PER FITNESS MODEL
# ==========================================

calculate_global_limits <- function(data, fitness_type) {
  cat(paste("\nCalculating global count range for", fitness_type, "fitness...\n"))
  
  all_counts <- list()
  
  for (evo_mode in c("ET-ET", "ER-ER", "ERpath-EThost", "ETpath-ERhost")) {
    data_subset <- data %>% 
      filter(evolution == evo_mode, fitness == fitness_type)
    
    if (nrow(data_subset) > 0) {
      counts <- calculate_hex_counts(data_subset, thin_factor = 1, nbins = 100)
      if (!is.null(counts)) {
        all_counts <- c(all_counts, list(counts))
      }
    }
  }
  
  if (length(all_counts) == 0) {
    return(list(limits = c(1, 1000), breaks = c(1, 10, 100, 1000)))
  }
  
  all_counts_vec <- unlist(all_counts)
  global_min <- max(1, min(all_counts_vec, na.rm = TRUE))
  global_max <- max(all_counts_vec, na.rm = TRUE)
  
  legend_breaks <- c(1, 10, 100, 1000, 10000)
  legend_breaks <- legend_breaks[legend_breaks >= global_min & legend_breaks <= global_max * 1.5]
  
  cat("Range:", global_min, "to", global_max, "\n")
  cat("Breaks:", paste(legend_breaks, collapse = ", "), "\n")
  
  return(list(limits = c(global_min, global_max), breaks = legend_breaks))
}

# ==========================================
# PLOTTING FUNCTION (BACK TO GEOM_HEX)
# ==========================================

plot_states_on_landscape <- function(data, 
                                     fitness_type = "regular",
                                     thin_factor = 1,
                                     nbins = 100,
                                     count_limits = c(1, 10000),
                                     count_breaks = c(1, 10, 100, 1000, 10000),
                                     show_x_axis = TRUE,
                                     show_y_axis = TRUE,
                                     show_x_title = TRUE,
                                     show_y_title = TRUE,
                                     show_nash = TRUE) {
  
  # Thin data
  data_thin <- data %>%
    filter(row_number() %% thin_factor == 0) %>%
    filter(!is.na(v), !is.na(s), v >= 0, v <= 1, s >= 0, s <= 1)
  
  if (nrow(data_thin) < 10) {
    warning("Insufficient data points")
    return(ggplot() + theme_void())
  }
  
  # Select fitness functions
  if (fitness_type == "minimal") {
    f_H <- f_H_minimal
    f_P <- f_P_minimal
    cache_name <- "br_cache_minimal"
  } else {
    f_H <- f_H_acute
    f_P <- f_P_acute
    cache_name <- "br_cache_regular"
  }
  
  # Create fitness grid
  v_grid <- seq(0, 1, length.out = 150)
  c_grid <- seq(0, 1, length.out = 150)
  
  fitness_grid <- expand.grid(v = v_grid, c = c_grid) %>%
    mutate(
      hostFit = f_H(c, v),
      pathFit = f_P(c, v)
    )
  
  # Calculate or retrieve cached best response curves
  if (!exists(cache_name, envir = .GlobalEnv)) {
    cat(paste("Calculating best response curves for", fitness_type, "...\n"))
    
    v_fine <- seq(0.01, 0.99, length.out = 500)
    c_fine <- seq(0.01, 0.99, length.out = 500)
    
    c_opt <- sapply(v_fine, function(v) {
      optimize(function(c) -f_H(c, v), c(0.01, 0.99))$minimum
    })
    
    v_opt <- sapply(c_fine, function(c) {
      optimize(function(v) -f_P(c, v), c(0.01, 0.99))$minimum
    })
    
    br_cache <- list(
      host = tibble(v = v_fine, c = c_opt),
      path = tibble(c = c_fine, v = v_opt)
    )
    
    assign(cache_name, br_cache, envir = .GlobalEnv)
  } else {
    br_cache <- get(cache_name, envir = .GlobalEnv)
  }
  
  # Calculate Nash equilibrium
  nash_point <- NULL
  if (show_nash) {
    v_seq <- seq(0.05, 0.95, length.out = 300)
    c_opt_seq <- sapply(v_seq, function(v) {
      optimize(function(c) -f_H(c, v), c(0.01, 0.99))$minimum
    })
    v_opt_seq <- sapply(c_opt_seq, function(c) {
      optimize(function(v) -f_P(c, v), c(0.01, 0.99))$minimum
    })
    
    idx <- which.min(abs(v_opt_seq - v_seq))
    nash_point <- tibble(v = v_seq[idx], c = c_opt_seq[idx])
  }
  
  # Mean state
  mean_state <- data_thin %>% 
    summarise(v = mean(v), s = mean(s))
  
  # Build plot
  p <- ggplot() +
    # Fitness contours
    geom_contour(
      data = fitness_grid, 
      aes(v, c, z = pathFit), 
      color = "lightcoral", alpha = 0.3, bins = 12, linewidth = 0.5
    ) +
    geom_contour(
      data = fitness_grid, 
      aes(v, c, z = hostFit), 
      color = "steelblue", alpha = 0.3, bins = 12, linewidth = 0.5
    ) +
    
    # Best response curves
    geom_line(
      data = br_cache$host, 
      aes(v, c), 
      color = "steelblue", linewidth = 1, linetype = "dashed"
    ) +
    geom_line(
      data = br_cache$path, 
      aes(v, c), 
      color = "lightcoral", linewidth = 1, linetype = "dashed"
    ) +
    
    # Realized states (hex density) - BACK TO ORIGINAL
    geom_hex(
      data = data_thin, 
      aes(x = v, y = s), 
      bins = nbins, 
      alpha = 0.7
    ) +
    scale_fill_viridis_c(
      option = "plasma", 
      name = "Count", 
      trans = "log10",
      limits = count_limits,
      breaks = count_breaks,
      labels = scales::comma,
      oob = scales::squish
    ) +
    
    # Nash equilibrium (black dot)
    {if (show_nash && !is.null(nash_point)) 
      geom_point(
        data = nash_point, 
        aes(x = v, y = c),
        size = 4, color = "black", shape = 16
      )
    } +
    
    # Mean state (gold X)
    geom_point(
      data = mean_state, 
      aes(x = v, y = s), 
      size = 4, color = "gold", shape = 4, stroke = 1.5
    ) +
    
    # Scales and theme
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), labels = c("0", "0.5", "1")) +
    coord_fixed() +
    mytheme +
    theme(legend.position = "none")
  
  # Conditionally show/hide axes
  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  if (!show_y_axis) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  if (!show_x_title) {
    p <- p + labs(x = NULL)
  } else {
    p <- p + labs(x = "v (virulence)")
  }
  if (!show_y_title) {
    p <- p + labs(y = NULL)
  } else {
    p <- p + labs(y = "c (clearance)")
  }
  
  return(p)
}

# ==========================================
# CREATE 2×4 COMBINED FIGURE
# ==========================================

create_combined_landscape_figure <- function(data,
                                             tag_position_x = 0.08,
                                             tag_position_y = 0.95,
                                             tag_size = 18) {
  
  cat("\n=== Creating combined 2×4 landscape figure ===\n")
  
  # Calculate global limits for each fitness type
  regular_limits <- calculate_global_limits(data, "regular")
  minimal_limits <- calculate_global_limits(data, "minimal")
  
  # Clear both caches
  if (exists("br_cache_minimal", envir = .GlobalEnv)) {
    rm(br_cache_minimal, envir = .GlobalEnv)
  }
  if (exists("br_cache_regular", envir = .GlobalEnv)) {
    rm(br_cache_regular, envir = .GlobalEnv)
  }
  
  # Define evolution modes
  evo_modes <- c("ET-ET", "ERpath-EThost", "ETpath-ERhost", "ER-ER")
  
  # Create all 8 plots
  plots <- list()
  
  # TOP ROW: Regular fitness (A-D)
  for (i in 1:4) {
    plots[[i]] <- plot_states_on_landscape(
      data = data %>% filter(evolution == evo_modes[i], fitness == "regular"),
      fitness_type = "regular",
      count_limits = regular_limits$limits,
      count_breaks = regular_limits$breaks,
      show_x_axis = FALSE,
      show_y_axis = (i == 1),
      show_x_title = FALSE,
      show_y_title = (i == 1),
      show_nash = TRUE
    )
  }
  
  # BOTTOM ROW: Minimal fitness (E-H)
  for (i in 1:4) {
    plots[[i + 4]] <- plot_states_on_landscape(
      data = data %>% filter(evolution == evo_modes[i], fitness == "minimal"),
      fitness_type = "minimal",
      count_limits = minimal_limits$limits,
      count_breaks = minimal_limits$breaks,
      show_x_axis = TRUE,
      show_y_axis = (i == 1),
      show_x_title = (i == 1),
      show_y_title = (i == 1),
      show_nash = TRUE
    )
  }
  
  # Combine with patchwork
  combined <- (plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) /
    (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]]) +
    plot_annotation(
      caption = "Black dot = Nash equilibrium | Gold X = mean state",
      tag_levels = 'A',
      theme = theme(
        plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 10)),
        plot.tag = element_text(size = tag_size, face = "bold")
      )
    ) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      legend.key.height = unit(2, "cm"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.tag.position = c(tag_position_x, tag_position_y),
      plot.tag = element_text(size = tag_size, face = "bold")
    )
  
  return(combined)
}

# ==========================================
# GENERATE FIGURE
# ==========================================

p_combined <- create_combined_landscape_figure(
  all_data,
  tag_position_x = 0.01,   # Adjust position
  tag_position_y = 1.05,   # Adjust position  
  tag_size = 18            # Bold, larger labels
)

print(p_combined)

# Save
ggsave(
  "landscape_combined_8panel.pdf", 
  p_combined, 
  width = 10, 
  height = 6.6
)

