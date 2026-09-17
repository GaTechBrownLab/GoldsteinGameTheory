# Test fitness functions 

suppressMessages({library(ggplot2); library(patchwork)})

## ---- Fitness functions -------------------------------------------------
# v = pathogen trait (virulence), c = host trait
# k scales the "best-response tracking" term for each player
WH <- function(c, v, k) c*(1-c)*(1-v) + k * c^2*(1-c)*v
WP <- function(c, v, k) v*(1-v)*(1-c) + k * v^2*(1-v)*c

## ---- Numerical best responses (argmax over [0,1], handles boundaries) --
own <- seq(0, 1, length.out = 1001)
opp <- seq(0, 1, length.out = 401)

host_BR <- function(k) {                       # c*(v): best host trait given v
  c_star <- sapply(opp, function(v) own[which.max(WH(own, v, k))])
  data.frame(v = opp, c = c_star)
}
path_BR <- function(k) {                       # v*(c): best pathogen trait given c
  v_star <- sapply(opp, function(c) own[which.max(WP(c, own, k))])
  data.frame(c = opp, v = v_star)
}

## ---- Nash equilibria = intersections of the two BR curves --------------
nash_points <- function(k) {
  Cfun <- approxfun(host_BR(k)$v, host_BR(k)$c, rule = 2)   # c as fn of v
  Vfun <- approxfun(path_BR(k)$c, path_BR(k)$v, rule = 2)   # v as fn of c
  phi  <- function(v) Vfun(Cfun(v)) - v
  vg   <- seq(1e-4, 1-1e-4, length.out = 4000)
  pv   <- phi(vg)
  idx  <- which(pv[-1] * pv[-length(pv)] < 0)
  roots <- numeric(0)
  for (i in idx) {
    r <- tryCatch(uniroot(phi, c(vg[i], vg[i+1]))$root, error = function(e) NA)
    if (!is.na(r)) roots <- c(roots, r)
  }
  if (length(roots) == 0) return(data.frame(v = numeric(0), c = numeric(0)))
  roots <- unique(round(roots, 4))
  data.frame(v = roots, c = Cfun(roots))
}

## ======================================================================
## FIGURE 1 — best-response geometry across k
## ======================================================================
ks <- c(0, 0.5, 1, 2, 4)
br_df <- do.call(rbind, lapply(ks, function(k) {
  rbind(transform(host_BR(k), player = "Host best response  c*(v)", k = k),
        transform(path_BR(k), player = "Pathogen best response  v*(c)", k = k))
}))
nash_df <- do.call(rbind, lapply(ks, function(k) {
  np <- nash_points(k); if (nrow(np)) transform(np, k = k) else NULL
}))
br_df$kf   <- factor(br_df$k,   labels = paste0("k = ", ks))
nash_df$kf <- factor(nash_df$k, labels = paste0("k = ", ks[ks %in% nash_df$k]))

mytheme <- theme_bw(base_size = 13) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
        strip.background = element_blank(), strip.text = element_text(size = 13),
        legend.title = element_blank(), legend.position = "bottom",
        panel.grid.minor = element_blank())

p1 <- ggplot(br_df, aes(v, c, colour = player)) +
  geom_path(linewidth = 1.1) +
  geom_point(data = nash_df, aes(v, c), inherit.aes = FALSE,
             colour = "black", fill = "gold", shape = 21, size = 3.2, stroke = 0.8) +
  facet_wrap(~ kf, nrow = 1) +
  scale_colour_manual(values = c("#1b6ca8", "#c0392b")) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "virulence, v", y = "clearence, c") +
  mytheme
ggsave("/home/claude/fig1_best_response.png", p1, width = 12, height = 4.2, dpi = 150)

## ======================================================================
## FIGURE 2 — fitness landscapes (W_H, W_P) at k = 0 vs k = 2, BR overlaid
## ======================================================================
grid <- expand.grid(v = seq(0,1,length.out = 161), c = seq(0,1,length.out = 161))
land <- function(k) {
  rbind(data.frame(grid, W = WH(grid$c, grid$v, k), who = "Host fitness  W_H", k = k),
        data.frame(grid, W = WP(grid$c, grid$v, k), who = "Pathogen fitness  W_P", k = k))
}
land_df <- rbind(land(0), land(2))
land_df$kf <- factor(land_df$k, labels = c("k = 0", "k = 2"))

# matching BR overlays
ov_host <- do.call(rbind, lapply(c(0,2), function(k) transform(host_BR(k), k=k, who="Host fitness  W_H")))
ov_path <- do.call(rbind, lapply(c(0,2), function(k) transform(path_BR(k), k=k, who="Pathogen fitness  W_P")))
ov_host$kf <- factor(ov_host$k, labels=c("k = 0","k = 2"))
ov_path$kf <- factor(ov_path$k, labels=c("k = 0","k = 2"))

p2 <- ggplot(land_df, aes(v, c, fill = W)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = W), colour = "white", alpha = 0.35, linewidth = 0.3, bins = 10) +
  geom_path(data = ov_host, aes(v, c), inherit.aes = FALSE, colour = "black", linewidth = 0.9) +
  geom_path(data = ov_path, aes(v, c), inherit.aes = FALSE, colour = "black",
            linetype = "22", linewidth = 0.9) +
  facet_grid(who ~ kf) +
  scale_fill_viridis_c(option = "magma") +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "pathogen trait  v", y = "host trait  c",
       title = "Fitness landscapes  (solid = host BR, dashed = pathogen BR)") +
  mytheme + theme(legend.position = "right", legend.title = element_text(size = 11))
ggsave("/home/claude/fig2_landscapes.png", p2, width = 9, height = 8, dpi = 150)

## ---- report Nash points to console ------------------------------------
cat("\nNash equilibria found:\n")
for (k in ks) {
  np <- nash_points(k)
  if (nrow(np)) for (i in seq_len(nrow(np)))
    cat(sprintf("  k=%.1f : (v*, c*) = (%.3f, %.3f)\n", k, np$v[i], np$c[i]))
  else cat(sprintf("  k=%.1f : none interior\n", k))
}
cat("done\n")




suppressMessages({library(ggplot2)})

## ===== Symmetric (·c) coevolutionary model ==============================
## v = pathogen virulence, c = host clearance
##   W_H = c(1-c)(1-v) + k c^2 (1-c) v
##   W_P = v(1-v)(1-c) + k v^2 (1-v) c
## -----------------------------------------------------------------------

WH <- function(c, v, k) c*(1-c)*(1-v) + k * c^2*(1-c)*v
WP <- function(c, v, k) v*(1-v)*(1-c) + k * v^2*(1-v)*c

## ---- Analytic partial derivatives (for FOC solve + IFT slopes) ---------
WH_c  <- function(c, v, k) (1-v)*(1-2*c) + k*v*(2*c - 3*c^2)        # dW_H/dc
WH_cc <- function(c, v, k) -2*(1-v)      + k*v*(2 - 6*c)           # d2W_H/dc2
WH_cv <- function(c, v, k) -(1-2*c)      + k*(2*c - 3*c^2)         # d2W_H/dc dv

WP_v  <- function(c, v, k) (1-c)*(1-2*v) + k*c*(2*v - 3*v^2)        # dW_P/dv
WP_vv <- function(c, v, k) -2*(1-c)      + k*c*(2 - 6*v)           # d2W_P/dv2
WP_vc <- function(c, v, k) -(1-2*v)      + k*(2*v - 3*v^2)         # d2W_P/dv dc

## ---- Numerical best responses (for plotting the curves) ----------------
own <- seq(0, 1, length.out = 1001)
opp <- seq(0, 1, length.out = 401)
host_BR <- function(k) data.frame(v = opp,
                                  c = sapply(opp, function(v) own[which.max(WH(own, v, k))]))
path_BR <- function(k) data.frame(c = opp,
                                  v = sapply(opp, function(c) own[which.max(WP(c, own, k))]))

## ---- Interior Nash: Newton on FOC system F = (W_H,c , W_P,v) = 0 -------
## seeded from the numerical BR intersection, polished analytically
nash_solve <- function(k) {
  Cf <- approxfun(host_BR(k)$v, host_BR(k)$c, rule = 2)
  Vf <- approxfun(path_BR(k)$c, path_BR(k)$v, rule = 2)
  phi <- function(v) Vf(Cf(v)) - v
  vg <- seq(1e-3, 1-1e-3, length.out = 2000); pv <- phi(vg)
  i  <- which(pv[-1]*pv[-length(pv)] < 0)[1]
  if (is.na(i)) return(c(c = NA, v = NA))
  v0 <- uniroot(phi, c(vg[i], vg[i+1]))$root
  x  <- c(Cf(v0), v0)                                  # (c, v) seed
  for (it in 1:50) {
    Fx <- c(WH_c(x[1], x[2], k), WP_v(x[1], x[2], k))
    J  <- matrix(c(WH_cc(x[1],x[2],k), WH_cv(x[1],x[2],k),
                   WP_vc(x[1],x[2],k), WP_vv(x[1],x[2],k)), 2, 2, byrow = TRUE)
    step <- solve(J, Fx); x <- x - step
    x <- pmin(pmax(x, 1e-9), 1 - 1e-9)
    if (max(abs(step)) < 1e-12) break
  }
  c(c = x[1], v = x[2])
}

## ---- Slopes via the implicit function theorem, at the Nash -------------
analyse <- function(k) {
  ns <- nash_solve(k); cs <- ns["c"]; vs <- ns["v"]
  host_slope <- -WH_cv(cs, vs, k) / WH_cc(cs, vs, k)   # dc*/dv  (host BR)
  path_slope <- -WP_vc(cs, vs, k) / WP_vv(cs, vs, k)   # dv*/dc  (pathogen BR)
  prod <- host_slope * path_slope
  data.frame(k = k, v_star = unname(vs), c_star = unname(cs),
             host_slope = unname(host_slope), path_slope = unname(path_slope),
             product = unname(prod), stable = abs(prod) < 1)
}

ks  <- c(0, 0.5, 1, 2, 4, 6)
tab <- do.call(rbind, lapply(ks, analyse))
print(round(transform(tab, stable = tab$stable),
            4), row.names = FALSE)

## ===== Plot: stability product vs k ====================================
mytheme <- theme_bw(base_size = 13) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
        panel.grid.minor = element_blank())

p <- ggplot(tab, aes(k, abs(product))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1, colour = "#1b6ca8") +
  geom_point(shape = 21, size = 3.5, stroke = 0.7) +
  #scale_fill_manual(values = c(`TRUE` = "gold", `FALSE` = "#c0392b"),
  #                  name = "|product| < 1") +
  #annotate("text", x = max(ks)*0.7, y = 1.05, label = "instability threshold",
   #        colour = "grey40", size = 4) +
  labs(x = "tracking strength  k",
       y = expression("|"*dc^"*"/dv %.% dv^"*"/dc*"|")) +
  mytheme
ggsave("/home/claude/stability_vs_k.png", p, width = 7, height = 5, dpi = 150)
cat("\nsaved stability_vs_k.png\n")