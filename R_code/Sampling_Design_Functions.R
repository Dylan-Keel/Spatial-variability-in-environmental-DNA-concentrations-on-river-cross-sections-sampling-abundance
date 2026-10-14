## Sample size calculations for estimating the mean eDNA concentration on a
## river cross section. Used by Main_Analysis_and_Figures.R (Figures 5, 7, 8
## and S14) and by Sampling_Design_App.R.
##
## Estimand: the mean eDNA concentration on the sampled cross section.
##
## Field model (stationary; lateral position only):
##   x_i = mu + u(s_i) + e_i
##   Var(u) = p sigma^2,  p = 1 - nugget / sill  (spatial proportion of variance)
##   Var(e) = (1 - p) sigma^2
##   Cor(u(s), u(t)) = c(|s - t|), Matern correlation with nu = 1 and
##   kappa = sqrt(8) / range (correlation about 0.14 at distance = range).
##
## The within-section variability measured in a pilot study, CV_site, has
## variance sigma_w^2 = sigma^2 (1 - p cbar_WW), where cbar_WW is the mean
## correlation between two random positions on a section of width W.
##
## Precision of the sample mean as an estimate of the section mean:
##   simple random placement:    CV(n) = CV_site / sqrt(n)                              (n_ind)
##   spatially stratified:       CV(n) = CV_site sqrt([p A(n, W) + (1 - p) / n] / (1 - p cbar_WW))   (n_strat)
##       A(n, W) = (1 / n^2) sum_ij c(s_i, s_j) - (2 / n) sum_i g(s_i) + cbar_WW,
##       s_i = (i - 1/2) W / n (stratum centres), g(s) = (1 / W) int_0^W c(|s - t|) dt;
##       A tends to zero with n, so evenly spaced samples converge on the section mean.
##
## For comparison, the design-effect adjustment for clustered samples gives the
## precision of the same sample mean as an estimate of the reach mean, where
## the deviation shared by all samples on a section is error:
##   CV(n) = CV_site sqrt((1 + (n - 1) rho) / (n (1 - rho))),  rho = p cbar_WW   (n_adj)
##   with the limit CV_site sqrt(rho / (1 - rho)) as n increases.
##
## All two-dimensional integrals reduce to one-dimensional integrals of c(u),
## tabulated once per correlation range.

matern_cor <- function(D, kappa) {
  x <- kappa * D
  ifelse(x < 1e-10, 1, x * besselK(x, 1))
}

# Tabulate c(u), G(u) = int_0^u c, H(u) = int_0^u t c(t) dt on a fine grid
matern_tables <- function(range_m, u_max = 200, du = 0.05) {
  kappa <- sqrt(8) / range_m
  u <- seq(0, u_max, by = du)
  cu <- matern_cor(u, kappa)
  trap <- function(y) c(0, cumsum((y[-1] + y[-length(y)]) / 2 * du))
  list(u = u, c = cu, G = trap(cu), H = trap(u * cu), range_m = range_m)
}
.interp <- function(tab, what, x) approx(tab$u, tab[[what]], xout = pmin(x, max(tab$u)), rule = 2)$y

# mean correlation between two uniform random positions on [0, W]
cbar_random <- function(W, tab) {
  (2 / W^2) * (W * .interp(tab, "G", W) - .interp(tab, "H", W))
}
# g(s) = mean correlation between position s and a uniform random position on [0, W]
g_point <- function(s, W, tab) (.interp(tab, "G", s) + .interp(tab, "G", W - s)) / W

# A(n, W) for stratum-centre positions (vector over n)
A_stratified <- function(n_vec, W, tab) {
  cWW <- cbar_random(W, tab)
  vapply(n_vec, function(n) {
    s <- W / n
    k <- seq_len(max(n - 1, 0))
    pair_mean <- if (n > 1) (1 + 2 * sum((1 - k / n) * .interp(tab, "c", k * s))) / n else 1
    si <- (seq_len(n) - 0.5) * s
    pair_mean - 2 * mean(g_point(si, W, tab)) + cWW
  }, numeric(1))
}

# A for any set of lateral positions on the section (the samples as collected)
A_positions <- function(xs, W, tab) {
  xs <- pmin(pmax(xs, 0), W)
  D  <- abs(outer(xs, xs, "-"))
  mean(.interp(tab, "c", as.vector(D))) - 2 * mean(g_point(xs, W, tab)) + cbar_random(W, tab)
}

# ---- CV(n) under each design / estimand -------------------------------------
cv_random      <- function(cv_site, n) cv_site / sqrt(n)
cv_stratified  <- function(cv_site, n, W, p, tab) {
  if (p <= 0) return(cv_site / sqrt(n))
  cWW <- cbar_random(W, tab)
  A   <- A_stratified(n, W, tab)
  cv_site * sqrt((p * A + (1 - p) / n) / (1 - p * cWW))
}
cv_positions   <- function(cv_site, xs, W, p, tab) {
  n <- length(xs)
  if (p <= 0) return(cv_site / sqrt(n))
  cv_site * sqrt((p * A_positions(xs, W, tab) + (1 - p) / n) / (1 - p * cbar_random(W, tab)))
}
icc_section    <- function(W, p, tab) p * cbar_random(W, tab)

# ---- n required ---------------------------------------------------------------
n_random <- function(cv_site, cv_target) ceiling((cv_site / cv_target)^2)

# First n at which CV(n) <= target, scanning a running minimum; NA if not within n_max
.first_crossing <- function(f_ratio, thresholds) {   # f_ratio = CV(n)/CV_site for n = 1..n_max
  fm <- cummin(f_ratio)
  vapply(thresholds, function(tt) { i <- which(fm <= tt)[1]; if (is.na(i)) NA_real_ else as.numeric(i) }, numeric(1))
}
n_stratified <- function(cv_site, cv_target, W, p, tab, n_max = 500) {
  f <- cv_stratified(1, seq_len(n_max), W, p, tab)
  .first_crossing(f, cv_target / cv_site)
}

# Grid of required n over CV_site x width.  design = "random" | "stratified" | "reach"
n_required_grid <- function(cv_site, width, cv_target, p, tab,
                            design = c("stratified", "random", "reach"), n_max = 500) {
  design <- match.arg(design)
  out <- lapply(width, function(W) {
    n <- switch(design,
      random     = { v <- n_random(cv_site, cv_target); v[v > n_max] <- NA; v },
      stratified = n_stratified(cv_site, cv_target, W, p, tab, n_max),
      reach      = n_reach(cv_site, cv_target, W, p, tab, n_max))
    data.frame(cv_site = cv_site, width = W, n = n)
  })
  res <- do.call(rbind, out)
  attr(res, "n_max") <- n_max
  res
}

# ---- Contour panel ------------------------------------------------------------
n_breaks_default <- c(0, 3, 5, 10, 15, 20, 30, 50, 100, 250, Inf)
n_band_labels <- function(br = n_breaks_default) {
  lo <- head(br, -1); hi <- tail(br, -1)
  ifelse(is.infinite(hi), paste0("> ", lo),
         ifelse(lo == 0, paste0("≤ ", hi), paste0(lo + 1, "–", hi)))
}

contour_panel <- function(grid, title, br = n_breaks_default, mark = NULL,
                          base_theme = NULL, others = NULL) {
  library(ggplot2)
  labs_b <- n_band_labels(br)
  grid$band <- cut(grid$n, breaks = br, labels = labs_b, right = TRUE, include.lowest = TRUE)
  na_part <- grid[is.na(grid$n), ]
  ## an off-canvas tile for every band so each one gets a legend key
  keys <- data.frame(cv_site = min(grid$cv_site) - 10, width = min(grid$width) - 10, band = factor(labs_b, levels = labs_b))
  g <- ggplot(grid, aes(x = cv_site, y = width)) +
    geom_tile(data = keys, aes(fill = band), key_glyph = "rect") +
    geom_raster(data = grid[!is.na(grid$n), ], aes(fill = band)) +
    scale_fill_manual(name = "Required samples (n)", limits = labs_b, drop = FALSE,
                      values = setNames(viridisLite::viridis(length(labs_b), direction = -1), labs_b)) +
    geom_contour(aes(z = n), breaks = br[is.finite(br) & br > 0],
                 colour = "white", linewidth = 0.3, alpha = 0.7, na.rm = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = range(grid$cv_site) + c(0, 0.02), ylim = range(grid$width)) +
    labs(x = expression(bold("Site CV  (" * sigma / mu * ")")),
         y = "Cross-section width (m)", subtitle = title)
  g <- g + (if (is.null(base_theme)) theme_bw() else base_theme) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(face = "bold"))
  if (nrow(na_part) > 0) {
    g <- g + geom_raster(data = na_part, fill = "grey85") +
      annotate("text", x = max(grid$cv_site), y = max(grid$width),
               label = paste0("unachievable within ", attr(grid, "n_max"), " samples"),
               hjust = 1.05, vjust = 1.8, size = 3, colour = "grey30")
  }
  if (!is.null(others)) {        # labelled points, e.g. the study sites
    g <- g + geom_point(data = others, aes(x = cv_site, y = width), shape = 21, size = 2.2,
                        fill = "white", colour = "black", stroke = 0.6, alpha = 0.85, inherit.aes = FALSE) +
      geom_text(data = others, aes(x = cv_site, y = width, label = label), colour = "white",
                size = 2.8, vjust = -0.9, inherit.aes = FALSE)
  }
  if (!is.null(mark)) {
    g <- g + annotate("point", x = mark$cv_site, y = mark$width, shape = 21,
                      size = 3.5, fill = "white", colour = "black", stroke = 1.2)
  }
  g
}
