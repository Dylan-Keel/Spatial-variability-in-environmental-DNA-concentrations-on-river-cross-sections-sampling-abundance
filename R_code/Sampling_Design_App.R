## Sampling design for eDNA on a river cross section
## An interactive version of Figures 5, 7, and 8, with the six Klamath cross
## sections of the study available as worked examples.
## Run from this folder with: shiny::runApp("Sampling_Design_App.R")

library(shiny)
library(ggplot2)
source("Sampling_Design_Functions.R", local = TRUE)

## ---- Klamath 2021 cross sections ---------------------------------------------
SITE_CODES <- c("KLM", "I5", "TOH", "BVR", "KMN", "SV")
SITE_W     <- c(45, 37, 37, 39, 34, 42.5)                 # wetted width (m), as in Figure 5
SPECIES    <- list(cs = list(label = "C. shasta",      p = 0.60, range = 50),
                   ot = list(label = "O. tshawytscha", p = 0.33, range = 14.4))
csd <- read.csv("../data/CSD.csv")
vol.cor <- function(c.well) { x <- round(c.well * csd$Elution.volum..uL. / 15 * csd$Volume..ml. / 1000, 0); x[is.na(x)] <- 0; x }
SAMPLES <- data.frame(site = match(csd$Site, SITE_CODES), x = csd$Distance.Left.Bank.m, depth = csd$D.M.,
                      cs = vol.cor(csd$C.shasta.c.well), ot = vol.cor(csd$Chinook.c.well))

## bed profile of a site: deepest sample on each vertical, scaled up a little, interpolated across the width
site_channel <- function(s, W) {
  bed <- tapply(s$depth, s$x, max)
  bed <- bed[bed > 0]
  xl  <- min(0, min(s$x) - 1); xr <- max(W, max(s$x) + 1)
  kx  <- c(xl, as.numeric(names(bed)), xr); kd <- c(0, bed * 1.25 + 0.1, 0)
  list(xl = xl, xr = xr, max_depth = max(kd), depth_fn = function(x) approx(kx, kd, xout = x, rule = 2)$y)
}
site_context <- function(id, sp) {
  s <- SAMPLES[SAMPLES$site == id, ]
  v <- s[[sp]]
  list(id = id, code = SITE_CODES[id], sp = sp, samples = s, vals = v, n = length(v), W = SITE_W[id],
       mean = mean(v), cv = sd(v) / mean(v), p = SPECIES[[sp]]$p, range = SPECIES[[sp]]$range,
       channel = site_channel(s, SITE_W[id]))
}
## the slider values a site loads; departing from any of them returns the page to a simulated section
site_inputs <- function(ctx) list(mean_conc = round(ctx$mean), cv_site = round(ctx$cv, 2), p = ctx$p,
                                  range_m = ctx$range, width_m = ctx$W, n_given = ctx$n)

## ---- simulated field ----------------------------------------------------------
## a correlated log-normal field in a channel whose bed is depth_fn(x), for the picture
matern_field <- function(xl, xr, depth_fn, max_depth, n_x = 35, n_y = 14, mean_conc, cv_site, p, range_m, seed = 42) {
  set.seed(seed)
  grid <- expand.grid(x = seq(xl, xr, length.out = n_x), y = seq(0, max_depth, length.out = n_y))
  grid$in_channel <- grid$y <= depth_fn(grid$x)
  sig2_log <- log(1 + cv_site^2)
  mu_log   <- log(mean_conc) - sig2_log / 2
  D     <- as.matrix(dist(grid[, c("x", "y")]))
  Cor   <- matern_cor(D, sqrt(8) / range_m)
  Sigma <- p * sig2_log * Cor + diag((1 - p) * sig2_log + 1e-6, nrow(grid))
  L <- tryCatch(t(chol(Sigma)), error = function(e) t(chol(Sigma + diag(1e-4, nrow(Sigma)))))
  z <- mu_log + as.vector(L %*% rnorm(nrow(grid)))
  grid$conc <- ifelse(grid$in_channel, exp(z), NA_real_)
  grid
}
## Gaussian field with mean mu and variance s2 on the channel grid (for the conditional simulation)
gauss_field <- function(xl, xr, depth_fn, max_depth, n_x = 35, n_y = 14, mu, s2, p, range_m, seed = 42) {
  set.seed(seed)
  grid <- expand.grid(x = seq(xl, xr, length.out = n_x), y = seq(0, max_depth, length.out = n_y))
  grid$in_channel <- grid$y <= depth_fn(grid$x)
  D     <- as.matrix(dist(grid[, c("x", "y")]))
  Sigma <- p * s2 * matern_cor(D, sqrt(8) / range_m) + diag((1 - p) * s2 + 1e-6 * s2, nrow(grid))
  L <- tryCatch(t(chol(Sigma)), error = function(e) t(chol(Sigma + diag(1e-4 * s2, nrow(Sigma)))))
  grid$g <- mu + as.vector(L %*% rnorm(nrow(grid)))
  grid
}
## Klamath site: kriging smoother of the measured concentrations, with the species' p and
## range, plus a simulated residual with the kriging variance (a conditional simulation:
## unconditional field minus its own kriging), truncated at zero. Re-simulating redraws the
## residual only.
conditional_field <- function(ctx, p, range_m, seed = 42) {
  s  <- ctx$samples; z <- ctx$vals; mu <- mean(z); s2 <- max(var(z), 1e-6); n <- length(z)
  ch <- ctx$channel
  grid <- gauss_field(ch$xl, ch$xr, ch$depth_fn, ch$max_depth, mu = mu, s2 = s2, p = p, range_m = range_m, seed = seed)
  kappa <- sqrt(8) / range_m
  obs <- cbind(s$x, s$depth)
  Sigma <- s2 * (p * matern_cor(as.matrix(dist(obs)), kappa) + (1 - p) * diag(n))
  ## value of the unconditional field at each sample: nearest grid cell inside the channel
  inside <- grid[grid$in_channel, ]
  near <- apply(obs, 1, function(o) which.min((inside$x - o[1])^2 + (inside$y - o[2])^2))
  s_obs <- inside$g[near]
  w1 <- solve(Sigma, z - mu); w2 <- solve(Sigma, s_obs - mu)
  K  <- s2 * p * matern_cor(sqrt(outer(grid$x, obs[, 1], "-")^2 + outer(grid$y, obs[, 2], "-")^2), kappa)
  pred <- mu + as.vector(K %*% w1); sim_k <- mu + as.vector(K %*% w2)
  grid$smooth <- ifelse(grid$in_channel, pmax(0, pred), NA_real_)
  grid$conc   <- ifelse(grid$in_channel, pmax(0, grid$g + pred - sim_k), NA_real_)
  grid
}
## colour scales: blues for simulated sections and O. tshawytscha, greens for C. shasta
PALETTE <- list(blue = c("#c6e2f5", "#08306b"), green = c("#e5f5e0", "#00441b"))
parabola <- function(W, max_depth = 3) function(x) max_depth * 4 * (x / W) * (1 - x / W)
sample_positions <- function(n, width, max_depth, strategy, seed = 123) {
  set.seed(seed)
  x_s <- if (strategy == "random") runif(n, 0, width) else (seq_len(n) - 0.5) * width / n
  x_n <- x_s / width
  data.frame(x = x_s, y = runif(n) * max_depth * 4 * x_n * (1 - x_n))
}
fmt_n   <- function(x) if (is.na(x) || is.infinite(x)) "not reachable" else as.character(ceiling(x))
fmt_num <- function(x) format(round(x), big.mark = ",")

hint <- function(...) tags$p(class = "hint", ...)
HINTS <- list(
  site      = "Loads the measured width, mean, site CV, and sample positions of that section, with the species' spatial structure from the fitted models.",
  mean_conc = "Mean eDNA concentration on the section. Only scales the colours of the simulated field; a Klamath site uses its measured samples instead.",
  cv_site   = "Coefficient of variation (SD ÷ mean) among single 1-L samples taken across the section, as measured in a pilot study. On the Klamath: 0.40–1.43.",
  p         = HTML("Share of the variance among samples that is spatially structured, 1 &minus; nugget/sill from a semivariogram. Klamath estimates: 0.60 (<i>C. shasta</i>), 0.33 (<i>O. tshawytscha</i>)."),
  range_m   = HTML("Distance at which the Mat&eacute;rn (&nu; = 1) correlation falls to about 0.14. Klamath estimates: about 50 m (<i>C. shasta</i>), about 14 m (<i>O. tshawytscha</i>)."),
  width_m   = "Wetted width of the cross-section. Samples spread over a wider section, relative to the correlation range, are less correlated.",
  strategy  = "Simple random: positions drawn at random across the width. Spatially stratified: one sample at the centre of each of n equal-width strata. As sampled: the verticals and depths used in 2021.",
  cv_target = "Precision wanted for the estimated section mean: its standard error ÷ the mean. A target CV is a precision criterion, not statistical power. Rules of thumb: 0.1 research, 0.2 management, 0.3 pilot work.",
  n_given   = "Number of 1-L samples collected on the section."
)
STRATEGY_SIM  <- c("Simple random" = "random", "Spatially stratified" = "stratified")
STRATEGY_SITE <- c(STRATEGY_SIM, "As sampled at the site (2021)" = "sampled")

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', sans-serif; }
    details.howto { background: #fffdf3; border: 1px solid #e8dfb8; border-radius: 4px; padding: 6px 14px; margin: 6px 0 10px 0; font-size: 13px; }
    details.howto summary { cursor: pointer; font-weight: bold; color: #7d6608; }
    details.howto ol { margin: 6px 0 4px 18px; line-height: 1.5; }
    .well h4 { color: #1a5276; border-bottom: 1px solid #d5e8f3; padding-bottom: 4px; margin-top: 18px; }
    .result-box { background: #eaf3fb; border-left: 4px solid #2980b9; border-radius: 4px; padding: 10px 14px; margin: 8px 0; font-size: 14px; }
    .result-box b { color: #1a5276; }
    .banner { display: inline-block; font-size: 13px; font-weight: bold; padding: 3px 12px; border-radius: 12px; margin: 8px 0 4px 0; }
    .banner.sim { background: #e5e8ea; color: #555; } .banner.site { background: #1a5276; color: white; }
    .panel-title { color: #1a5276; font-weight: bold; font-size: 15px; margin-bottom: 6px; }
    .interp { background: white; border: 1px solid #d5e8f3; border-left: 4px solid #1a5276; border-radius: 4px; padding: 10px 14px; margin: 6px 0 14px 0; font-size: 13.5px; line-height: 1.55; }
    .interp h3 { color: #1a5276; font-size: 15px; margin: 0 0 6px 0; } .interp ul { margin: 2px 0 4px 18px; }
    .interp .ok { color: #1e8449; font-weight: bold; } .interp .short { color: #c0392b; font-weight: bold; }
    .hint { font-size: 11.5px; color: #666; line-height: 1.35; margin: -8px 0 12px 0; }
    .cite { font-size: 12.5px; color: #555; margin: -6px 0 8px 0; }
    details.defs { background: #f7fafd; border: 1px solid #d5e8f3; border-radius: 4px; padding: 8px 14px; margin: 10px 0; font-size: 13px; }
    details.defs summary { cursor: pointer; font-weight: bold; color: #1a5276; }
    details.defs dt { font-weight: bold; margin-top: 6px; }
    details.defs dd { margin-left: 0; }
  ")),
  tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function () {
      var h = document.getElementById('howto'); if (!h) return;
      try { if (localStorage.getItem('edna_howto_closed') === '1') h.removeAttribute('open'); } catch (e) {}
      h.addEventListener('toggle', function () { try { localStorage.setItem('edna_howto_closed', h.open ? '0' : '1'); } catch (e) {} });
    });"))),
  titlePanel("eDNA River Cross-Section Sampling Design"),
  tags$p(class = "cite", HTML("Companion to Keel, D. J., N. A. Som, M. J. Henderson, and A. P. Kinziger. In prep. <i>Spatial variability in environmental DNA concentrations on river cross sections: implications for sampling design and abundance estimation</i>.")),
  tags$details(class = "howto", id = "howto", open = NA,
    tags$summary("How to use this page"),
    tags$ol(
      tags$li(HTML("<b>Choose a cross section.</b> <i>Simulated</i>: set the site CV, spatial structure, and width yourself, from a pilot study or from the Klamath values given in the hints. <i>Klamath site</i>: load one of the six sections of Keel et al. (2026) and read what the 2021 sampling achieved for either species.")),
      tags$li(HTML("<b>Set the precision you need.</b> The target CV is the standard error of the estimated section mean divided by the mean. The right-hand column gives the number of 1-L samples for that target, and the CV reached with a given n.")),
      tags$li(HTML("<b>Choose placement.</b> Simple random uses the manuscript's formula, n = (CV<sub>site</sub> / CV<sub>target</sub>)&sup2;. Spatially stratified spreads samples evenly across the width and needs fewer samples when concentrations are autocorrelated. For a Klamath site, <i>as sampled</i> uses the positions actually sampled.")),
      tags$li(HTML("<b>Plan in general</b> with the contour maps at the bottom, which give n against site CV and width for any section.")),
      tags$li("Changing the site CV, spatial structure, width, mean, n, or placement returns the page to a simulated section; the target CV and the contour controls do not.")
    )),
  sidebarLayout(
    sidebarPanel(width = 3,
      h4("Cross section"),
      radioButtons("mode", NULL, choices = c("Simulated" = "sim", "Klamath site from the paper" = "site"), selected = "sim"),
      conditionalPanel("input.mode == 'site'",
        selectInput("site_sel", "Site", choices = c("1 (KLM, upstream)" = 1, "2 (I5)" = 2, "3 (TOH)" = 3, "4 (BVR)" = 4, "5 (KMN)" = 5, "6 (SV, downstream)" = 6), selected = 3),
        selectInput("sp_sel", "Species", choices = c("C. shasta" = "cs", "O. tshawytscha" = "ot")),
        hint(HINTS$site)),
      h4("Concentration"),
      numericInput("mean_conc", "Mean (copies / L)", value = 1000, min = 1, step = 100), hint(HINTS$mean_conc),
      sliderInput("cv_site", "Site CV (single-sample CV)", min = 0.10, max = 2.00, value = 0.80, step = 0.01, ticks = FALSE), hint(HINTS$cv_site),
      h4("Spatial structure"),
      sliderInput("p", HTML("Spatial proportion of variance, <i>p</i>"), min = 0, max = 0.95, value = 0.50, step = 0.01, ticks = FALSE), hint(HINTS$p),
      sliderInput("range_m", "Correlation range (m)", min = 1, max = 100, value = 30, step = 0.1, ticks = FALSE), hint(HINTS$range_m),
      sliderInput("width_m", "Cross-section width (m)", min = 10, max = 150, value = 40, step = 0.5, ticks = FALSE), hint(HINTS$width_m),
      uiOutput("icc_box"),
      h4("Sample placement"),
      radioButtons("strategy", NULL, choices = STRATEGY_SIM, selected = "stratified"),
      hint(HINTS$strategy),
      hr(),
      h4("Solve for n"),
      sliderInput("cv_target", "Target CV (precision goal)", min = 0.05, max = 1.00, value = 0.20, step = 0.01, ticks = FALSE), hint(HINTS$cv_target),
      h4("Solve for CV"),
      sliderInput("n_given", "Number of samples (n)", min = 1, max = 60, value = 9, step = 1, ticks = FALSE), hint(HINTS$n_given),
      hr(),
      actionButton("resim", "↻  Re-simulate field", class = "btn-primary btn-block", width = "100%")
    ),
    mainPanel(width = 9,
      uiOutput("banner"),
      fluidRow(
        column(8, uiOutput("cs_title"), plotOutput("cross_section", height = "240px")),
        column(4, br(), div(class = "panel-title", "Samples needed for the target CV"), uiOutput("n_required"),
               br(), div(class = "panel-title", "Precision (CV of the estimated mean) with the given n"), uiOutput("cv_achieved"),
               HTML('<div class="result-box" style="border-left-color:#7f8c8d; background:#f4f6f7">
  <b>Which mean are you estimating?</b><br>
  <b>Section mean</b> (n<sub>ind</sub>, n<sub>strat</sub>): the mean concentration on this cross section on the day it is sampled. Variation that every sample on the section shares is part of what you are measuring, so adding samples keeps improving precision.<br>
  <b>Reach mean</b> (n<sub>adj</sub>): using this one section to stand for a longer reach of river. The section-wide deviation is then error that more samples on the same section cannot remove, so precision levels off at CV<sub>site</sub>&radic;(&rho;/(1&minus;&rho;)). Only sampling more sections reduces it.
</div>'))
      ),
      uiOutput("site_interp"),
      hr(),
      div(class = "panel-title", "Precision of the section mean against sample size"),
      plotOutput("cv_curve", height = "280px"),
      hr(),
      div(class = "panel-title", "Samples required against site CV and cross-section width"),
      fluidRow(
        column(4, radioButtons("contour_layout", NULL,
                               choices = c("2 panels: random vs. stratified at the chosen p" = "two",
                                           "4 panels: random, then p = 0.25, 0.50, 0.75" = "four"), selected = "two")),
        column(4, sliderInput("n_max_contour", "Cap on n (grey = more than this)", min = 50, max = 500, value = 300, step = 50, ticks = FALSE)),
        column(4, uiOutput("contour_settings"))
      ),
      plotOutput("contour_plot", height = "480px"),
      hint("Each panel maps the number of 1-L samples needed to reach the target CV for the section mean, against site CV (x) and cross-section width (y), at the correlation range set in the sidebar. The white point marks the current site CV and width; in the site view the six Klamath sites are shown for the chosen species."),
      tags$details(class = "defs", tags$summary("Definitions and formulas"),
            tags$dl(
              tags$dt("Site CV"), tags$dd("Standard deviation of single 1-L sample concentrations on the section divided by their mean. This is what a pilot study measures."),
              tags$dt("Target CV"), tags$dd("Standard error of the estimated section mean divided by the mean."),
              tags$dt("Simple random placement"),
              tags$dd(HTML("CV(n) = CV<sub>site</sub> / &radic;n, so n = (CV<sub>site</sub> / CV<sub>target</sub>)&sup2; (Cochran 1977). This is the manuscript's formula.")),
              tags$dt("Spatially stratified placement"),
              tags$dd(HTML("Samples at the centres of n equal strata across the width. Evenly spaced samples converge on the section average, so precision improves faster than 1/&radic;n when the field is autocorrelated. The gain depends on <i>p</i> and on the correlation range relative to the width; it disappears when the range is much larger or much smaller than the section.")),
              tags$dt("As sampled"),
              tags$dd("The same calculation with the positions actually sampled across the width. The precision model is one-dimensional, so paired surface and depth samples on one vertical count as one position."),
              tags$dt(HTML("Intraclass correlation &rho;")),
              tags$dd(HTML("&rho; = <i>p</i> &times; <i>c&#772;</i>, where <i>c&#772;</i> is the mean Mat&eacute;rn correlation between two random positions on the section. It is shown in the sidebar for the current width and range.")),
              tags$dt("Reach mean (design effect)"),
              tags$dd(HTML("If the target is the mean of the reach rather than of the sampled section, the deviation shared by all samples on the section is error: CV(n) = CV<sub>site</sub> &radic;[(1 + (n &minus; 1)&rho;) / (n(1 &minus; &rho;))], which levels off at CV<sub>site</sub> &radic;(&rho; / (1 &minus; &rho;)). Only more sections reduce it.")),
              tags$dt("Simulated field"),
              tags$dd(HTML("A log-normal field with the chosen mean and site CV, a share <i>p</i> of its log-variance spatially correlated with the chosen range. The simulated channel is 3 m deep at the thalweg.")),
              tags$dt("Klamath site field"),
              tags$dd(HTML("A kriging smoother of log(concentration + 1) through the measured samples, with the species' <i>p</i> and range, plus a simulated residual drawn with the kriging variance (a conditional simulation). It honours the samples where they were taken, within the nugget, and has realistic error between them; re-simulate redraws only the error. The channel follows the sampled verticals.")),
              tags$dt("Klamath sites"),
              tags$dd(HTML("The six cross sections of Keel et al. (2026), sampled in 2021 and numbered from upstream to downstream. Site CV, mean, width, and sample positions are the measured values; <i>p</i> and the range are the species-level estimates from the fitted spatial models (<i>C. shasta</i>: p = 0.60, range about 50 m; <i>O. tshawytscha</i>: p = 0.33, range about 14 m)."))
            ))
    )
  )
)

server <- function(input, output, session) {
  tab <- reactive(matern_tables(input$range_m))
  sim_seed <- reactiveVal(42)
  observeEvent(input$resim, sim_seed(sample.int(1e6, 1)))

  ## ---- simulated section or a Klamath site ----
  mode <- reactiveVal("sim")
  site <- reactiveVal(NULL)
  enter_site <- function() {
    ctx <- site_context(as.integer(input$site_sel), input$sp_sel)
    site(ctx); mode("site"); sim_seed(1000 + ctx$id)
    v <- site_inputs(ctx)
    updateNumericInput(session, "mean_conc", value = v$mean_conc)
    for (nm in c("cv_site", "p", "range_m", "width_m", "n_given")) updateSliderInput(session, nm, value = v[[nm]])
    updateRadioButtons(session, "strategy", choices = STRATEGY_SITE, selected = "sampled")
  }
  leave_site <- function() {
    keep <- if (isTRUE(input$strategy == "sampled")) "stratified" else input$strategy
    mode("sim"); site(NULL)
    updateRadioButtons(session, "mode", selected = "sim")
    updateRadioButtons(session, "strategy", choices = STRATEGY_SIM, selected = keep)
  }
  observeEvent(input$mode, { if (input$mode == "site") enter_site() else if (mode() == "site") leave_site() }, ignoreInit = TRUE)
  observeEvent(c(input$site_sel, input$sp_sel), { if (input$mode == "site") enter_site() }, ignoreInit = TRUE)
  ## any departure from the loaded site values returns to a simulated section
  for (nm in c("mean_conc", "cv_site", "p", "range_m", "width_m", "n_given")) local({
    id <- nm
    observeEvent(input[[id]], {
      ctx <- site()
      if (mode() == "site" && !is.null(ctx) && !isTRUE(all.equal(as.numeric(input[[id]]), site_inputs(ctx)[[id]], tolerance = 1e-6))) leave_site()
    }, ignoreInit = TRUE)
  })
  observeEvent(input$strategy, { if (mode() == "site" && input$strategy != "sampled") leave_site() }, ignoreInit = TRUE)
  in_site <- reactive(mode() == "site" && !is.null(site()))

  output$banner <- renderUI({
    if (in_site()) { s <- site(); HTML(sprintf('<div class="banner site">Klamath site %d (%s) &middot; <i>%s</i> &middot; sampled 2021</div>', s$id, s$code, SPECIES[[s$sp]]$label)) }
    else HTML('<div class="banner sim">Simulated cross section</div>')
  })
  output$cs_title <- renderUI({
    if (in_site()) div(class = "panel-title", sprintf("Site %d: concentration smoothed from the %d measured samples, with simulated error between them", site()$id, site()$n))
    else div(class = "panel-title", "Simulated river cross-section")
  })
  output$icc_box <- renderUI({
    rho <- icc_section(input$width_m, input$p, tab())
    HTML(sprintf('<div class="result-box" style="font-size:13px"><b>Implied intraclass correlation</b> for this width and range:<br>&rho; = p &times; c&#772; = %.2f</div>', rho))
  })

  output$cross_section <- renderPlot({
    if (in_site()) { s <- site(); ch <- s$channel } else { ch <- list(xl = 0, xr = input$width_m, max_depth = 3, depth_fn = parabola(input$width_m)) }
    gd <- if (in_site()) conditional_field(s, input$p, input$range_m, seed = sim_seed())
          else matern_field(ch$xl, ch$xr, ch$depth_fn, ch$max_depth, mean_conc = input$mean_conc, cv_site = input$cv_site,
                            p = input$p, range_m = input$range_m, seed = sim_seed())
    pal <- if (in_site() && s$sp == "cs") PALETTE$green else PALETTE$blue
    outline <- data.frame(x = seq(ch$xl, ch$xr, length.out = 300)); outline$y <- -ch$depth_fn(outline$x)
    g <- ggplot(gd[!is.na(gd$conc), ], aes(x = x, y = -y, fill = conc)) +
      geom_tile() +
      geom_line(data = outline, aes(x = x, y = y), inherit.aes = FALSE, colour = "grey20", linewidth = 1.2) +
      annotate("segment", x = ch$xl, xend = ch$xr, y = 0, yend = 0, colour = "grey20", linewidth = 1.2)
    if (in_site()) {
      pts <- data.frame(x = s$samples$x, y = -s$samples$depth, conc = s$vals)
      g <- g + geom_point(data = pts, aes(x = x, y = y, fill = conc), shape = 21, size = 3.5, colour = "black", stroke = 1.2)
    } else {
      sp <- sample_positions(input$n_given, input$width_m, 3, input$strategy)
      g <- g + geom_point(data = sp, aes(x = x, y = -y), inherit.aes = FALSE, shape = 21, size = 3.5, fill = "white", colour = "black", stroke = 1.3)
    }
    lim <- if (in_site()) c(0, max(s$vals)) else NULL
    g + scale_fill_gradient(low = pal[1], high = pal[2], name = "Conc.\n(copies/L)", labels = scales::comma, limits = lim, oob = scales::squish) +
      scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x, " m")) +
      scale_y_continuous(expand = c(0.02, 0), labels = abs) +
      labs(x = "Distance from left bank", y = "Depth (m)") +
      theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"))
  })

  output$n_required <- renderUI({
    n_r <- n_random(input$cv_site, input$cv_target)
    n_s <- n_stratified(input$cv_site, input$cv_target, input$width_m, input$p, tab(), n_max = 2000)
    n_h <- n_reach(input$cv_site, input$cv_target, input$width_m, input$p, tab(), n_max = 2000)
    HTML(paste0('<div class="result-box">',
                '<b>Simple random, n<sub>ind</sub>:</b> ', fmt_n(n_r), '<br>',
                '<b>Spatially stratified, n<sub>strat</sub>:</b> ', fmt_n(n_s), '<br>',
                '<small style="color:#555">Reach mean, autocorrelation-adjusted n<sub>adj</sub>: ', fmt_n(n_h), '</small></div>'))
  })

  output$cv_achieved <- renderUI({
    n <- input$n_given
    sampled <- if (in_site()) sprintf('<b>As sampled (n = %d):</b> %.3f<br>', site()$n, cv_positions(input$cv_site, site()$samples$x, input$width_m, input$p, tab())) else ""
    HTML(paste0('<div class="result-box">', sampled,
                '<b>Simple random:</b> ', round(cv_random(input$cv_site, n), 3), '<br>',
                '<b>Spatially stratified:</b> ', round(cv_stratified(input$cv_site, n, input$width_m, input$p, tab()), 3), '<br>',
                '<small style="color:#555">Reach mean, design effect: ', round(cv_reach(input$cv_site, n, input$width_m, input$p, tab()), 3), '</small></div>'))
  })

  ## ---- what the 2021 sampling achieved at the chosen site ----
  output$site_interp <- renderUI({
    if (!in_site()) return(NULL)
    s <- site(); sp <- SPECIES[[s$sp]]; n <- s$n; cvS <- input$cv_site; W <- input$width_m
    verts <- sort(unique(s$samples$x)); spacing <- median(diff(verts))
    n_surf <- sum(s$samples$depth == 0); n_deep <- n - n_surf
    cv_ind <- cv_random(cvS, n); cv_as <- cv_positions(cvS, s$samples$x, W, input$p, tab())
    n_ind <- n_random(cvS, input$cv_target); n_str <- n_stratified(cvS, input$cv_target, W, input$p, tab(), n_max = 2000)
    zeros <- sum(s$vals == 0); top3 <- sum(sort(s$vals, decreasing = TRUE)[1:3]) / sum(s$vals)
    se_cv <- cvS / sqrt(2 * n) * sqrt(1 + 2 * cvS^2)
    verdict <- if (cv_as <= input$cv_target)
      sprintf('<span class="ok">This meets the target CV of %.2f.</span> The smallest design that would also meet it is about n<sub>strat</sub> = %s evenly spaced samples (n<sub>ind</sub> = %d if treated as independent).', input$cv_target, fmt_n(n_str), n_ind)
    else sprintf('<span class="short">This misses the target CV of %.2f.</span> Reaching it would take about n<sub>strat</sub> = %s evenly spaced samples (n<sub>ind</sub> = %d if treated as independent)%s.', input$cv_target, fmt_n(n_str), n_ind,
                 if (is.na(n_str)) "" else if (n_str > n) sprintf(", %d more than were taken", n_str - n) else ", fewer than were taken, once the layout is spread evenly")
    species_note <- if (s$sp == "ot")
      "In the fitted model <i>O. tshawytscha</i> concentration fell with depth, so a surface-only design would overestimate the section mean. Mixing surface and depth samples, as done here, matters more for this species than spreading them across the width."
    else "No covariate explained <i>C. shasta</i> concentration; its variation was spatial autocorrelation (p = 0.60, range about 50 m) and noise. Spreading samples across the width is the main lever for precision."
    HTML(paste0('<div class="interp">',
      sprintf('<h3>Site %d (%s), <i>%s</i>: what the 2021 sampling of Keel et al. (2026) achieved</h3>', s$id, s$code, sp$label),
      sprintf('<p><b>Design.</b> %d one-liter samples on %d verticals about %s m apart across a %s m wide section: %d surface grabs and %d samples at depth. Measured mean %s copies/L, site CV %.2f%s. The three highest samples hold %d%% of the summed concentration.</p>',
              n, length(verts), format(round(spacing, 1)), format(s$W), n_surf, n_deep, fmt_num(s$mean), cvS,
              if (zeros > 0) sprintf("; %d of the %d samples had no detectable <i>%s</i>", zeros, n, sp$label) else "", round(top3 * 100)),
      sprintf('<p><b>Precision reached.</b> Treating the samples as independent, the CV of the estimated section mean is %.3f (= %.2f / &radic;%d). Using the positions actually sampled and the species\' autocorrelation (p = %.2f, range %s m), it is %.3f. %s</p>',
              cv_ind, cvS, n, input$p, format(round(input$range_m, 1)), cv_as, verdict),
      sprintf('<p><b>What this section says.</b> %s</p>', species_note),
      '<p><b>Limits of this reading.</b></p><ul>',
      sprintf('<li>The site CV is estimated from these %d samples and is itself uncertain (standard error roughly &plusmn;%.2f); the required n scales with its square.</li>', n, se_cv),
      '<li><i>p</i> and the range are species-level estimates pooled across the six sites, not estimates for this section.</li>',
      '<li>The samples come from one day. The numbers describe the section mean on that day, not day-to-day change or the reach mean (see n<sub>adj</sub>).</li>',
      '<li>The precision model is one-dimensional across the width, so paired surface and depth samples on a vertical count as one position; the as-sampled CV is conservative and the independence value is the optimistic bound.</li>',
      '<li>The site CV includes ddPCR measurement error and the volume filtered, which differed among samples.</li>',
      '</ul></div>'))
  })

  output$cv_curve <- renderPlot({
    n_seq <- 1:60
    df <- rbind(
      data.frame(n = n_seq, CV = cv_random(input$cv_site, n_seq), type = "Simple random (independence)"),
      data.frame(n = n_seq, CV = cv_stratified(input$cv_site, n_seq, input$width_m, input$p, tab()), type = "Spatially stratified"),
      data.frame(n = n_seq, CV = cv_reach(input$cv_site, n_seq, input$width_m, input$p, tab()), type = "Reach mean (autocorrelation-adjusted)"))
    df$type <- factor(df$type, levels = c("Simple random (independence)", "Spatially stratified", "Reach mean (autocorrelation-adjusted)"))
    y_max <- min(input$cv_site * 1.3, 2.5)
    g <- ggplot(df, aes(x = n, y = CV, colour = type, linetype = type)) +
      geom_hline(yintercept = input$cv_target, linetype = "dashed", colour = "grey30") +
      annotate("text", x = 59, y = input$cv_target + y_max * 0.03, label = sprintf("Target CV = %.2f", input$cv_target), hjust = 1, size = 3.5, colour = "grey25") +
      geom_vline(xintercept = input$n_given, linetype = "dotted", colour = "grey55") +
      geom_line(linewidth = 0.95)
    if (in_site()) {
      pt <- data.frame(n = site()$n, CV = cv_positions(input$cv_site, site()$samples$x, input$width_m, input$p, tab()))
      g <- g + geom_point(data = pt, aes(x = n, y = CV), inherit.aes = FALSE, shape = 23, size = 4, fill = "#1a5276", colour = "white") +
        annotate("text", x = pt$n + 1, y = pt$CV + y_max * 0.06, label = sprintf("as sampled in 2021 (n = %d)", pt$n), hjust = 0, size = 3.5, colour = "#1a5276")
    }
    g + scale_colour_manual(values = c("#0072B2", "#009E73", "#D55E00"), name = NULL) +
      scale_linetype_manual(values = c("solid", "22", "dotted"), name = NULL) +
      scale_x_continuous(limits = c(1, 60), breaks = seq(0, 60, 10), expand = expansion(mult = c(0, 0.01))) +
      coord_cartesian(ylim = c(0, y_max)) +
      labs(x = "Number of 1-L samples (n)", y = "CV of the estimated mean") +
      theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"), legend.position = "bottom")
  })

  contour_grids <- reactive({
    cv_seq <- seq(0.10, 1.50, by = 0.05); w_seq <- seq(5, 150, by = 5)
    mk <- function(p, design) n_required_grid(cv_seq, w_seq, input$cv_target, p, tab(), design, n_max = input$n_max_contour)
    if (input$contour_layout == "two") {
      list(list(grid = mk(0, "random"), title = "A)  Simple random placement (any width)"),
           list(grid = mk(input$p, "stratified"), title = sprintf("B)  Spatially stratified, p = %.2f", input$p)))
    } else {
      c(list(list(grid = mk(0, "random"), title = "A)  Simple random placement (any width)")),
        Map(function(p, L) list(grid = mk(p, "stratified"), title = sprintf("%s)  Spatially stratified, p = %.2f", L, p)),
            c(0.25, 0.5, 0.75), c("B", "C", "D")))
    }
  })
  output$contour_settings <- renderUI({
    HTML(sprintf('<div class="result-box" style="font-size:13px"><b>Target CV:</b> %.2f<br><b>Range:</b> %s m<br><b>Marker:</b> site CV = %.2f, width = %s m</div>',
                 input$cv_target, format(round(input$range_m, 1)), input$cv_site, format(round(input$width_m, 1))))
  })
  output$contour_plot <- renderPlot({
    gl <- contour_grids()
    mark <- list(cv_site = input$cv_site, width = input$width_m)
    others <- if (in_site()) { sp <- site()$sp
      do.call(rbind, lapply(1:6, function(i) { v <- SAMPLES[[sp]][SAMPLES$site == i]; data.frame(cv_site = sd(v) / mean(v), width = SITE_W[i], label = i) })) } else NULL
    panels <- lapply(gl, function(g) contour_panel(g$grid, g$title, mark = mark, others = others))
    ggpubr::ggarrange(plotlist = panels, ncol = 2, nrow = if (length(panels) == 2) 1 else 2, common.legend = TRUE, legend = "right")
  })
}

shinyApp(ui, server)
