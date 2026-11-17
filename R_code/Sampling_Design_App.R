## Sampling design tool for eDNA on a river cross section
## A step-by-step interactive version of Figures 5 and 7: how many 1-L
## samples are needed to estimate the mean eDNA concentration of a cross section.
## Run from this folder with: shiny::runApp("Sampling_Design_App.R")

library(shiny)
library(ggplot2)
source("Sampling_Design_Functions.R", local = TRUE)

N_MAX <- 500

## irregular channel bed: an off-centre thalweg with bars and pools, zero depth at both banks; u = x / width
bed_profile <- function(seed, max_depth = 3) {
  set.seed(seed + 7)
  skew <- runif(1, 0.6, 1.6); a <- c(runif(1, 0.15, 0.35), runif(1, 0.08, 0.18)); k <- c(sample(1:2, 1), sample(3:5, 1)); f <- runif(2, 0, 2 * pi)
  raw <- function(u) sin(pi * u^skew) * pmax(0.25, 1 + a[1] * sin(2 * pi * k[1] * u + f[1]) + a[2] * sin(2 * pi * k[2] * u + f[2]))
  top <- max(raw(seq(0, 1, length.out = 201)))
  function(u) max_depth * pmax(0, raw(pmin(1, pmax(0, u)))) / top
}
## a correlated log-normal field in a stylised channel, relative to the section mean (for the picture)
matern_field <- function(width, max_depth = 3, n_x = 40, n_y = 14, cv_site, p, range_m, seed = 42) {
  grid <- expand.grid(x = seq(0, width, length.out = n_x), y = seq(0, max_depth, length.out = n_y))
  grid$in_channel <- grid$y <= bed_profile(seed, max_depth)(grid$x / width)
  set.seed(seed)
  sig2  <- log(1 + cv_site^2)
  D     <- as.matrix(dist(grid[, c("x", "y")]))
  Sigma <- p * sig2 * matern_cor(D, sqrt(8) / range_m) + diag((1 - p) * sig2 + 1e-6, nrow(grid))
  L <- tryCatch(t(chol(Sigma)), error = function(e) t(chol(Sigma + diag(1e-4, nrow(Sigma)))))
  grid$rel <- ifelse(grid$in_channel, exp(-sig2 / 2 + as.vector(L %*% rnorm(nrow(grid)))), NA_real_)
  grid
}
fmt_n <- function(x) if (is.na(x) || x > N_MAX) paste("more than", N_MAX) else as.character(x)
CVT <- "CV<sub>target</sub>"; CVS <- "CV<sub>site</sub>"

## a parameter label with its technical definition on hover
tip <- function(label, value_id, definition)
  HTML(sprintf('%s <span class="v">%s</span><span class="tip" tabindex="0">i<span class="tt">%s</span></span>',
               label, as.character(textOutput(value_id, inline = TRUE)), definition))
lay  <- function(...) tags$p(class = "lay", HTML(paste0(...)))
klam <- function(text) tags$span(class = "klam", HTML(text))
step <- function(num, title, ...) div(class = "card", div(class = "step-head", span(class = "num", num), h2(title)), ...)
slider <- function(id, label, min, max, value, step) sliderInput(id, label, min = min, max = max, value = value, step = step, ticks = FALSE, width = "100%")

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Arial, sans-serif; font-size: 13.5px; background: #f5f7fa; color: #2b2b2b; }
    .container-fluid { max-width: 1250px; }
    .top { background: #1a5276; color: white; padding: 10px 20px; margin: 0 -15px 14px -15px; }
    .top h1 { font-size: 17px; font-weight: bold; margin: 0; } .top p { font-size: 12.5px; color: #d5e8f3; margin: 2px 0 0 0; }
    .card { background: white; border: 1px solid #d5e8f3; border-radius: 6px; padding: 12px 14px; margin-bottom: 12px; }
    .step-head { display: flex; align-items: center; gap: 8px; margin-bottom: 8px; }
    .num { background: #1a5276; color: white; border-radius: 50%; width: 22px; height: 22px; display: inline-flex; align-items: center; justify-content: center; font-size: 12.5px; font-weight: bold; flex: none; }
    .step-head h2 { font-size: 14.5px; color: #1a5276; margin: 0; font-weight: bold; }
    .lay { font-size: 12.5px; color: #444; line-height: 1.4; margin: 2px 0 6px 0; }
    .klam { font-size: 11.5px; color: #1e6b52; background: #eaf6f1; border-radius: 10px; padding: 1px 8px; display: inline-block; margin-bottom: 8px; }
    .v { color: #1a5276; } .form-group { margin-bottom: 2px; } .irs--shiny .irs-min, .irs--shiny .irs-max, .irs--shiny .irs-single { display: none; }
    .pilot { display: flex; align-items: flex-end; gap: 8px; flex-wrap: wrap; font-size: 12.5px; color: #444; margin: 0 0 6px 0; }
    .pilot label { font-weight: normal; font-size: 12.5px; margin-bottom: 0; } .pilot .form-control { height: 28px; padding: 2px 6px; font-size: 12.5px; }
    .chip { border: 1px solid #2980b9; background: white; color: #1a5276; border-radius: 16px; padding: 4px 11px; font-size: 12.5px; margin: 0 4px 6px 0; }
    .tip { position: relative; display: inline-flex; align-items: center; justify-content: center; width: 15px; height: 15px; border-radius: 50%; background: #d5e8f3; color: #1a5276; font-size: 10.5px; font-weight: bold; cursor: help; margin-left: 4px; }
    .tip .tt { display: none; position: absolute; left: 20px; top: -6px; z-index: 20; width: 280px; background: #1f2d3a; color: white; font-weight: normal; font-size: 12px; line-height: 1.45; border-radius: 5px; padding: 8px 10px; text-align: left; }
    .tip:hover .tt, .tip:focus .tt { display: block; }
    .radio label { font-size: 12.5px; line-height: 1.35; } .radio b { color: #1a5276; font-size: 13.5px; }
    .answer { background: #1a5276; color: white; border-radius: 6px; padding: 12px 16px; margin-bottom: 12px; }
    .answer .big { font-size: 20px; font-weight: bold; } .answer .sub { font-size: 12.5px; color: #d5e8f3; margin-top: 3px; }
    .panel-title { color: #1a5276; font-weight: bold; font-size: 14px; } .panel-sub { font-size: 12px; color: #555; margin: 2px 0 4px 0; }
    .check { font-size: 13px; border-left: 4px solid #2980b9; background: #eaf3fb; border-radius: 4px; padding: 6px 10px; margin-top: 6px; }
    .ok { color: #1e8449; font-weight: bold; } .short { color: #c0392b; font-weight: bold; }
    details { background: white; border: 1px solid #d5e8f3; border-radius: 6px; padding: 8px 14px; margin-bottom: 8px; font-size: 13px; line-height: 1.5; }
    details summary { cursor: pointer; font-weight: bold; color: #1a5276; }
    details.step5 { margin-bottom: 12px; padding: 10px 14px; } details.step5 summary { font-size: 14.5px; }
    details.step5 summary .num { margin: 0 6px 0 4px; } details.step5 .in { max-width: 420px; }
    .skip .num { background: #9fb3c2; } .skip h2 { color: #6b7f8c; } .skip .step-head { margin-bottom: 0; }
  "))),
  div(class = "top", h1("eDNA River Cross-Section Sampling Design Tool"),
      p("How many 1-L water samples do you need to estimate the mean eDNA concentration on a river cross section? Work through the steps.")),

  ## steps 1 and 2
  fluidRow(
    column(6, step(1, "How precise must the estimate be?",
      actionButton("cv01", HTML("0.1 &nbsp;Research"), class = "chip"), actionButton("cv02", HTML("0.2 &nbsp;Management"), class = "chip"),
      actionButton("cv03", HTML("0.3 &nbsp;Pilot study"), class = "chip"),
      slider("cv_target", tip(CVT, "v_cvt", "CV<sub>target</sub> = standard error of the estimated cross-section mean &divide; the mean. It is a precision criterion, not statistical power. A 95% interval is roughly the mean &plusmn; 2 &times; CV<sub>target</sub>."), 0.05, 0.60, 0.20, 0.01),
      uiOutput("lay_target"), klam("Keel et al. reached 0.07&ndash;0.25, median 0.10"))),
    column(6, step(2, "How variable is your site?",
      slider("cv_site", tip(CVS, "v_cvs", "CV<sub>site</sub> = standard deviation &divide; mean of eDNA concentration among single 1-L samples taken across the cross section. It includes field heterogeneity and laboratory measurement error."), 0.10, 2.00, 0.55, 0.01),
      div(class = "pilot", span(style = "width:100%", "or from a pilot study:"),
          div(numericInput("pilot_mean", "mean", value = 1000, min = 0, width = "95px")), div(numericInput("pilot_sd", "SD", value = 550, min = 0, width = "95px")),
          HTML('<span class="tip" tabindex="0">i<span class="tt">Enter the mean and standard deviation of concentration among your pilot samples (any units, for example copies per liter). CV<sub>site</sub> = SD &divide; mean. Changing CV<sub>site</sub> updates the SD for the same mean.</span></span>')),
      lay("How much single 1-L samples differ from one another. <b>Higher CV<sub>site</sub> &rarr; more 1-L samples.</b>"),
      klam("Keel et al.: 0.40&ndash;1.43, median 0.55")))),

  ## step 3 and its graph
  fluidRow(
    column(4, step(3, "How will you place the samples?",
      radioButtons("strategy", NULL, selected = "random", choiceValues = c("random", "stratified"), choiceNames = list(
        HTML("<b>Simple random</b><br>Sample points picked at random across the width, or samples you cannot space evenly. Steps 1 and 2 are all you need."),
        HTML("<b>Spatially stratified</b><br>Split the width into n equal strips and take one sample at the centre of each, bank to bank. The same or fewer 1-L samples than simple random. Continue to step 4."))))),
    column(8, div(class = "card", div(class = "panel-title", textOutput("rand_title", inline = TRUE)),
      div(class = "panel-sub", HTML("The number of 1-L samples needed is where the curve drops below the dashed CV<sub>target</sub> line.")),
      plotOutput("curve_random", height = "240px")))),

  ## step 4 and its graph (spatially stratified only)
  conditionalPanel("input.strategy == 'stratified'", fluidRow(
    column(4, step(4, "Describe the cross section",
      slider("width_m", tip("Width (m)", "v_w", "Wetted width of the cross section, bank to bank."), 5, 150, 38, 1),
      klam("Keel et al.: 34&ndash;45 m, median 38 m"),
      slider("p", tip("Patchiness", "v_p", "Spatial proportion of variance, p = 1 &minus; nugget &divide; sill of a semivariogram of pilot samples: the share of the variance among samples that is spatially structured rather than random noise."), 0, 0.95, 0.50, 0.01),
      lay("Share of the variation that comes in patches rather than random noise. <b>Higher patchiness &rarr; fewer 1-L samples.</b>"), klam("Keel et al.: 0.33 Chinook salmon, 0.60 <i>C. shasta</i>"),
      slider("range_m", tip("Patch size (m)", "v_r", "Correlation range of a Mat&eacute;rn (&nu; = 1) spatial model: the distance at which the correlation between two samples falls to about 0.14. Read it from a semivariogram of pilot samples."), 1, 100, 30, 1),
      lay("Distance over which samples stay alike. <b>Fewest 1-L samples when patches are about a tenth to a fifth of the width; more when they are smaller or larger.</b>"), klam("Keel et al.: 14 m Chinook salmon, 50 m <i>C. shasta</i>"))),
    column(8, div(class = "card", div(class = "panel-title", textOutput("strat_title", inline = TRUE)),
      div(class = "panel-sub", "The same graph for evenly spaced samples, with simple random for comparison."),
      plotOutput("curve_strat", height = "330px"))))),
  conditionalPanel("input.strategy == 'random'",
    div(class = "card skip", div(class = "step-head", span(class = "num", 4), h2(HTML("Describe the cross section <span style='font-weight:normal'>&mdash; only needed for spatially stratified placement</span>"))))),

  uiOutput("answer"),

  ## step 5
  tags$details(class = "step5", tags$summary(span(class = "num", 5), "Exploring other sample sizes"),
    div(class = "in",
      slider("n_given", tip("Number of 1-L samples", "v_n", "Number of 1-L water samples collected on the cross section on one visit. The tool returns the CV of the estimated mean that this many samples achieve."), 1, 100, 9, 1),
      klam("Keel et al.: 28&ndash;35 per section, median 33"), uiOutput("check"))),

  conditionalPanel("input.strategy == 'stratified'",
    div(class = "card", div(class = "panel-title", "What a cross section like this looks like", actionButton("resim", "↻ New simulation", class = "btn-xs", style = "float:right")),
        div(class = "panel-sub", HTML("A simulated cross section with your CV<sub>site</sub>, patchiness, and patch size. Colour is concentration relative to the section mean.")),
        plotOutput("cross_section", height = "230px"))),
  div(class = "card", div(class = "panel-title", HTML("Planning maps: 1-L samples needed for any width and CV<sub>site</sub>")),
      div(class = "panel-sub", uiOutput("map_sub", inline = TRUE)),
      plotOutput("maps", height = "430px")),
  tags$details(tags$summary("How the numbers are calculated"),
    tags$p(HTML("<b>Simple random.</b> CV(n) = CV<sub>site</sub> &divide; &radic;n, so n = (CV<sub>site</sub> &divide; CV<sub>target</sub>)&sup2; (Cochran 1977).")),
    tags$p(HTML("<b>Spatially stratified.</b> Concentration across the width is modelled as a share <i>p</i> of spatially correlated variation (Mat&eacute;rn correlation, &nu; = 1, with the chosen range) plus a share 1 &minus; <i>p</i> of independent noise. The CV of the mean of n samples at the centres of n equal strips is CV<sub>site</sub> &times; &radic;[(<i>p</i> A(n) + (1 &minus; <i>p</i>)/n) &divide; (1 &minus; <i>p</i> c&#772;)], where A(n) is the variance of the difference between the sample mean and the true section mean of the correlated part, and c&#772; is the mean correlation between two random positions on the section. The tool reports the smallest n at which this CV meets CV<sub>target</sub>.")),
    tags$p(HTML("<b>Simulation.</b> A log-normal field with the chosen CV<sub>site</sub>, patchiness, and patch size in a randomly drawn channel up to 3 m deep. It is an illustration; the sample numbers do not depend on the bed shape or on the simulation."))),
  tags$details(tags$summary("What this tool does not cover"), tags$ul(
    tags$li("It concerns the mean of one cross section on one visit, not day-to-day change and not how well one section represents a longer reach."),
    tags$li(HTML("CV<sub>site</sub> from a small pilot is itself uncertain, and the samples needed scale with its square. Round up.")),
    tags$li("Patterns with depth are not modelled. In Keel et al. Chinook salmon eDNA declined with depth, so sample both surface and depth."),
    tags$li("Detection limits and the volume filtered per sample are outside its scope."))),
  tags$details(tags$summary("Where the Keel et al. values come from"),
    tags$p(HTML("The green tags give the values measured in one study, Keel et al. (in prep.): six cross sections of the Klamath River, California, sampled during one week in April 2021 with 28 to 35 one-liter samples each, for <i>Ceratonova shasta</i> and Chinook salmon (<i>Oncorhynchus tshawytscha</i>) eDNA. CV<sub>site</sub> and width are the range and median across sites and species; patchiness and patch size are one estimate per species from spatial models. They describe those sections on those days, not the Klamath River in general, and are offered as a reference point, not as defaults for other sites.")),
    tags$p(HTML("Keel, D. J., N. A. Som, M. J. Henderson, and A. P. Kinziger. In prep. <i>Spatial variability in environmental DNA concentrations on river cross sections: implications for sampling design and abundance estimation.</i>")))
)

server <- function(input, output, session) {
  tab <- reactive(matern_tables(input$range_m))
  sim_seed <- reactiveVal(42)
  observeEvent(input$resim, sim_seed(sample.int(1e6, 1)))
  observeEvent(input$cv01, updateSliderInput(session, "cv_target", value = 0.1))
  observeEvent(input$cv02, updateSliderInput(session, "cv_target", value = 0.2))
  observeEvent(input$cv03, updateSliderInput(session, "cv_target", value = 0.3))
  output$v_cvt <- renderText(sprintf("%.2f", input$cv_target)); output$v_cvs <- renderText(sprintf("%.2f", input$cv_site))
  output$v_w <- renderText(input$width_m); output$v_p <- renderText(sprintf("%.2f", input$p))
  output$v_r <- renderText(input$range_m); output$v_n <- renderText(input$n_given)

  ## step 2: CV_site and the pilot mean and SD keep each other in step
  observeEvent(input$cv_site, { m <- input$pilot_mean
    if (isTruthy(m) && m > 0 && (!isTruthy(input$pilot_sd) || abs(input$pilot_sd / m - input$cv_site) > 0.006)) updateNumericInput(session, "pilot_sd", value = signif(m * input$cv_site, 4)) }, ignoreInit = TRUE)
  observeEvent(c(input$pilot_mean, input$pilot_sd), { m <- input$pilot_mean; s <- input$pilot_sd
    if (isTruthy(m) && isTruthy(s) && m > 0 && s > 0) { cv <- min(2, max(0.1, round(s / m, 2)))
      if (abs(cv - input$cv_site) > 0.006) updateSliderInput(session, "cv_site", value = cv) } }, ignoreInit = TRUE)

  n_rand  <- reactive(ceiling((input$cv_site / input$cv_target)^2 - 1e-9))
  n_strat <- reactive(n_stratified(input$cv_site, input$cv_target, input$width_m, input$p, tab(), n_max = N_MAX))
  strat   <- reactive(input$strategy == "stratified")

  output$lay_target <- renderUI(lay(sprintf("Your estimate of the mean will typically be within about <b>%d%%</b> of the true value. <b>Lower %s &rarr; more 1-L samples.</b>", round(input$cv_target * 100), CVT)))
  output$rand_title  <- renderText(paste("Simple random: n =", fmt_n(n_rand())))
  output$strat_title <- renderText(paste("Spatially stratified: n =", fmt_n(n_strat())))
  output$map_sub <- renderUI(HTML(sprintf("At your %s%s. The white dot is your site.", CVT, if (strat()) ", patchiness, and patch size" else "")))
  output$answer <- renderUI({
    n <- if (strat()) n_strat() else n_rand()
    div(class = "answer",
        div(class = "big", sprintf("Collect %s 1-L samples, %s", if (is.na(n) || n > N_MAX) paste("more than", N_MAX) else paste("at least", n),
                                   if (strat()) "evenly spaced across the width" else "placed at random")),
        div(class = "sub", HTML(sprintf("for %s = %.2f at %s = %.2f.%s", CVT, input$cv_target, CVS, input$cv_site,
                                        if (strat()) sprintf(" Placed at random you would need %s.", fmt_n(n_rand())) else ""))))
  })
  output$check <- renderUI({
    n <- input$n_given; cr <- cv_random(input$cv_site, n)
    cg <- if (strat()) cv_stratified(input$cv_site, n, input$width_m, input$p, tab()) else cr
    HTML(sprintf('<div class="check">%d 1-L sample%s, %s: CV of the estimated mean = <b>%.2f</b>, %s%s</div>',
                 n, if (n == 1) "" else "s", if (strat()) "evenly spaced" else "placed at random", cg,
                 if (cg <= input$cv_target) sprintf('<span class="ok">meets %s</span>', CVT) else sprintf('<span class="short">misses %s</span>', CVT),
                 if (strat()) sprintf('<br><small style="color:#555">Placed at random: %.2f</small>', cr) else ""))
  })
  ## step 5 sits in a closed drop-down, so its outputs must not wait to be visible
  outputOptions(output, "check", suspendWhenHidden = FALSE); outputOptions(output, "v_n", suspendWhenHidden = FALSE)

  ## precision curves: one graph per placement
  curve_plot <- function(design) {
    is_s <- design == "stratified"; nr <- n_rand(); ns <- if (strat()) n_strat() else NA
    x_max <- min(300, max(30, ceiling(max(nr, ns, na.rm = TRUE) * 1.4 / 10) * 10)); n_seq <- 1:x_max
    df  <- data.frame(n = n_seq, CV = cv_random(input$cv_site, n_seq), type = "Simple random")
    pts <- data.frame(n = nr, CV = input$cv_target, type = "Simple random", v = -0.9, h = -0.15)
    if (is_s) {
      df  <- rbind(df, data.frame(n = n_seq, CV = cv_stratified(input$cv_site, n_seq, input$width_m, input$p, tab()), type = "Spatially stratified"))
      pts <- rbind(pts, data.frame(n = n_strat(), CV = input$cv_target, type = "Spatially stratified", v = 1.9, h = 1.15))
    }
    pts <- pts[!is.na(pts$n) & pts$n <= x_max, ]
    ggplot(df, aes(x = n, y = CV, colour = type, linetype = type, linewidth = type)) +
      geom_hline(yintercept = input$cv_target, linetype = "dashed", colour = "grey30") +
      annotate("text", x = x_max, y = input$cv_target, label = sprintf("CV[target] == '%.2f'", input$cv_target), parse = TRUE, hjust = 1, vjust = -0.4, size = 3.8, colour = "grey25") +
      geom_line() +
      geom_point(data = pts, size = 3.5, show.legend = FALSE) +
      geom_text(data = pts, aes(label = paste("n =", n), vjust = v, hjust = h), size = 4, show.legend = FALSE) +
      scale_colour_manual(values = c("Simple random" = "#0072B2", "Spatially stratified" = "#009E73"), name = NULL) +
      scale_linetype_manual(values = c("Simple random" = "solid", "Spatially stratified" = "22"), name = NULL) +
      scale_linewidth_manual(values = c("Simple random" = if (is_s) 0.8 else 1.5, "Spatially stratified" = 1.5), name = NULL) +
      scale_x_continuous(limits = c(1, x_max), expand = expansion(mult = c(0, 0.01))) +
      coord_cartesian(ylim = c(0, min(input$cv_site * 1.1, 2.2))) +
      labs(x = "Number of 1-L samples (n)", y = "CV of the estimated mean") +
      theme_bw(base_size = 13) + theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"), legend.position = if (is_s) "bottom" else "none")
  }
  output$curve_random <- renderPlot(curve_plot("random"))
  output$curve_strat  <- renderPlot(curve_plot("stratified"))

  output$cross_section <- renderPlot({
    W <- input$width_m
    gd <- matern_field(W, cv_site = input$cv_site, p = input$p, range_m = input$range_m, seed = sim_seed())
    outline <- data.frame(x = seq(0, W, length.out = 300)); outline$y <- -bed_profile(sim_seed())(outline$x / W)
    top <- max(2, quantile(gd$rel, 0.98, na.rm = TRUE))
    ggplot(gd[!is.na(gd$rel), ], aes(x = x, y = -y, fill = rel)) +
      geom_tile() +
      geom_line(data = outline, aes(x = x, y = y), inherit.aes = FALSE, colour = "grey20", linewidth = 1.2) +
      annotate("segment", x = 0, xend = W, y = 0, yend = 0, colour = "grey20", linewidth = 1.2) +
      scale_fill_gradientn(colours = c("#eef5fb", "#9ecae1", "#3182bd", "#08306b"), limits = c(0, top), oob = scales::squish,
                           name = "× section\nmean", labels = function(x) paste0(x, "×")) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0.02, 0), labels = abs) +
      labs(x = "Distance from left bank (m)", y = "Depth (m)") +
      theme_bw(base_size = 13) + theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"))
  })

  ## planning maps: width on the x axis, CV_site on the y axis, n range written in each colour block
  output$maps <- renderPlot({
    cv_seq <- seq(0.10, 2.00, by = 0.025); w_seq <- seq(5, 150, by = 2.5)
    mk <- function(p, design) n_required_grid(cv_seq, w_seq, input$cv_target, p, tab(), design, n_max = N_MAX)
    mark <- list(cv_site = input$cv_site, width = input$width_m)
    br <- c(0, 3, 5, 10, 15, 20, 30, 50, 100, 250, N_MAX)
    a <- contour_panel(mk(0, "random"), "Simple random (width has no effect)", br = br, mark = mark, flip = TRUE, band_text = TRUE) + labs(y = expression(bold(CV[site])))
    if (!strat()) return(a + theme(legend.position = "none"))
    b <- contour_panel(mk(input$p, "stratified"), sprintf("Spatially stratified, patchiness %.2f, patch size %s m", input$p, input$range_m), br = br, mark = mark, flip = TRUE, band_text = TRUE) + labs(y = expression(bold(CV[site])))
    ggpubr::ggarrange(a, b, ncol = 2, legend = "none")
  })
}

shinyApp(ui, server)
