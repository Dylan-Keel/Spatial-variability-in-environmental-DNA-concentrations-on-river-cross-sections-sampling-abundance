#### Spatial variability in environmental DNA concentrations on river cross sections ####
## Keel, D. J., N. A. Som, M. J. Henderson, and A. P. Kinziger
##
## Complete analysis for the manuscript: data preparation, GLMMs, sample size
## calculations, results tables, manuscript figures, and supporting information
## figures. Run from the R_code folder:
##   Rscript Main_Analysis_and_Figures.R
##
## Run time is about 10 minutes. INLA fits use a single thread because
## multithreaded INLA fits are not numerically reproducible. The model
## selection cross-validation (100 repeats of 10-fold cross-validation over 16
## candidate models, 16,000 fits per species) is not run by default
## (RUN.KFOLD = FALSE); its saved results are read from data/. When enabled it
## runs in parallel on all but two cores (approximately 2 hours per species on
## 16 cores), writes one file per repeat, and can be resumed.
##
## Sections
##   0. setup
##   1. data
##   2. C. shasta GLMMs
##   3. O. tshawytscha GLMMs
##   4. estimating mean eDNA concentration on a cross section
##   5. results tables
##   6. manuscript figures
##   7. supporting information figures
##   8. save

#### 0. setup ####

library(INLA)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(sp)
library(gstat)
library(mgcv)
library(inlatools)   # plot methods for the dispersion and distribution checks
library(INLAutils)   # residual plots

## Helper functions
source("HighstatLibV13.R")               # Mypairs (Zuur et al. 2009)
source("INLA_plotting_functions.R")      # inlaPredict and inlaPlot
source("inlatools_modified.R")           # dispersion_check and distribution_check for stacked models
source("Sampling_Design_Functions.R")    # sample size under random and stratified placement

## Seed for the simulation-based model checks. INLA fits use one thread for
## reproducibility; the Poisson model of the C. shasta counts does not run on a
## single thread and is refit on two (its result is stable across runs).
## Parallelism is applied across fits in the cross-validation.
set.seed(20230126)
SEED <- 20230126L
inla.setOption(num.threads = "1:1")
fit.inla <- function(...) tryCatch(inla(...), error = function(e) {
  message("inla could not run this fit on one thread; refitting on two"); inla(..., num.threads = "2:1") })
N.WORKERS <- max(1, parallel::detectCores() - 2)
RUN.KFOLD <- FALSE      # TRUE repeats the cross-validation instead of reading the saved results

## Output folders
dir.create("../results", showWarnings = FALSE)
for (d in c("../figures", "../figures/pdf", "../figures/supporting_information", "../figures/supporting_information/pdf"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
## Base-graphics diagnostics (pairs plots, meshes) are written to one PDF
pdf("../results/diagnostic_plots.pdf", width = 9, height = 7)

#### 1. data ####

## One row per water sample with the ddPCR results
CSD2 <- read.csv("../data/CSD.csv")

## Copies per liter of water from copies per reaction, elution volume, and
## volume filtered. Upper and lower are the ddPCR Poisson confidence limits.
vol.cor <- function(c.well) replace_na(round(c.well * CSD2$Elution.volum..uL. / 15 * CSD2$Volume..ml. / 1000, 0), 0)
CSD2$C.shasta.vol.cor.conc       <- vol.cor(CSD2$C.shasta.c.well)
CSD2$Chinook.vol.cor.conc        <- vol.cor(CSD2$Chinook.c.well)
CSD2$C.shasta.vol.cor.conc.upper <- vol.cor(CSD2$C.shasta.c.well.upper)
CSD2$Chinook.vol.cor.conc.upper  <- vol.cor(CSD2$Chinook.c.well.upper)
CSD2$C.shasta.vol.cor.conc.lower <- vol.cor(CSD2$C.shasta.c.well.lower)
CSD2$Chinook.vol.cor.conc.lower  <- vol.cor(CSD2$Chinook.c.well.lower)

## Collinearity among the continuous covariates
MyVar <- c("D.thal.m", "DistanceFromShore.m.z", "D.M.")
Mypairs(CSD2[, MyVar])

## Covariates scaled to mean zero and unit standard deviation; the scaling
## parameters are kept for back-transformation
CSD2.scale <- dplyr::select(CSD2, all_of(MyVar))
names(CSD2.scale) <- paste(MyVar, ".z", sep = "")
scale.par <- do.call(data.frame, list(mean = apply(CSD2.scale, 2, mean),
                                      sd   = apply(CSD2.scale, 2, sd),
                                      min  = apply(CSD2.scale, 2, min),
                                      max  = apply(CSD2.scale, 2, max)))
CSD2.scale <- as.data.frame(scale(CSD2.scale))
## CSD.csv already contains the scaled depth (D.M..z) used by the models; only
## new columns are added
CSD2.scale <- CSD2.scale[, setdiff(names(CSD2.scale), names(CSD2)), drop = FALSE]
CSD2 <- cbind(CSD2, CSD2.scale)

## Collinearity among all covariates. Depth and sample type, and distance from
## shore and velocity, are collinear and are not included in the same model.
CSD2$Distance.Left.Bank.m[CSD2$Distance.Left.Bank.m == 0] <- 0.01
Mypairs(CSD2[, c("DistanceFromShore.m.z", "D.M.", "V.M.", "Type01", "SiteN", "D.thal.m")])

## Site numbered 1-6 from upstream to downstream, depth squared, and the
## responses as integer counts of copies per liter
CSD2$SiteN  <- recode_factor(as.factor(CSD2$Site), KLM = 1, I5 = 2, TOH = 3, BVR = 4, KMN = 5, SV = 6)
CSD2$Site   <- as.numeric(CSD2$SiteN)
CSD2$Depth2 <- as.numeric(CSD2$D.M..z) * as.numeric(CSD2$D.M..z)
y.cs <- CSD2$C.shasta.vol.cor.conc
y.ot <- CSD2$Chinook.vol.cor.conc

## Sample coordinates in km. The river is represented as a straight line with
## distance downstream on the x axis and distance from the left bank on the y
## axis, so between-site distances are to scale.
xy <- with(CSD2, data.frame(Sample.ID, Meters / 1000 / 111.111, Distance.Left.Bank.m / 1000 / 111.111))
coordinates(xy) <- c("Meters.1000.111.111", "Distance.Left.Bank.m.1000.111.111")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
utm <- data.frame(spTransform(xy, CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")))
names(utm) <- c("Sample.ID", "Xkm", "Ykm", "drop")
CSD2 <- merge(CSD2, utm)

## Site summaries reported in the text
site.summary <- CSD2 %>% group_by(Site) %>%
  summarise(n = n(),
            mean.cs = mean(C.shasta.vol.cor.conc), sd.cs = sd(C.shasta.vol.cor.conc), cv.cs = sd.cs / mean.cs,
            mean.ot = mean(Chinook.vol.cor.conc),  sd.ot = sd(Chinook.vol.cor.conc),  cv.ot = sd.ot / mean.ot,
            width.m = max(Distance.Left.Bank.m), .groups = "drop")
print(site.summary, width = Inf)

## Fixed effects of the candidate models and the random effect of site
all.form <- Copies ~ Depth2 + D.M..z + DistanceFromShore.m.z.z + Type01 + V.M..z + D.thal.m.z + f(Site, model = "iid")

## Design matrix for the INLA stacks (main effects only)
make.Xm <- function(form, data) {
  terms <- gsub(" ", "", as.character(form[3]))
  terms <- unlist(strsplit(terms, "\\+"))
  terms <- terms[-grep("iid", terms)]
  terms <- c(terms, "Site")
  Xm <- model.matrix(formula(paste("~", paste(terms, collapse = "+"))), data = data)
  as.data.frame(Xm)
}

## The sixteen candidate models. Each includes the spatial field w; the random
## effect of site is added at fitting.
candidate.forms <- list(
  sf.d2.d.dm = Copies ~ -1 + Depth2 + D.M..z + DistanceFromShore.m.z.z + f(w, model = spde),
  sf.d2.d.dt = Copies ~ -1 + Depth2 + D.M..z + D.thal.m.z + f(w, model = spde),
  sf.d2.d.v  = Copies ~ -1 + Depth2 + D.M..z + V.M..z + f(w, model = spde),
  sf.d2.d    = Copies ~ -1 + Depth2 + D.M..z + f(w, model = spde),
  sf.d       = Copies ~ -1 + D.M..z + f(w, model = spde),
  sf.d.v     = Copies ~ -1 + D.M..z + V.M..z + f(w, model = spde),
  sf.d.dt    = Copies ~ -1 + D.M..z + D.thal.m.z + f(w, model = spde),
  sf.d.dm    = Copies ~ -1 + DistanceFromShore.m.z.z + D.M..z + f(w, model = spde),
  sf.dm      = Copies ~ -1 + DistanceFromShore.m.z.z + f(w, model = spde),
  sf.dt      = Copies ~ -1 + D.thal.m.z + f(w, model = spde),
  sf.v       = Copies ~ -1 + V.M..z + f(w, model = spde),
  sf.v.t     = Copies ~ -1 + V.M..z + Type01 + f(w, model = spde),
  sf.dm.t    = Copies ~ -1 + Type01 + DistanceFromShore.m.z.z + f(w, model = spde),
  sf.dt.t    = Copies ~ -1 + Type01 + D.thal.m.z + f(w, model = spde),
  sf.t       = Copies ~ -1 + Type01 + f(w, model = spde),
  sf.null    = Copies ~ -1 + f(w, model = spde)
)
model.names <- c("Depth-Squared + Depth + Distance From Shore",
                 "Depth-Squared + Depth + Distance From Thalweg",
                 "Depth-Squared + Depth + Velocity",
                 "Depth-Squared + Depth",
                 "Depth",
                 "Depth + Velocity",
                 "Depth + Distance From Thalweg",
                 "Depth + Distance From Shore",
                 "Distance From Shore",
                 "Distance From Thalweg",
                 "Velocity",
                 "Velocity + Sample Type",
                 "Distance From Shore + Sample Type",
                 "Distance From Thalweg + Sample Type",
                 "Sample Type",
                 "NULL")
n.covariates <- c(3, 3, 3, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 0)

## Spatial negative binomial GLMM with site as an iid random effect
fit.spatial <- function(form, stack) {
  fit.inla(formula = update(form, . ~ . + f(Site, model = "iid")),
       family = "nbinomial",
       data = inla.stack.data(stack),
       control.compute = list(waic = TRUE, config = TRUE, return.marginals.predictor = TRUE),
       control.predictor = list(link = 1, compute = TRUE, A = inla.stack.A(stack)))
}

## Pseudo-R2 of fitted versus observed counts
pseudo.r2 <- function(y, fitted) 1 - sum((y - fitted)^2) / sum((y - mean(y))^2)

## Pearson residuals
pearson.resid <- function(y, fitted) (y - fitted) / sqrt(fitted)

## Semivariogram of residuals at the scale of a cross section
resid.variogram <- function(resid, X, Y, cutoff) {
  MyData <- data.frame(D1 = as.numeric(resid), X = as.numeric(X), Y = as.numeric(Y))
  coordinates(MyData) <- c("X", "Y")
  gstat::variogram(D1 ~ X + Y, data = MyData, cutoff = cutoff, cressie = TRUE)
}

## Summary of the cross-validation results: the number of folds in which each
## model had the highest validation R2, and the mean validation R2
## (transferability).
selection.table <- function(metawaic) {
  msr <- matrix(nrow = 10, ncol = 100)
  for (k in 1:100) for (i in 1:10) {
    v <- as.numeric(lapply(metawaic, function(x) x[i, k]))
    msr[i, k] <- if (all(is.na(v))) NA else which.max(v)
  }
  times.selected <- tabulate(msr[!is.na(msr)], nbins = 16)
  transferability <- sapply(metawaic, function(x) mean(unlist(x), na.rm = TRUE))
  out <- data.frame(model = model.names, n.covariates = n.covariates,
                    times.selected = times.selected, weight = times.selected / 1000 * 100,
                    transferability = transferability)
  ## Most parsimonious model(s) within 0.05 of the highest transferability
  close <- out$transferability >= max(out$transferability) - 0.05
  out$selected <- close & out$n.covariates == min(out$n.covariates[close])
  out[order(-out$transferability), ]
}

#### 2. C. shasta GLMMs ####

## Non-spatial models with a relatively full set of covariates to select the
## count distribution: Poisson, negative binomial, and zero-inflated versions
CSD2$Copies <- y.cs
base.form.d2 <- Copies ~ Depth2 + D.M..z + DistanceFromShore.m.z.z + f(Site, model = "iid")
fit.nonspatial <- function(family) {
  fit.inla(base.form.d2, family = family, data = CSD2,
       control.compute = list(waic = TRUE, cpo = TRUE, po = TRUE, config = TRUE, return.marginals.predictor = TRUE),
       control.predictor = list(link = 1, compute = TRUE))
}
nonspat.cs       <- fit.nonspatial("poisson")
nonspat.nbn.cs   <- fit.nonspatial("nbinomial")
nonspat.zip.cs   <- fit.nonspatial("zeroinflatedpoisson1")
nonspat.zinbn.cs <- fit.nonspatial("zeroinflatednbinomial1")

## Overdispersion checks (Figures S1 and S2). The Poisson model is
## overdispersed; the negative binomial model is not.
set.seed(SEED); dc.pois.cs <- dispersion_check(nonspat.cs, y.cs)
set.seed(SEED); dc.nbn.cs  <- dispersion_check(nonspat.nbn.cs, y.cs)
c(data = dc.pois.cs$data, model = mean(dc.pois.cs$model), P.overdispersed = mean(dc.pois.cs$data > dc.pois.cs$model))
c(data = dc.nbn.cs$data,  model = mean(dc.nbn.cs$model),  P.overdispersed = mean(dc.nbn.cs$data > dc.nbn.cs$model))

## Spatial structure in the Pearson residuals of the negative binomial model.
## The downstream coordinate is divided by 50 so that between-site distances do
## not dominate within-site distances; the semivariogram is cut off at 60 m to
## show autocorrelation at the scale of a cross section (Figures S3 and S4).
resid.nbn.cs <- pearson.resid(y.cs, nonspat.nbn.cs$summary.fitted.values[, "mean"])
V.nbn.cs <- resid.variogram(resid.nbn.cs, utm$Xkm / 50, utm$Ykm, cutoff = 0.06)
nugget.cs <- 9600; sill.cs <- 24200; range.cs <- 0.05   # read from the semivariogram (Figure S3)
spatial.proportion.cs <- 1 - nugget.cs / sill.cs            # 60.3 percent of the residual variance is spatial

## Mesh and SPDE for the spatial random field. The prior range is the
## semivariogram range (50 m).
Loc.cs <- cbind(utm$Xkm / 50, utm$Ykm)
RangeGuess.cs <- 0.05
Hull.cs <- inla.nonconvex.hull(Loc.cs, convex = -0.085)
MaxEdge.cs <- RangeGuess.cs / 6
mesh.cs <- inla.mesh.2d(boundary = Hull.cs, max.edge = c(0.01, 10000) * MaxEdge.cs, cutoff = MaxEdge.cs, max.n = 1500)
mesh.cs$n
plot(mesh.cs, asp = 1, main = "C. shasta mesh"); points(Loc.cs, col = 2, pch = 1, cex = 1)
A.cs <- inla.spde.make.A(mesh.cs, loc = Loc.cs)
spde.cs <- inla.spde2.pcmatern(mesh.cs, prior.range = c(RangeGuess.cs, 0.05), prior.sigma = c(1, 0.05))
w.index.cs <- inla.spde.make.index("w", n.spde = spde.cs$n.spde)

## INLA stack: response, covariates, and projector matrix
Xm.cs <- make.Xm(all.form, CSD2)
stack.cs <- inla.stack(remove.unused = TRUE, tag = "Fit",
                       data = list(Copies = y.cs),
                       A = list(1, 1, A.cs),
                       effects = list(Intercept = rep(1, nrow(CSD2)), Xm = Xm.cs[, -1], w = w.index.cs))
fitIndex.cs <- inla.stack.index(stack.cs, tag = "Fit")$data

## Full model with the spatial field, used to confirm that the field accounts
## for the residual spatial autocorrelation (Figure S5)
spde <- spde.cs
spat.full.cs <- fit.spatial(candidate.forms$sf.d2.d.dm, stack.cs)
set.seed(SEED); dc.full.cs <- dispersion_check(spat.full.cs, y.cs)
resid.full.cs <- pearson.resid(y.cs, spat.full.cs$summary.fitted.values[fitIndex.cs, "mean"])
V.full.cs <- resid.variogram(resid.full.cs, utm$Xkm, utm$Ykm, cutoff = 0.06)   # Figure S5

## Model selection by 100 repeats of 10-fold cross-validation over the sixteen
## candidate models. Saved results are read from data/; RUN.KFOLD = TRUE
## repeats the cross-validation (parallel, resumable). CV.CORRECTED = FALSE
## reproduces the fold definition used for the saved results; TRUE fits each
## model to nine folds and scores it on the held-out fold.
CV.CORRECTED <- FALSE
one.repeat <- function(p, y, mesh, Loc, spde, w.index, Xm, forms, corrected) {
  set.seed(SEED + p)
  n <- length(y)
  g <- sample(cut(seq(n), n * cumsum(c(0, rep(0.1, 10))), labels = paste0("val", 1:10)))
  r2 <- matrix(NA_real_, 10, length(forms))
  for (k in 1:10) {
    index_val   <- which(g == paste0("val", k))
    index_train <- if (corrected) which(g != paste0("val", k)) else which(g != paste0("train", k))
    A.train <- inla.spde.make.A(mesh = mesh, loc = Loc[index_train, ])
    A.val   <- inla.spde.make.A(mesh = mesh, loc = Loc[index_val, ])
    stack.train <- inla.stack(remove.unused = TRUE, tag = "train", data = list(Copies = y[index_train]),
                              A = list(1, 1, A.train),
                              effects = list(Intercept = rep(1, length(index_train)), Xm = Xm[index_train, -1], w = w.index))
    stack.val <- inla.stack(remove.unused = TRUE, tag = "val",
                            data = list(Copies = if (corrected) rep(NA_real_, length(index_val)) else y[index_val]),
                            A = list(1, 1, A.val),
                            effects = list(Intercept = rep(1, length(index_val)), Xm = Xm[index_val, -1], w = w.index))
    join.stack <- inla.stack(stack.train, stack.val)
    index_inla_val <- inla.stack.index(join.stack, "val")$data
    for (i in seq_along(forms)) {
      model <- tryCatch(fit.spatial(forms[[i]], join.stack), error = function(e) NULL)
      if (is.null(model)) next
      r2[k, i] <- pseudo.r2(y[index_val], model$summary.fitted.values$mean[index_inla_val])
    }
  }
  r2
}
kfold.cv <- function(tag, y, mesh, Loc, spde, w.index, Xm) {
  out.dir <- "../results/cross_validation"; dir.create(out.dir, showWarnings = FALSE)
  done <- as.integer(sub(".*_rep(\\d+)\\.rds$", "\\1", list.files(out.dir, pattern = paste0("^", tag, "_rep\\d+\\.rds$"))))
  to.run <- setdiff(1:100, done)
  if (length(to.run)) {
    cl <- parallel::makeCluster(N.WORKERS)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(INLA)); inla.setOption(num.threads = "1:1"); NULL })
    parallel::clusterExport(cl, c("one.repeat", "fit.spatial", "pseudo.r2", "candidate.forms", "SEED", "CV.CORRECTED",
                                  "out.dir", "tag", "y", "mesh", "Loc", "spde", "w.index", "Xm"), envir = environment())
    parallel::parLapply(cl, to.run, function(p) {
      assign("spde", spde, envir = globalenv())    # formulas reference the spde object by name
      r2 <- one.repeat(p, y, mesh, Loc, spde, w.index, Xm, candidate.forms, CV.CORRECTED)
      saveRDS(r2, file.path(out.dir, sprintf("%s_rep%03d.rds", tag, p)))
      p
    })
  }
  reps <- lapply(1:100, function(p) readRDS(file.path(out.dir, sprintf("%s_rep%03d.rds", tag, p))))
  lapply(seq_along(candidate.forms), function(i) sapply(reps, function(r) r[, i]))
}
if (RUN.KFOLD) {
  C.shasta.metawaic <- kfold.cv("C.shasta", y.cs, mesh.cs, Loc.cs, spde.cs, w.index.cs, Xm.cs)
  save(C.shasta.metawaic, file = "../data/C.shasta.metar.new.rda")
}
load("../data/C.shasta.metar.new.rda")
selection.cs <- selection.table(C.shasta.metawaic)
selection.cs

## The null model is the most parsimonious model within 0.05 of the highest
## transferability and is the selected model
Cs.inla.null <- fit.spatial(candidate.forms$sf.null, stack.cs)
summary(Cs.inla.null)

## Goodness of fit of the selected model: dispersion, posterior predictive
## check, standardized residuals, and pseudo-R2 (Figures S6 and S7)
set.seed(SEED); dc.null.cs <- dispersion_check(Cs.inla.null, y.cs)
ds.null.cs <- distribution_check(Cs.inla.null, y.cs, seed = SEED)
fitted.cs  <- Cs.inla.null$summary.fitted.values[fitIndex.cs, ]
resid.null.cs <- pearson.resid(y.cs, fitted.cs$mean)
r2.cs <- pseudo.r2(y.cs, fitted.cs$mean)     # 0.49

## Spatial field: range, kappa, and the implied Matern correlation between
## samples on a cross section (Figures S8 and S9)
SpatField.cs <- inla.spde2.result(inla = Cs.inla.null, name = "w", spde = spde.cs, do.transfer = TRUE)
Kappa.cs   <- inla.emarginal(function(x) x, SpatField.cs$marginals.kappa[[1]])
Sigma_u.cs <- inla.emarginal(function(x) sqrt(x), SpatField.cs$marginals.variance.nominal[[1]])
Range.cs   <- inla.emarginal(function(x) x, SpatField.cs$marginals.range.nominal[[1]])
matern <- function(d, kappa) { out <- (kappa * d) * besselK(kappa * d, 1); out[d == 0] <- 1; out }
cor.cs <- data.frame(d.km = seq(0, 2, length = 10000)) %>% mutate(Cor.M = matern(d.km, Kappa.cs))
matern.at.cs <- matern(c(0.005, 0.02, 0.04), Kappa.cs)     # 5, 20, and 40 m apart
c(Kappa = Kappa.cs, Range = Range.cs, Sigma_u = Sigma_u.cs, matern.at.cs)

#### 3. O. tshawytscha GLMMs ####

## Same workflow as for C. shasta. Distances are not rescaled, the prior range
## is smaller, and the random effect of site has a penalized complexity prior.
CSD2$Copies <- y.ot
pcprior <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))
base.form.d2 <- Copies ~ Depth2 + D.M..z + DistanceFromShore.m.z.z + f(Site, model = "iid", hyper = pcprior)
nonspat.ot       <- fit.nonspatial("poisson")
nonspat.nbn.ot   <- fit.nonspatial("nbinomial")
nonspat.zip.ot   <- fit.nonspatial("zeroinflatedpoisson1")
nonspat.zinbn.ot <- fit.nonspatial("zeroinflatednbinomial1")
set.seed(SEED); dc.pois.ot <- dispersion_check(nonspat.ot, y.ot)
set.seed(SEED); dc.nbn.ot  <- dispersion_check(nonspat.nbn.ot, y.ot)
c(data = dc.pois.ot$data, model = mean(dc.pois.ot$model), P.overdispersed = mean(dc.pois.ot$data > dc.pois.ot$model))
c(data = dc.nbn.ot$data,  model = mean(dc.nbn.ot$model),  P.overdispersed = mean(dc.nbn.ot$data > dc.nbn.ot$model))

## Residual spatial autocorrelation. About one third of the residual variance
## is spatial, over 0 to 60 m.
resid.nbn.ot <- pearson.resid(y.ot, nonspat.nbn.ot$summary.fitted.values[, "mean"])
V.nbn.ot <- resid.variogram(resid.nbn.ot, utm$Xkm, utm$Ykm, cutoff = 0.05)
nugget.ot <- 110; sill.ot <- 165; range.ot <- 0.03
spatial.proportion.ot <- 1 - nugget.ot / sill.ot            # 33.3 percent

## Mesh and SPDE
Loc.ot <- cbind(utm$Xkm, utm$Ykm)
RangeGuess.ot <- 0.03
Hull.ot <- inla.nonconvex.hull(Loc.ot, convex = -0.085)
MaxEdge.ot <- RangeGuess.ot / 6
mesh.ot <- inla.mesh.2d(boundary = Hull.ot, max.edge = c(0.01, 10000) * MaxEdge.ot, cutoff = MaxEdge.ot, max.n = 1500)
mesh.ot$n
plot(mesh.ot, asp = 1, main = "O. tshawytscha mesh"); points(Loc.ot, col = 2, pch = 1, cex = 1)
A.ot <- inla.spde.make.A(mesh.ot, loc = Loc.ot)
spde.ot <- inla.spde2.pcmatern(mesh.ot, prior.range = c(RangeGuess.ot / 10, 0.005), prior.sigma = c(.1, .05))
w.index.ot <- inla.spde.make.index("w", n.spde = spde.ot$n.spde)

Xm.ot <- make.Xm(all.form, CSD2)
stack.ot <- inla.stack(remove.unused = TRUE, tag = "Fit",
                       data = list(Copies = y.ot),
                       A = list(1, 1, A.ot),
                       effects = list(Intercept = rep(1, nrow(CSD2)), Xm = Xm.ot[, -1], w = w.index.ot))
fitIndex.ot <- inla.stack.index(stack.ot, tag = "Fit")$data

## Full model with the spatial field. Its predicted depth and distance to
## shore effects are shown in Figure 6.
spde <- spde.ot
spat.full.ot <- fit.spatial(candidate.forms$sf.d2.d.dm, stack.ot)
set.seed(SEED); dc.full.ot <- dispersion_check(spat.full.ot, y.ot)
resid.full.ot <- pearson.resid(y.ot, spat.full.ot$summary.fitted.values[fitIndex.ot, "mean"])
V.full.ot <- resid.variogram(resid.full.ot, utm$Xkm, utm$Ykm, cutoff = 0.05)
## Linear combinations are built from the fixed effects only; the site column
## is excluded from the data passed to inlaPredict
modrando <- "Site"
GRSpred <- inlaPredict(spat.full.ot, mdl.stack = stack.ot, spatialMod = TRUE,
                       orig.dat = CSD2[, setdiff(names(CSD2), "Site")], data.pts = 192)

## Model selection
if (RUN.KFOLD) {
  Chinook.metawaic <- kfold.cv("Chinook", y.ot, mesh.ot, Loc.ot, spde.ot, w.index.ot, Xm.ot)
  save(Chinook.metawaic, file = "../data/Chinook.metar.new.rda")
}
load("../data/Chinook.metar.new.rda")
selection.ot <- selection.table(Chinook.metawaic)
selection.ot

## Depth and sample type are the most parsimonious models within 0.05 of the
## top model. Sample type is a proxy for depth, so the depth model is selected.
## The sample type, null, and top (distance to shore + sample type) models are
## fit for comparison.
Xm2.ot <- make.Xm(Copies ~ D.M..z + Type01 + DistanceFromShore.m.z.z + f(Site, model = "iid"), CSD2)
stack2.ot <- inla.stack(remove.unused = TRUE, tag = "Fit",
                        data = list(Copies = y.ot),
                        A = list(1, 1, A.ot),
                        effects = list(Intercept = rep(1, nrow(CSD2)), Xm = Xm2.ot[, -1], w = w.index.ot))
fitIndex2.ot <- inla.stack.index(stack2.ot, tag = "Fit")$data
Chin.inla.d    <- fit.spatial(candidate.forms$sf.d, stack2.ot)
Chin.inla.t    <- fit.spatial(candidate.forms$sf.t, stack2.ot)
Chin.inla.null <- fit.spatial(candidate.forms$sf.null, stack2.ot)
Chin.inla.top  <- fit.spatial(candidate.forms$sf.dm.t, stack2.ot)
summary(Chin.inla.d)

## Goodness of fit of the selected model
set.seed(SEED); dc.d.ot <- dispersion_check(Chin.inla.d, y.ot)
ds.d.ot <- distribution_check(Chin.inla.d, y.ot, seed = SEED)
fitted.ot      <- Chin.inla.d$summary.fitted.values[fitIndex2.ot, ]
fitted.t.ot    <- Chin.inla.t$summary.fitted.values[fitIndex2.ot, ]
fitted.null.ot <- Chin.inla.null$summary.fitted.values[fitIndex2.ot, ]
resid.d.ot <- pearson.resid(y.ot, fitted.ot$mean)
r2.ot      <- pseudo.r2(y.ot, fitted.ot$mean)        # 0.45
r2.t.ot    <- pseudo.r2(y.ot, fitted.t.ot$mean)
r2.null.ot <- pseudo.r2(y.ot, fitted.null.ot$mean)   # random effects only
c(depth = r2.ot, type = r2.t.ot, null = r2.null.ot,
  depth.explains = r2.ot - r2.null.ot, type.explains = r2.t.ot - r2.null.ot)

## Depth effect: proportional change in concentration per one standard
## deviation (48 cm) increase in depth, exp(beta) - 1
beta.depth <- Chin.inla.d$summary.fixed["D.M..z", ]
depth.effect <- c(mean = 1 - exp(beta.depth$mean), lower = 1 - exp(beta.depth$`0.025quant`), upper = 1 - exp(beta.depth$`0.975quant`))
depth.effect
sd(CSD2$D.M.)
## Distance to shore in the top model: posterior probability of a negative
## effect
beta.shore.top <- Chin.inla.top$summary.fixed["DistanceFromShore.m.z.z", ]
P.shore.negative <- inla.pmarginal(0, Chin.inla.top$marginals.fixed$DistanceFromShore.m.z.z)
P.shore.negative

## Spatial field (Figures S10 and S11)
SpatField.ot <- inla.spde2.result(inla = Chin.inla.d, name = "w", spde = spde.ot, do.transfer = TRUE)
Kappa.ot   <- inla.emarginal(function(x) x, SpatField.ot$marginals.kappa[[1]])
Sigma_u.ot <- inla.emarginal(function(x) sqrt(x), SpatField.ot$marginals.variance.nominal[[1]])
Range.ot   <- inla.emarginal(function(x) x, SpatField.ot$marginals.range.nominal[[1]])
cor.ot <- data.frame(d.km = seq(0, 0.04, length = 100)) %>% mutate(Cor.M = matern(d.km, Kappa.ot))
matern.at.ot <- matern(c(0.005, 0.02, 0.04), Kappa.ot)
c(Kappa = Kappa.ot, Range = Range.ot, Sigma_u = Sigma_u.ot, matern.at.ot)

#### 4. estimating mean eDNA concentration on a cross section ####

## Number of 1 L samples needed to estimate the mean concentration on a cross
## section to a target precision, n = (CV_site / CV_target)^2 (Cochran 1977).
## Site CVs are computed from the samples on each section, so this is the
## precision of the section mean under simple random placement.
CV.targets <- c(0.1, 0.2, 0.3, 0.4, 0.5)
n.required <- site.summary %>% select(Site, cv.cs, cv.ot) %>%
  tidyr::crossing(CV = CV.targets) %>%
  mutate(n.cs = ceiling((cv.cs / CV)^2), n.ot = ceiling((cv.ot / CV)^2))
n.required %>% filter(CV == 0.2) %>% summarise(across(c(n.cs, n.ot), list(min = min, max = max, mean = mean)))

## Evenly spaced (spatially stratified) samples give higher precision than
## random placement when the field is autocorrelated. The gain depends on the
## spatial proportion of variance (p, from the semivariogram) and on the
## correlation range relative to the section width. C. shasta uses the
## semivariogram range (50 m), since its fitted field range reflects the
## between-site scale; O. tshawytscha uses the range of its fitted field.
range.ot.m <- sqrt(8) / Kappa.ot * 1000
species.par <- list("C. shasta"      = list(p = spatial.proportion.cs, range_m = range.cs * 1000),
                    "O. tshawytscha" = list(p = spatial.proportion.ot, range_m = range.ot.m))
matern.tab <- lapply(species.par, function(s) matern_tables(s$range_m))

## Sample sizes for a target standard error (Cochran finite population form),
## reported in the Supporting Information
cochran.n <- function(x, N, SE) { B <- 2 * SE; V <- B^2 / 4; s2 <- var(x); ceiling(s2 / (V + s2 / N)) }
sample.size.supp <- CSD2 %>% group_by(Site) %>%
  summarise(N = n(),
            mean.cs = mean(C.shasta.vol.cor.conc), var.cs = var(C.shasta.vol.cor.conc), CV.cs = sqrt(var.cs) / mean.cs * 100,
            n.cs.SE15000 = cochran.n(C.shasta.vol.cor.conc, N, 15000), n.cs.SE3000 = cochran.n(C.shasta.vol.cor.conc, N, 3000),
            mean.ot = mean(Chinook.vol.cor.conc), var.ot = var(Chinook.vol.cor.conc), CV.ot = sqrt(var.ot) / mean.ot * 100,
            n.ot.SE214 = cochran.n(Chinook.vol.cor.conc, N, 214), n.ot.SE107 = cochran.n(Chinook.vol.cor.conc, N, 107),
            .groups = "drop")

#### 5. results tables ####

## Molecular results reported in the text: field blanks and the proportion of
## samples above the limits of detection and quantification
ddpcr <- read.csv("../data/MASTER.FINAL.ddPCRdata.Cshasta2021.csv")
field.blanks <- ddpcr[grep("FB", ddpcr$Sample.ID), c("Sample.ID", "Chinook.c.well", "C.shasta.c.well")]
LOD <- c(cs = 6.3, ot = 11.8); LOQ <- c(cs = 18, ot = 16)     # copies per reaction, from the standard-curve analysis (Supporting Information)
molecular.summary <- data.frame(
  species = c("C. shasta", "O. tshawytscha"),
  LOD = LOD, LOQ = LOQ,
  pct.above.LOD = c(mean(CSD2$C.shasta.c.well >= LOD["cs"]), mean(CSD2$Chinook.c.well >= LOD["ot"])) * 100,
  pct.above.LOQ = c(mean(CSD2$C.shasta.c.well >= LOQ["cs"]), mean(CSD2$Chinook.c.well >= LOQ["ot"])) * 100,
  max.copies.L  = c(max(y.cs), max(y.ot)),
  mean.copies.L = c(mean(y.cs), mean(y.ot)))
molecular.summary

## Table 1
table1 <- data.frame(model = model.names,
                     transferability.cs = round(selection.cs$transferability[match(model.names, selection.cs$model)], 2),
                     selected.cs = selection.cs$selected[match(model.names, selection.cs$model)],
                     transferability.ot = round(selection.ot$transferability[match(model.names, selection.ot$model)], 2),
                     selected.ot = selection.ot$selected[match(model.names, selection.ot$model)])
table1 <- table1[order(-table1$transferability.cs), ]
write.csv(table1, "../results/Table_1_model_selection.csv", row.names = FALSE)
write.csv(site.summary, "../results/site_summary.csv", row.names = FALSE)
write.csv(molecular.summary, "../results/molecular_summary.csv", row.names = FALSE)
write.csv(n.required, "../results/sample_size_by_site_and_CV.csv", row.names = FALSE)
write.csv(sample.size.supp, "../results/sample_size_by_site_SE_targets.csv", row.names = FALSE)

## Selected model summaries and spatial field parameters
model.summaries <- bind_rows(
  Cs.inla.null$summary.fixed    %>% tibble::rownames_to_column("term") %>% mutate(species = "C. shasta", model = "NULL", type = "fixed"),
  Cs.inla.null$summary.hyperpar %>% tibble::rownames_to_column("term") %>% mutate(species = "C. shasta", model = "NULL", type = "hyperparameter"),
  Chin.inla.d$summary.fixed     %>% tibble::rownames_to_column("term") %>% mutate(species = "O. tshawytscha", model = "Depth", type = "fixed"),
  Chin.inla.d$summary.hyperpar  %>% tibble::rownames_to_column("term") %>% mutate(species = "O. tshawytscha", model = "Depth", type = "hyperparameter"),
  Chin.inla.top$summary.fixed   %>% tibble::rownames_to_column("term") %>% mutate(species = "O. tshawytscha", model = "Distance From Shore + Sample Type", type = "fixed"),
  spat.full.ot$summary.fixed    %>% tibble::rownames_to_column("term") %>% mutate(species = "O. tshawytscha", model = "Depth-Squared + Depth + Distance From Shore", type = "fixed")
) %>% select(species, model, type, term, everything())
write.csv(model.summaries, "../results/model_summaries.csv", row.names = FALSE)
spatial.field <- data.frame(
  species = c("C. shasta", "O. tshawytscha"),
  kappa = c(Kappa.cs, Kappa.ot), range.km = c(Range.cs, Range.ot), sigma.u = c(Sigma_u.cs, Sigma_u.ot),
  spatial.proportion = c(spatial.proportion.cs, spatial.proportion.ot),
  cor.5m = c(matern.at.cs[1], matern.at.ot[1]), cor.20m = c(matern.at.cs[2], matern.at.ot[2]), cor.40m = c(matern.at.cs[3], matern.at.ot[3]),
  pseudo.R2 = c(r2.cs, r2.ot))
write.csv(spatial.field, "../results/spatial_field.csv", row.names = FALSE)

#### 6. manuscript figures ####

## Common figure style
COL_CS <- "#1b7837"; COL_OT <- "#2166ac"
COL_SPECIES <- c("C. shasta" = COL_CS, "O. tshawytscha" = COL_OT)
COL_SITES <- c("1" = "#D55E00", "2" = "#E69F00", "3" = "#009E73", "4" = "#0072B2", "5" = "#CC79A7", "6" = "#56B4E9")
COL_CV   <- c("0.1" = "#1b7837", "0.2" = "#762a83", "0.3" = "#d6604d")
theme_ms <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.6),
          axis.text = element_text(colour = "black", size = base_size - 1),
          axis.title = element_text(face = "bold", size = base_size),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size, hjust = 0),
          legend.position = "bottom", legend.box = "vertical",
          legend.title = element_text(face = "bold", size = base_size - 1),
          legend.text = element_text(size = base_size - 2),
          legend.key.width = grid::unit(1.4, "lines"),
          plot.subtitle = element_text(face = "bold"))
}
theme_set(theme_ms())
sp_lab <- as_labeller(c("C. shasta" = "italic('C. shasta')", "O. tshawytscha" = "italic('O. tshawytscha')"), label_parsed)
FIG_W1 <- 3.5; FIG_W15 <- 5.5; FIG_W2 <- 7.0     # widths in inches: 89, 140 and 178 mm (journal maximum 180 mm)
save_fig <- function(plot, name, w = FIG_W2, h = 5, dir = "../figures") {
  tryCatch({
    force(plot)
    ggsave(file.path(dir, "pdf", paste0(name, ".pdf")), plot, width = w, height = h, device = cairo_pdf, bg = "white")
    ggsave(file.path(dir, paste0(name, ".png")), plot, width = w, height = h, dpi = 600, bg = "white")
    cat("saved", name, "\n")
  }, error = function(e) message("could not make ", name, ": ", conditionMessage(e)))
}

## Long format of the concentrations by species
CSD2$SiteF <- factor(CSD2$Site, levels = 1:6)
long <- bind_rows(
  CSD2 %>% transmute(Site = SiteF, species = "C. shasta",      conc = C.shasta.vol.cor.conc),
  CSD2 %>% transmute(Site = SiteF, species = "O. tshawytscha", conc = Chinook.vol.cor.conc)) %>%
  mutate(species = factor(species, levels = names(COL_SPECIES)))
site.stats <- long %>% group_by(species, Site) %>%
  summarise(cv = sd(conc) / mean(conc), .groups = "drop") %>%
  left_join(site.summary %>% transmute(Site = factor(Site, levels = 1:6), W = width.m + 1), by = "Site")

## Figure 2: eDNA concentration by site
fig2 <- ggplot(long, aes(x = Site, y = conc / 1000, fill = species)) +
  geom_violin(colour = "black", trim = FALSE, alpha = 0.9, linewidth = 0.4) +
  geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 1.8, fill = "white") +
  facet_wrap(~species, ncol = 1, scales = "free_y", labeller = sp_lab) +
  scale_fill_manual(values = COL_SPECIES, guide = "none") +
  labs(x = "Site", y = "eDNA concentration (thousands of copies / L)")
save_fig(fig2, "Figure_2", w = FIG_W15, h = 6)

## Figure 3: covariate distributions and eDNA concentration against each covariate
cov.rows <- list(list(var = "D.M.", lab = "Depth (m)"),
                 list(var = "DistanceFromShore.m", lab = "Distance to shore (m)"),
                 list(var = "D.thal.m", lab = "Distance to thalweg (m)"),
                 list(var = "V.M.", lab = "Velocity (m / s)"))
pt.theme <- theme_ms(11) + theme(plot.margin = margin(4, 8, 2, 4), plot.title.position = "plot",
                                 plot.title = element_text(face = "bold", size = 11, hjust = 0, margin = margin(b = 1)))
ylab.conc <- expression(10^3 ~ "copies / L")
panel.hist <- function(var, lab) ggplot(CSD2, aes(x = .data[[var]])) +
  geom_histogram(bins = 12, colour = "black", fill = "grey80", linewidth = 0.3) + labs(x = lab, y = "Count") + pt.theme
panel.scatter <- function(var, lab, ycol, col) ggplot(CSD2, aes(x = .data[[var]], y = .data[[ycol]] / 1000)) +
  geom_point(colour = col, size = 1.6, alpha = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black", fill = "grey70", linewidth = 0.6, span = 0.9) +
  labs(x = lab, y = ylab.conc) + pt.theme
type.df <- CSD2 %>% mutate(Sampler = factor(ifelse(Type == "G", "Surface", "Depth"), levels = c("Surface", "Depth")))
panel.type.bar <- ggplot(type.df, aes(x = Sampler)) + geom_bar(colour = "black", fill = "grey80", linewidth = 0.3) +
  labs(x = "Sampler type", y = "Count") + pt.theme
panel.type.box <- function(ycol, col) ggplot(type.df, aes(x = Sampler, y = .data[[ycol]] / 1000)) +
  geom_boxplot(fill = col, alpha = 0.6, outlier.size = 1, linewidth = 0.4) + labs(x = "Sampler type", y = ylab.conc) + pt.theme
col1 <- c(lapply(cov.rows, function(r) panel.hist(r$var, r$lab)), list(panel.type.bar))
col2 <- c(lapply(cov.rows, function(r) panel.scatter(r$var, r$lab, "C.shasta.vol.cor.conc", COL_CS)), list(panel.type.box("C.shasta.vol.cor.conc", COL_CS)))
col3 <- c(lapply(cov.rows, function(r) panel.scatter(r$var, r$lab, "Chinook.vol.cor.conc", COL_OT)), list(panel.type.box("Chinook.vol.cor.conc", COL_OT)))
tag <- function(plots, letters) Map(function(p, tg) p + ggtitle(tg), plots, letters)
column <- function(plots, header) annotate_figure(ggarrange(plotlist = plots, ncol = 1, nrow = 5),
                                                  top = text_grob(header, face = "bold.italic", size = 12))
fig3 <- ggarrange(column(tag(col1, paste0(LETTERS[1:5], "1")), "Covariate"),
                  column(tag(col2, paste0(LETTERS[1:5], "2")), "Ceratonova shasta"),
                  column(tag(col3, paste0(LETTERS[1:5], "3")), "Oncorhynchus tshawytscha"), ncol = 3, nrow = 1)
save_fig(fig3, "Figure_3", w = FIG_W2, h = 11.5)

## Figure 4: eDNA concentration on the six cross sections, scaled by site and
## species. Depth samples were collected at 60 percent of the vertical depth,
## so the bed elevation is D.M. / 0.6.
bed <- CSD2 %>% filter(Type01 == 1) %>% transmute(Site = SiteF, x = Distance.Left.Bank.m, y = D.M. / 0.6) %>%
  bind_rows(CSD2 %>% group_by(Site = SiteF) %>% reframe(x = c(0, max(Distance.Left.Bank.m) + 1), y = 0)) %>%
  arrange(Site, x)
xs.panel <- function(conc.col, palette, legend.title, strips) {
  pts <- CSD2 %>% transmute(Site = SiteF, x = Distance.Left.Bank.m, y = D.M. + 0.1, conc = .data[[conc.col]]) %>%
    group_by(Site) %>% mutate(z = as.numeric(scale(conc))) %>% ungroup()
  ggplot() +
    geom_line(data = bed, aes(x = x, y = y), colour = "black", linewidth = 1) +
    geom_point(data = pts, aes(x = x, y = y, fill = z), shape = 21, size = 2.6, colour = "black", stroke = 0.3) +
    scale_y_reverse() +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, palette), name = legend.title,
                         guide = guide_colourbar(barwidth = grid::unit(9, "lines"), title.vjust = 0.9)) +
    facet_wrap(~Site, ncol = 1, scales = "free", strip.position = "right",
               labeller = labeller(Site = function(x) paste("Site", x))) +
    labs(x = NULL, y = NULL) +
    theme_ms(12) + theme(strip.text.y = if (strips) element_text(angle = 0) else element_blank(),
                         legend.position = "top", legend.title = element_text(face = "plain"))
}
fig4 <- ggarrange(xs.panel("C.shasta.vol.cor.conc", "Greens", "Site-scaled copies / L", strips = FALSE) + ggtitle(expression(italic("C. shasta"))),
                  xs.panel("Chinook.vol.cor.conc", "Blues", NULL, strips = TRUE) + ggtitle(expression(italic("O. tshawytscha"))),
                  ncol = 2, nrow = 1, widths = c(1, 1.08))
fig4 <- annotate_figure(fig4, bottom = text_grob("Distance from left bank (m)", face = "bold", size = 12),
                        left = text_grob("Depth (m)", face = "bold", size = 12, rot = 90))
save_fig(fig4, "Figure_4", w = FIG_W2, h = 9.5)

## Figure 5: precision of the section mean against number of samples at each
## site, under random (n_ind) and spatially stratified (n_strat) placement
n.seq <- 1:60
curves5 <- site.stats %>% rowwise() %>% do({
  s <- .; par <- species.par[[as.character(s$species)]]; tb <- matern.tab[[as.character(s$species)]]
  bind_rows(data.frame(species = s$species, Site = s$Site, n = n.seq, CV = cv_random(s$cv, n.seq), design = "Simple random"),
            data.frame(species = s$species, Site = s$Site, n = n.seq, CV = cv_stratified(s$cv, n.seq, s$W, par$p, tb), design = "Spatially stratified"))
}) %>% ungroup() %>%
  mutate(species = factor(species, levels = names(COL_SPECIES)), Site = factor(Site, levels = 1:6),
         design = factor(design, levels = c("Simple random", "Spatially stratified")))
ref.lines <- data.frame(CV = c(0.1, 0.2, 0.3), label = paste("Target CV =", c(0.1, 0.2, 0.3)))
fig5 <- ggplot(curves5, aes(x = n, y = CV, colour = Site, linetype = design)) +
  geom_hline(data = ref.lines, aes(yintercept = CV), colour = "grey50", linetype = "dashed", linewidth = 0.5) +
  geom_text(data = ref.lines, aes(x = 60, y = CV + 0.025, label = label), inherit.aes = FALSE, hjust = 1, size = 3, colour = "grey30") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~species, ncol = 1, scales = "free_y", labeller = sp_lab) +
  scale_colour_manual(values = COL_SITES, name = "Site") +
  scale_linetype_manual(values = c("Simple random" = "solid", "Spatially stratified" = "22"), name = "Sample placement",
                        labels = c(expression("Simple random, independence (" * n[ind] * ")"), expression("Spatially stratified (" * n[strat] * ")"))) +
  scale_x_continuous(limits = c(1, 60), breaks = seq(0, 60, 10), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
  labs(x = "Number of 1-L samples (n)", y = "Precision of the section mean (CV of the estimated mean)") +
  guides(colour = guide_legend(order = 1, nrow = 1), linetype = guide_legend(order = 2))
save_fig(fig5, "Figure_5", w = FIG_W2, h = 7.5)

## Figure 6: fitted versus observed concentrations for the selected models,
## and the multiplicative effects of depth and distance to shore on
## O. tshawytscha from the full spatial model
fit.theme <- theme_ms(12) + theme(plot.title.position = "plot", plot.margin = margin(4, 10, 4, 4),
                                  plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 2)))
p6a <- ggplot(data.frame(fitted = fitted.cs$mean, observed = y.cs), aes(x = fitted / 1000, y = observed / 1000)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40") +
  geom_point(colour = COL_CS, size = 1.8, alpha = 0.8) +
  coord_equal(xlim = c(0, 200), ylim = c(0, 200)) +
  annotate("text", x = 15, y = 185, label = sprintf("R^2 == %.2f", r2.cs), parse = TRUE, hjust = 0, size = 4) +
  labs(x = expression(bold("Fitted (" * 10^3 ~ "copies / L)")), y = expression(bold("Observed (" * 10^3 ~ "copies / L)"))) + fit.theme
p6b <- ggplot(data.frame(fitted = fitted.ot$mean, observed = y.ot), aes(x = log(fitted), y = log(observed))) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40") +
  geom_point(colour = COL_OT, size = 1.8, alpha = 0.8, na.rm = TRUE) +
  coord_equal(xlim = c(2.8, 8.5), ylim = c(2.8, 8.5)) +
  annotate("text", x = 3.1, y = 8.2, label = sprintf("R^2 == %.2f", r2.ot), parse = TRUE, hjust = 0, size = 4) +
  labs(x = "log(fitted copies / L)", y = "log(observed copies / L)") + fit.theme
effect.plot <- function(pred, plotVar, xlab) {
  lc  <- pred$model$summary.lincomb.derived
  idx <- grep(plotVar, pred$newRowNames)
  pd  <- data.frame(x  = pred$newData[idx, grep(plotVar, names(pred$newData), fixed = TRUE)[1]],
                    mn = exp(lc$mean[idx]), lo = exp(lc$`0.025quant`[idx]), hi = exp(lc$`0.975quant`[idx]))
  ggplot(pd, aes(x = x)) +
    geom_hline(yintercept = 1, colour = "grey60", linetype = "dotted") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = COL_OT, alpha = 0.2) +
    geom_line(aes(y = mn), colour = COL_OT, linewidth = 0.9) +
    labs(x = xlab, y = "Multiplicative effect") + fit.theme
}
p6c <- effect.plot(GRSpred, "D.M..z", "Scaled depth")
p6d <- effect.plot(GRSpred, "DistanceFromShore.m.z.z", "Scaled distance to shore")
fig6 <- ggarrange(p6a + ggtitle("A"), p6b + ggtitle("B"), p6c + ggtitle("C"), p6d + ggtitle("D"), ncol = 2, nrow = 2, align = "hv")
save_fig(fig6, "Figure_6", w = FIG_W2, h = 7)

## Figure 7: the general case, n = (sigma / (CV_target x mean))^2, as a lookup
## chart in measured units: the number of samples needed for a target CV of
## 0.1, 0.2, or 0.3, given the mean and standard deviation of eDNA concentration
## from a pilot study. Points are the six cross sections of this study.
mean.g <- exp(seq(log(100), log(200000), length.out = 300))
sd.g   <- exp(seq(log(100), log(200000), length.out = 300))
pts7 <- site.summary %>%
  transmute(Site = factor(Site, levels = 1:6), `C. shasta` = sd.cs, `O. tshawytscha` = sd.ot, mean.cs, mean.ot) %>%
  tidyr::pivot_longer(c(`C. shasta`, `O. tshawytscha`), names_to = "species", values_to = "sd") %>%
  mutate(mean = ifelse(species == "C. shasta", mean.cs, mean.ot), species = factor(species, levels = names(COL_SPECIES)))
lookup.panel <- function(cv_t, title) {
  g <- tidyr::expand_grid(mean = mean.g, sd = sd.g) %>% mutate(n = ceiling((sd / (cv_t * mean))^2))
  g$n[g$n > 500] <- NA
  labs_b <- n_band_labels(n_breaks_default)
  g$band <- cut(g$n, breaks = n_breaks_default, labels = labs_b, right = TRUE, include.lowest = TRUE)
  keys <- data.frame(mean = 1, sd = 1, band = factor(labs_b, levels = labs_b))
  ggplot(g, aes(x = mean, y = sd)) +
    geom_tile(data = keys, aes(fill = band), key_glyph = "rect") +
    geom_raster(data = g[!is.na(g$n), ], aes(fill = band)) +
    geom_raster(data = g[is.na(g$n), ], fill = "grey85") +
    scale_fill_manual(name = "Required 1-L samples (n)", limits = labs_b, drop = FALSE,
                      values = setNames(viridisLite::viridis(length(labs_b), direction = -1), labs_b)) +
    geom_point(data = pts7, aes(x = mean, y = sd, shape = species), fill = "white", colour = "black", size = 2.4, stroke = 0.7) +
    geom_text(data = pts7, aes(x = mean, y = sd, label = Site), colour = "black", size = 2.6, vjust = -0.9, show.legend = FALSE) +
    scale_shape_manual(values = c("C. shasta" = 21, "O. tshawytscha" = 24), name = NULL,
                       labels = c(expression(italic("C. shasta")), expression(italic("O. tshawytscha")))) +
    scale_x_log10(labels = scales::label_comma(), expand = c(0, 0)) +
    scale_y_log10(labels = scales::label_comma(), expand = c(0, 0)) +
    coord_cartesian(xlim = range(mean.g), ylim = range(sd.g)) +
    labs(x = "Mean eDNA concentration (copies / L)", y = expression(bold("Standard deviation, " * sigma * " (copies / L)")), subtitle = title) +
    theme_ms(11) + theme(panel.grid = element_blank(), plot.subtitle = element_text(face = "bold")) +
    guides(fill = guide_legend(order = 1, nrow = 2, title.position = "top"), shape = guide_legend(order = 2, override.aes = list(size = 2.6)))
}
fig7 <- ggarrange(lookup.panel(0.1, "A)  Target CV 0.1 (research)"),
                  lookup.panel(0.2, "B)  Target CV 0.2 (management)"),
                  lookup.panel(0.3, "C)  Target CV 0.3 (pilot work)"),
                  ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
save_fig(fig7, "Figure_7", w = FIG_W15, h = 11)

## Figure 8: samples required against site CV and section width, under random
## placement and under stratified placement at three spatial proportions of
## variance
CONTOUR.RANGE <- 30
tab.c <- matern_tables(CONTOUR.RANGE)
cv.g <- seq(0.10, 1.50, by = 0.01); w.g <- seq(5, 150, by = 1)
grids8 <- c(list(list(g = n_required_grid(cv.g, w.g, 0.2, 0, tab.c, "random", 500), t = "A)  Simple random (any width)")),
            Map(function(p, L) list(g = n_required_grid(cv.g, w.g, 0.2, p, tab.c, "stratified", 500),
                                    t = sprintf("%s)  Stratified, p = %.2f", L, p)),
                c(0.25, 0.5, 0.75), c("B", "C", "D")))
fig8 <- ggarrange(plotlist = lapply(grids8, function(x) contour_panel(x$g, x$t, base_theme = theme_ms(11))),
                  ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
fig8 <- annotate_figure(fig8, top = text_grob(sprintf("Target CV 0.20 for the section mean, correlation range %d m", CONTOUR.RANGE),
                                              face = "bold", size = 12))
save_fig(fig8, "Figure_8", w = FIG_W2, h = 7.6)

#### 7. supporting information figures ####
## Supporting Information figures are numbered S1, S2, ... as required by the
## journal. The corresponding appendix labels of the supplement draft are
## listed in figures/figure_captions.md.

supp <- "../figures/supporting_information"
## Figures S1-S7: model assessment for C. shasta
save_fig(plot(dc.pois.cs) + theme_ms(), "Figure_S1", w = FIG_W15, h = 4, dir = supp)
save_fig(plot(dc.nbn.cs) + theme_ms(),  "Figure_S2", w = FIG_W15, h = 4, dir = supp)
semivariogram <- function(V, ylim, nugget = NULL, sill = NULL, range = NULL) {
  g <- ggplot(V, aes(x = dist, y = gamma)) + geom_point() +
    geom_smooth(span = 0.9, se = FALSE, colour = COL_CS, linewidth = 0.8, method = "loess", formula = y ~ x) +
    coord_cartesian(ylim = ylim) + labs(x = "Distance (km)", y = "Semivariance") + theme_ms()
  if (!is.null(nugget)) g <- g +
    geom_hline(yintercept = nugget, colour = "red") + annotate("text", x = 0.04, y = nugget * 0.95, label = "Nugget", colour = "red") +
    geom_hline(yintercept = sill, colour = "darkgreen") + annotate("text", x = 0.04, y = sill * 1.05, label = "Sill", colour = "darkgreen") +
    geom_vline(xintercept = range, colour = "brown") + annotate("text", x = range * 1.04, y = nugget, label = "Estimated range", colour = "brown", angle = -90)
  g
}
save_fig(semivariogram(as.data.frame(V.nbn.cs), c(0, 30000), nugget.cs, sill.cs, range.cs), "Figure_S3", w = FIG_W15, h = 4.5, dir = supp)
bubble <- function(resid, X, Y) {
  d <- data.frame(resid = resid, X = X, Y = Y, size = 5 * abs(resid) / max(resid), sign = ifelse(resid > 0, "positive", "negative"))
  ggplot(d, aes(x = X, y = Y, size = size, colour = sign)) + geom_point(shape = 1) +
    scale_colour_manual(values = c(positive = "red", negative = "black"), name = "Residual") + scale_size_identity() +
    labs(x = "Downstream distance (km / 50)", y = "Distance from left bank (km)") + theme_ms()
}
save_fig(bubble(resid.nbn.cs, utm$Xkm / 50, utm$Ykm), "Figure_S4", w = FIG_W2, h = 4, dir = supp)
save_fig(semivariogram(as.data.frame(V.full.cs), c(0, 12500)), "Figure_S5", w = FIG_W15, h = 4.5, dir = supp)
save_fig(ggplot_inla_residuals2(Cs.inla.null, y.cs, CI = TRUE, method = NA) + theme_ms() + theme(plot.margin = margin(5.5, 16, 5.5, 5.5)), "Figure_S6", w = FIG_W15, h = 4.5, dir = supp)
save_fig(plot(ds.null.cs) + theme_ms(), "Figure_S7", w = FIG_W15, h = 4.5, dir = supp)

## Figures S8-S11: posterior densities of range and kappa for the selected models
posterior.plot <- function(marg, xlab) ggplot(as.data.frame(marg), aes(x = x, y = y)) + geom_line(linewidth = 0.8) +
  labs(x = xlab, y = "Posterior density") + theme_ms()
save_fig(posterior.plot(SpatField.cs$marginals.range.nominal[[1]], "Nominal range (km)"), "Figure_S8", w = FIG_W1 * 1.3, h = 3.5, dir = supp)
save_fig(posterior.plot(SpatField.cs$marginals.kappa[[1]], expression(kappa)), "Figure_S9", w = FIG_W1 * 1.3, h = 3.5, dir = supp)
save_fig(posterior.plot(SpatField.ot$marginals.range.nominal[[1]], "Nominal range (km)"), "Figure_S10", w = FIG_W1 * 1.3, h = 3.5, dir = supp)
save_fig(posterior.plot(SpatField.ot$marginals.kappa[[1]], expression(kappa)), "Figure_S11", w = FIG_W1 * 1.3, h = 3.5, dir = supp)

## Figure S12: difference between paired surface and depth samples of
## O. tshawytscha at each vertical
bardata <- CSD2 %>% group_by(Distance.Left.Bank.m, Site, Type01) %>% summarise(Chin = mean(Chinook.vol.cor.conc), .groups = "drop_last") %>%
  summarise(diff = Chin[1] - Chin[2], .groups = "drop") %>% drop_na() %>% mutate(Site = factor(Site, levels = 1:6))
pct.surface <- bardata %>% group_by(Site) %>% summarise(label = paste0(round(mean(diff > 0) * 100, 1), "% surface > depth"), .groups = "drop")
figD1 <- ggplot(bardata, aes(x = Distance.Left.Bank.m, y = diff / 1000)) +
  geom_hline(yintercept = 0, linewidth = 0.4) + geom_col(fill = COL_OT, width = 2) +
  geom_text(data = pct.surface, aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.5, size = 3) +
  facet_wrap(~Site, ncol = 2, scales = "free", labeller = labeller(Site = function(x) paste("Site", x))) +
  labs(x = "Distance from left bank (m)", y = expression(bold("Surface" ~ "−" ~ "depth (" * 10^3 ~ "copies / L)"))) + theme_ms()
save_fig(figD1, "Figure_S12", w = FIG_W2, h = 7, dir = supp)

## Figure S13: C. shasta concentration on cross sections 1 and 3 with a tensor
## spline surface (mgcv)
xs.spline <- function(site.n) {
  d <- subset(CSD2, Site == site.n)
  g <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M., Distance.Left.Bank.m, k = 3), data = d)
  new <- expand.grid(Distance.Left.Bank.m = seq(0, max(d$Distance.Left.Bank.m), 0.1), D.M. = seq(0, max(d$D.M.) / 0.6, 0.02))
  new$z <- predict(g, newdata = new)
  bed <- d %>% filter(Type01 == 1) %>% transmute(x = Distance.Left.Bank.m, y = D.M. / 0.6) %>%
    bind_rows(data.frame(x = c(0, max(d$Distance.Left.Bank.m) + 1), y = 0)) %>% arrange(x)
  ggplot() +
    geom_tile(data = new, aes(x = Distance.Left.Bank.m, y = D.M., fill = z)) +
    geom_contour(data = new, aes(x = Distance.Left.Bank.m, y = D.M., z = z), colour = "grey30", linewidth = 0.3) +
    geom_point(data = d, aes(x = Distance.Left.Bank.m, y = D.M., fill = as.numeric(scale(C.shasta.vol.cor.conc))), shape = 21, size = 3) +
    geom_line(data = bed, aes(x = x, y = y), linewidth = 1) +
    scale_y_reverse() + scale_fill_viridis_c(name = "Site-scaled concentration (z-score)") +
    labs(x = "Distance from left bank (m)", y = "Depth (m)", subtitle = paste("Site", site.n)) + theme_ms(11)
}
figD2 <- ggarrange(xs.spline(1), xs.spline(3), ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
save_fig(figD2, "Figure_S13", w = FIG_W2, h = 7, dir = supp)

#### 8. save ####

dev.off()
save(CSD2, utm, site.summary, site.stats, long, scale.par,
     selection.cs, selection.ot, table1, molecular.summary, n.required, sample.size.supp,
     fitted.cs, fitted.ot, r2.cs, r2.ot, r2.null.ot, GRSpred,
     Kappa.cs, Range.cs, Sigma_u.cs, cor.cs, matern.at.cs, SpatField.cs,
     Kappa.ot, Range.ot, Sigma_u.ot, cor.ot, matern.at.ot, SpatField.ot,
     dc.pois.cs, dc.nbn.cs, dc.null.cs, ds.null.cs, V.nbn.cs, V.full.cs, resid.nbn.cs,
     dc.pois.ot, dc.nbn.ot, dc.d.ot, ds.d.ot, V.nbn.ot, V.full.ot, resid.nbn.ot,
     nugget.cs, sill.cs, range.cs, nugget.ot, sill.ot, range.ot, spatial.proportion.cs, spatial.proportion.ot,
     species.par, model.summaries, spatial.field,
     file = "../results/analysis_cache.rda")
writeLines(capture.output(sessionInfo()), "../results/sessionInfo.txt")
cat("done\n")
