options(repos = c(CRAN = "https://cloud.r-project.org"))

# Use user library so we don't need admin rights
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))
cat("Installing to:", user_lib, "\n")

# Standard CRAN packages
cran_pkgs <- c(
  'doSNOW', 'doParallel', 'foreach',
  'tidyverse', 'tictoc', 'fields', 'ggplot2',
  'sp', 'lattice', 'plotly', 'geometry', 'viridis',
  'gstat', 'remotes', 'kableExtra', 'splancs',
  'sf', 'terra'  # modern replacements that support rgdal functions
)

installed <- rownames(installed.packages())
to_install <- cran_pkgs[!cran_pkgs %in% installed]
if (length(to_install) > 0) {
  cat("Installing CRAN packages:", paste(to_install, collapse=", "), "\n")
  install.packages(to_install, dependencies = TRUE)
} else {
  cat("All CRAN packages already installed.\n")
}

# doMPI requires MPI — skip if MPI not available
tryCatch(install.packages("doMPI"), error = function(e) cat("doMPI skipped:", e$message, "\n"))

# INLA — special repo
if (!"INLA" %in% rownames(installed.packages())) {
  cat("Installing INLA...\n")
  install.packages("INLA",
    repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
    dependencies = TRUE)
}

# rgdal — retired from CRAN; install from archive
if (!"rgdal" %in% rownames(installed.packages())) {
  cat("Installing rgdal from archive...\n")
  tryCatch(
    remotes::install_version("rgdal", version = "1.6-7",
      repos = "https://cloud.r-project.org"),
    error = function(e) cat("rgdal archive install failed:", e$message, "\n")
  )
}

# GitHub packages
library(remotes)
gh_pkgs <- list(
  INLAOutputs  = "oswaldosantos/INLAOutputs",
  inlatools    = "inbo/inlatools",
  INLAutils    = "timcdlucas/INLAutils",
  ggregplot    = "gfalbery/ggregplot",
  inlamesh3d   = "inlabru-org/inlamesh3d"
)

for (pkg in names(gh_pkgs)) {
  if (!pkg %in% rownames(installed.packages())) {
    cat("Installing", pkg, "from GitHub...\n")
    tryCatch(remotes::install_github(gh_pkgs[[pkg]]),
             error = function(e) cat(pkg, "failed:", e$message, "\n"))
  }
}

cat("\n--- Final status ---\n")
all_pkgs <- c(cran_pkgs, 'doMPI', 'INLA', 'rgdal',
              'INLAOutputs', 'inlamesh3d', 'inlatools', 'INLAutils', 'ggregplot')
installed2 <- rownames(installed.packages())
cat("Installed:", paste(all_pkgs[all_pkgs %in% installed2], collapse=", "), "\n")
cat("Still missing:", paste(all_pkgs[!all_pkgs %in% installed2], collapse=", "), "\n")
