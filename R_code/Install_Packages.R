## Installs the packages used by Main_Analysis_and_Figures.R.
## INLA comes from its own repository and two packages come from GitHub.

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c("dplyr", "tidyr", "purrr", "stringr", "ggplot2", "ggpubr", "sp", "gstat",
          "mgcv", "assertthat", "tibble", "remotes", "shiny", "scales")
missing <- cran[!cran %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, dependencies = TRUE)

if (!"INLA" %in% rownames(installed.packages())) {
  install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
                   dependencies = TRUE)
}

github <- c(inlatools = "inbo/inlatools", INLAutils = "timcdlucas/INLAutils")
for (pkg in names(github)) {
  if (!pkg %in% rownames(installed.packages())) remotes::install_github(github[[pkg]])
}
