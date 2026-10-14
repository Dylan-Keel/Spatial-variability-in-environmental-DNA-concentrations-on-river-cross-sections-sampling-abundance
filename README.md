# Spatial variability in environmental DNA concentrations on river cross sections: implications for sampling design and abundance estimation

This repository holds the data, code, and metadata for the manuscript above (Keel, Som, Henderson, and Kinziger; prepared for Molecular Ecology Resources). It will be the home for resources related to the publication.

# Project Abstract

High variance across repeated measurements of environmental DNA (eDNA) concentrations may be partially responsible for the mixed efficacy of eDNA concentration as a proxy for aquatic organism abundance or biomass in streams. While environmental factors such as temperature and stream discharge are commonly used to explain variability in eDNA concentrations, few studies have assessed variability of eDNA concentration across river cross sections and rarely do eDNA studies incorporate spatially stratified sampling to account for variability within a site. We used droplet digital PCR to measure the concentration of extra-organismal eDNA from Chinook salmon (*Oncorhynchus tshawytscha*) and the organismal eDNA from their abundant and pathogenic myxozoan parasite, *Ceratonova shasta*, in water samples collected from six cross sections of a river (basin area ≅ 40,000 km²). Concentrations of eDNA were heterogeneously distributed on river cross sections varying from non-detectable levels to 257,222 copies per liter. Within cross sections the coefficient of variation (CV) of eDNA concentration spanned from 0.40 to 1.43. We found that to estimate the mean eDNA concentration at a river cross section with a level of precision (target CV = 0.2), the number of spatially stratified water samples varied by species and cross section. For *C. shasta*, an average of 15 one-liter samples were necessary, while *O. tshawytscha* required an average of 11 one-liter samples. A generalized linear mixed effects model fit to the *C. shasta* data did not explain a meaningful amount of variability in eDNA concentration, but the model fit to the *O. tshawytscha* data identified a pattern where eDNA concentrations decreased with sample depth.

# What is here

```
R_code/     Main_Analysis_and_Figures.R   the complete analysis, every figure, and the results tables
            Sampling_Design_Functions.R   sample size under random and spatially stratified placement
            Sampling_Design_App.R         interactive version of the sampling design figures (Shiny)
            Sampling_Design_App.html      the same, as a stand-alone web page
            HighstatLibV13.R, INLA_plotting_functions.R, inlatools_modified.R   helper code
data/       the 192 water samples with ddPCR results (CSD.csv), field data, and the saved
            cross-validation results
figures/    Figures 1 to 8 as 600 dpi PNG (vector PDF copies in figures/pdf/), the Supporting
            Information figures S1 to S13 in figures/supporting_information/, and the captions
results/    results tables, model summaries, Moran's I of the residuals, the analysis cache, and the session information
```

# Installation

Analysis was conducted with R version 4.5.3 and R-INLA (see `results/sessionInfo.txt` for every package version). To rerun it, clone the repository, open R in the `R_code` folder, install the packages listed at the top of `Main_Analysis_and_Figures.R` (INLA comes from https://inla.r-inla-download.org/R/stable; INLAutils and inlatools come from GitHub), and run the script. It takes 10 to 15 minutes. Each INLA fit runs on a single thread because INLA's multithreaded fits give slightly different posteriors on every run. On one thread the fits reproduce the values in `results/` to the precision reported in the manuscript; the full *O. tshawytscha* spatial model can settle on one of two optima that differ in the fourth significant digit. The figure section can be rerun on its own from `results/analysis_cache.rda`.

The model selection in the manuscript used 100 repeats of 10-fold cross-validation over sixteen candidate models, 16,000 spatial fits per species. Each fit takes about 90 seconds on one thread, so the full cross-validation is roughly 800 CPU-hours. The saved results are in `data/` and are what the script reads. Set `RUN.KFOLD <- TRUE` at the top of the script to repeat it: it runs in parallel on all but two of the machine's cores, writes one file per repeat to `results/cross_validation/`, and can be stopped and resumed.

**Sampling design tool.** An interactive version of Figures 5, 7, and 8 is at
https://htmlpreview.github.io/?https://raw.githubusercontent.com/Dylan-Keel/Spatial-variability-in-environmental-DNA-concentrations-on-river-cross-sections-sampling-abundance/main/R_code/Sampling_Design_App.html
(or download `R_code/Sampling_Design_App.html` and open it in a browser; the R version runs with `shiny::runApp("Sampling_Design_App.R")` from the `R_code` folder). It is meant for planning: given a pilot estimate of the site CV, the coefficient of variation among single 1-L samples on a cross section, together with its width and the spatial structure of concentration, it gives the number of 1-L samples needed to estimate the mean eDNA concentration of that cross section to a target CV, under random and spatially stratified placement, and the precision reached by a given n. The six Klamath cross sections of the study can be loaded to see what the 2021 sampling achieved. The tool addresses the mean of one cross section on one day; it says nothing about how many cross sections are needed to represent a reach, about detection limits, or about day-to-day change.

The analysis was first run in 2023 (R 4.4.2). The code was reorganized and rerun on R 4.5.3 in 2026 after three packages it used were retired from CRAN; the results in `results/` are from that rerun and reproduce the original values.

Metadata for the data files are in `data/metadata/`: one FGDC (CSDGM) record per file, with a definition, units, and range or code list for every column, validated against the FGDC schema. They open in the USGS MetadataWizard (https://github.com/DOI-USGS/fort-pymdwizard) and accompany the Dryad data release.

# Contribute

Please contact the corresponding authors of the paper for information about this repository or about collaborating on future work.

# Code References

Cochran, W. G. (1977). Sampling techniques (3rd ed.). John Wiley & Sons.

Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B, 73(4), 423–498.

Onkelinx, T. (2024). inlatools: diagnostic tools for INLA models. R package. https://github.com/inbo/inlatools

Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations (with discussion). Journal of the Royal Statistical Society, Series B, 71(2), 319–392.

Zuur, A. F., Ieno, E. N., Walker, N., Saveliev, A. A., & Smith, G. M. (2009). Mixed effects models and extensions in ecology with R. Springer New York. https://doi.org/10.1007/978-0-387-87458-6
