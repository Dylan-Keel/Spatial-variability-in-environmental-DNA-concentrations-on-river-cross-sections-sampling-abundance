# Instructions

Run `Main_Analysis_and_Figures.R` from this folder to replicate the analyses and create the figures and tables presented in the manuscript and its supplement. The script is in eight sections: setup, data, the *C. shasta* GLMMs, the *O. tshawytscha* GLMMs, estimating the mean eDNA concentration on a cross section, results tables, manuscript figures, and supplemental figures. Each section is commented.

`Sampling_Design_Functions.R` holds the sample size calculations for random and spatially stratified placement of samples on a cross section (Figures 5, 7, and 8). `Sampling_Design_App.R` is a step-by-step tool built on the same calculations: the user sets a target CV, the site CV from a pilot study, the width and spatial structure of the cross section, and how samples are placed, and it returns the number of 1-L samples needed, the precision reached by a given number, a simulated cross section, and planning maps of n against width and site CV. Run it with `shiny::runApp("Sampling_Design_App.R")`, or open `Sampling_Design_App.html` in a browser for a version that needs no R; the two implement the same calculations.

`HighstatLibV13.R` (Zuur et al. 2009), `INLA_plotting_functions.R`, and `inlatools_modified.R` (adapted from the inlatools package) are sourced by the main script for data exploration, model predictions, and model checks.

To clone this repository to your local machine, run the following command in your terminal or Git Bash:

git clone https://github.com/Dylan-Keel/Spatial-variability-in-environmental-DNA-concentrations-on-river-cross-sections-sampling-abundance.git

Install the packages listed at the top of the main script before running it.
