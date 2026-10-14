# Figures

Files are named as the journal asks (Figure_1, Figure_2, ...). Each figure is in `figures/` as a 600 dpi PNG and in `figures/pdf/` as a vector PDF; widths are 89, 140 or 178 mm. Figure 1 is assembled by hand (map, sampling design sketch and photographs) and is not produced by the code.

**Figure 1.** A) Klamath River Basin, showing sampling sites, major tributaries, the study area, and the former Iron Gate Dam. B) Cross section sampling design for collection of water samples, where *b*ᵢ is the distance to shore to the center of the sampling vertical (gridded rectangles), *d*ᵢ is the depth of the sampling vertical, black stars are surface sample locations, and gray stars are depth sample locations. C) Sounding instruments, and cataraft mounted crane and winch for the deployment of dual US DH-59 isokinetic samplers. D) US DH-59 isokinetic sampler and sample jar used for collection of samples at depth.

**Figure 2.** Violin plots of eDNA concentration (thousands of copies per liter) by site for *C. shasta* (top) and *O. tshawytscha* (bottom). White dots are the median and boxes are the interquartile range. Note that the y-axis scales differ between species.

**Figure 3.** Frequency distributions of depth (m), distance to shore (m), distance to thalweg (m), water velocity (m/s), and sampler type (column 1), and eDNA concentration against each covariate for *C. shasta* (column 2) and *O. tshawytscha* (column 3). Lines are locally estimated scatterplot smoothers (LOESS) with 95 percent confidence bands and make no assumptions about underlying structure.

**Figure 4.** Environmental DNA concentrations of *C. shasta* (left) and *O. tshawytscha* (right) at the six cross sections. Concentrations are scaled by site and species to a mean of zero and a standard deviation of one. The thick black line is the measured river bottom.

**Figure 5.** Precision of the estimated mean eDNA concentration on each of the six cross sections (CV of the mean) for 1 to 60 one-liter samples, from the observed variability at each site. Solid lines follow *n*_ind = (CV_site / CV_target)² and apply to samples placed at random across the section (independence). Dashed lines (*n*_strat) apply to samples spaced evenly across the section (spatially stratified), using the spatial proportion of variance and correlation range estimated for each species and the measured width of each section. Horizontal lines mark target CVs of the estimated mean of 0.1, 0.2, and 0.3.

**Figure 6.** A) Fitted versus observed *C. shasta* eDNA concentrations (thousands of copies per liter) with a 1:1 line. B) Fitted versus observed *O. tshawytscha* eDNA concentrations (natural log of copies per liter) with a 1:1 line. C) The multiplicative effect of scaled depth (mean of zero, standard deviation of one) on *O. tshawytscha* eDNA concentration. D) The multiplicative effect of scaled distance to shore on *O. tshawytscha* eDNA concentration. Bands in C and D are 95 percent credible intervals from the full spatial model.

**Figure 7.** The general case, *n*_ind = (CV_site / CV_target)², where the number of one-liter samples needed on a cross section (*n*) depends on the coefficient of variation of eDNA concentration among samples on the section (site CV) and the precision wanted for the estimated mean (target CV). The three target CVs represent values suitable for research (0.1), management (0.2), and preliminary or pilot work (0.3). Points are the six cross sections of this study at a target CV of 0.2 (circles, *C. shasta*; triangles, *O. tshawytscha*).

**Figure 8.** Number of one-liter samples needed to estimate the mean eDNA concentration on a cross section to a target CV of 0.2, as a function of the site CV and the width of the cross section. A) Samples placed at random, for which the width of the section does not matter. B to D) Samples spaced evenly across the section, when 25, 50, or 75 percent of the variance in concentration is spatially structured (*p*), with a Matérn correlation range of 30 m. Spreading samples over a wider section, relative to the correlation range, makes them less correlated and fewer are needed.

# Supporting Information

Molecular Ecology publishes the supplement as Supporting Information with figures numbered S1, S2, ... Files are in `figures/supporting_information/` (PNG) and `figures/supporting_information/pdf/`. The letters in brackets are the appendix numbers used in the current supplement draft.

**Figure S1 (A1).** Simulated probability density of overdispersion from 1,000 samples of the full non-spatial Poisson GLMM for *C. shasta*. The red line is the dispersion of the data.

**Figure S2 (A2).** As S1, for the full non-spatial negative binomial GLMM.

**Figure S3 (A3).** Semivariance of the Pearson residuals of the non-spatial negative binomial GLMM for *C. shasta* over distance (km, with the downstream axis shrunk by a factor of 50). Sill 24,200, nugget 9,600, and an estimated range of 50 m; 60.3 percent of the residual variance is spatial.

**Figure S4 (A4).** Spatial Pearson residuals of the non-spatial negative binomial GLMM for *C. shasta* plotted in pseudo-space. Red circles are positive residuals, black circles negative, and circle size is proportional to the residual.

**Figure S5 (A5).** Semivariance of the Pearson residuals of the full spatial negative binomial GLMM for *C. shasta*. No sill is evident, so residual spatial autocorrelation is taken to be negligible.

**Figure S6 (A6).** Standardized residuals of the selected *C. shasta* model against fitted values.

**Figure S7 (A7).** Posterior predictive check of the selected *C. shasta* model: the ratio of the observed to the expected empirical cumulative distribution of counts, with its 95 percent interval, from 1,000 posterior simulations.

**Figure S8 and S9 (C1, C2).** Posterior densities of the nominal range (km) and κ of the spatial field of the selected *C. shasta* model.

**Figure S10 and S11 (C3, C4).** Posterior densities of the nominal range (km) and κ of the spatial field of the selected *O. tshawytscha* model.

**Figure S12 (D, "Figure 9" in the draft).** Difference in *O. tshawytscha* eDNA concentration (thousands of copies per liter) between paired surface and depth samples at each vertical of the six cross sections, with the percentage of pairs at each site in which the surface sample was higher.

**Figure S13 (D, "Figure 11" in the draft).** *C. shasta* eDNA concentration on cross sections one and three, with a two-dimensional tensor spline surface (R package mgcv) showing the cline at site three. Points are site-scaled measured concentrations; the thick black line is the river bottom.

**Figure S14 (new).** Figure 5 under a different target of inference: the mean concentration of the reach rather than of the sampled cross section. Solid lines are as in Figure 5. Dotted lines (*n*_adj, the autocorrelation-adjusted sample size) treat the samples on a section as one cluster drawn from the reach, CV(*n*) = CV_site √[(1 + (*n* − 1)ρ) / (*n*(1 − ρ))], where ρ is the intraclass correlation among samples on the section. The deviation shared by all samples on a section cannot be removed by adding samples to it, so precision levels off; only additional sections reduce it.

Appendix B of the draft (limit of detection and quantification, Figures B1 to B6) comes from a separate analysis of the standard curves and is not produced by this code.
