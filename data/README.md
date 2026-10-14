# Data

`CSD.csv` is the analysis dataset: one row for each of the 192 water samples, with the sampling location on the cross section (distance from the left bank, depth, distance to shore, distance to thalweg, velocity, sampler type), the water volume filtered and the elution volume, and the ddPCR results for *C. shasta* and *O. tshawytscha* (Poisson-corrected copies per reaction with their limits). The main script converts these to copies per liter.

`2021FIELD_DATA_Klamath.2.csv` and `2021FIELD_DATA_Klamath_Discharge.csv` are the field sheets from the six cross sections, including the depth-integrated samples and the discharge measurements. `MASTER.FINAL.ddPCRdata.Cshasta2021.csv` is the complete ddPCR output, including the no-template controls and field blanks.

`C.shasta.metar.new.rda` and `Chinook.metar.new.rda` are the saved results of the 100 repeats of 10-fold cross-validation used for model selection (validation R² for each of the sixteen candidate models, in each fold and repeat).

Column-level metadata, compiled with the USGS MetadataWizard, are distributed with the data release.
