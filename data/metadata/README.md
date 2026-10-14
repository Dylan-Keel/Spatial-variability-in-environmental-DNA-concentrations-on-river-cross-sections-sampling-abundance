# Metadata for the data files

Readable version of the FGDC (CSDGM) metadata records in this folder. The XML files are the formal records, submitted with the Dryad data release and readable with the USGS MetadataWizard; `Metadata.xlsx` holds the same content as a workbook. All three are generated from the same source.

## Dataset

| Item | Value |
|---|---|
| Title | Spatial variability in environmental DNA concentrations on river cross sections: implications for sampling design and abundance estimation |
| Authors | Keel, Dylan J.; Som, Nicholas A.; Henderson, Mark J.; Kinziger, Andrew P. |
| Abstract | Environmental DNA (eDNA) concentrations of the myxozoan parasite Ceratonova shasta (organismal eDNA, actinospores) and of Chinook salmon Oncorhynchus tshawytscha (extra-organismal eDNA) were measured by droplet digital PCR in 192 one-liter water samples collected on six cross sections of the Klamath River, California, 16 to 22 April 2021. On each cross section 15 to 19 verticals were sampled at the surface (grab) and at 60 percent of depth (US DH-59 isokinetic sampler), with water velocity, depth, and discharge measured at each vertical. The data support an analysis of within-cross-section variability of eDNA concentration, spatial generalized linear mixed models, and sample size calculations for estimating the mean eDNA concentration of a river cross section. |
| Purpose | To quantify spatial variability of eDNA concentration across river cross sections and to inform the number and placement of water samples needed to estimate the mean eDNA concentration of a cross section to a chosen precision. |
| Time period | 16 to 22 April 2021 (ground condition) |
| Place | Klamath River, Siskiyou County, California, USA; six cross sections between the Klamathon Bridge and Seiad Valley |
| Bounding coordinates | West -123.25, East -122.40, North 41.98, South 41.78 (approximate) |
| Sites | KLM = site 1, just downstream of the Klamathon Bridge; I5 = site 2, just downstream of the Interstate 5 bridge; TOH = site 3, just downstream of the Tree of Heaven Campground; BVR = site 4, just upstream of the Beaver Creek confluence; KMN = site 5, just downstream of the Kinsman fish trapping site; SV = site 6, in Seiad Valley just upstream of the Grider Creek confluence |
| Attribute accuracy | ddPCR reactions with fewer than 10,000 accepted droplets or without clear separation of positive and negative droplets were rerun (less than 10 percent of reactions). A random subset of extracts (25 for C. shasta, 23 for O. tshawytscha) was run twice to assess measurement error. Limits of detection were 6.3 (C. shasta) and 11.8 (O. tshawytscha) copies per reaction; samples below the limit of detection are recorded as zero in the analysis. No-template controls and field blanks were below the limit of detection. |
| Completeness | All 192 surface and depth water samples from the six cross sections are present in CSD.csv. The field sheet also lists depth-integrated samples and field blanks that were not used in the published analysis. |
| Field methods | Water samples (target 1 L) were collected 16 to 22 April 2021 at 15 to 19 verticals on each of six cross sections of the Klamath River, one at the surface (Whirl-Pak bag) and one at 60 percent of depth (US DH-59 isokinetic sampler). Velocity (Price AA meter at 20 and 60 percent of depth), depth (sounding reel), and distance from the bank (tape on a static line) were recorded at each vertical; discharge was computed by the midsection method. |
| Laboratory methods | Samples were filtered in the field within 10 hours (47 mm, 5 um cellulose nitrate). DNA was extracted by acetone dissolution of the filter followed by the Qiagen DNeasy Blood and Tissue protocol with a 100 uL final elution. Each 22 uL ddPCR reaction (Bio-Rad QX200) contained 15 uL of extract and species-specific primers and probes for C. shasta (HEX) and O. tshawytscha (FAM); Poisson-corrected copies per reaction were estimated by Bio-Rad QX Manager 1.2. |
| Concentration | Copies per liter of water were computed as copies per reaction x elution volume (uL) / 15 uL x volume filtered (mL) / 1000, in R_code/Main_Analysis_and_Figures.R of the repository. |
| Contact | Dylan J. Keel, Resource Environmental Solutions, LLC, Sacramento, CA 95818; +1 (831) 227-1168; dkeel@res.us |
| Code and data | https://github.com/Dylan-Keel/Spatial-variability-in-environmental-DNA-concentrations-on-river-cross-sections-sampling-abundance |
| Data release | Dryad DOI to be assigned |
| Metadata standard | FGDC Content Standard for Digital Geospatial Metadata (FGDC-STD-001-1998); full records are the XML files in this folder |
| Metadata date | 20261014 |

## CSD.csv

Analysis dataset: one row per one-liter water sample (192 rows), with its position on the cross section, hydraulics, volume filtered, elution volume, and ddPCR results for C. shasta and O. tshawytscha.

| Column | Definition | Units | Range or codes |
|---|---|---|---|
| `Sample.ID` | Water sample identifier: site code, sample type (G, surface grab; P, point sample at depth), and vertical number. |  | Unique text identifier |
| `Chinook.c.ul` | O. tshawytscha target DNA concentration in the ddPCR reaction, Poisson-corrected, reported by Bio-Rad QX Manager. | copies per microliter of reaction | 0 to 27.4 |
| `Chinook.c.uL.lower` | Lower bound of the 95 percent Poisson confidence interval of Chinook.c.ul. | copies per microliter of reaction | 0 to 24.9 |
| `Chinook.c.uL.upper` | Upper bound of the 95 percent Poisson confidence interval of Chinook.c.ul. | copies per microliter of reaction | 0.255 to 30 |
| `Chinook.c.well` | O. tshawytscha target DNA copies per 20 microliter ddPCR reaction (Chinook.c.ul times 20). | copies per reaction | 0 to 549 |
| `Chinook.c.well.lower` | Lower bound of the 95 percent Poisson confidence interval of Chinook.c.well. | copies per reaction | 0 to 498 |
| `Chinook.c.well.upper` | Upper bound of the 95 percent Poisson confidence interval of Chinook.c.well. | copies per reaction | 5.1 to 600 |
| `C.shasta.c.uL` | C. shasta target DNA concentration in the ddPCR reaction, Poisson-corrected, reported by Bio-Rad QX Manager. | copies per microliter of reaction | 0 to 1837 |
| `C.shasta.c.uL.lower` | Lower bound of the 95 percent Poisson confidence interval of C.shasta.c.uL. | copies per microliter of reaction | 0 to 1798 |
| `C.shasta.c.uL.upper` | Upper bound of the 95 percent Poisson confidence interval of C.shasta.c.uL. | copies per microliter of reaction | 0.178 to 1878 |
| `C.shasta.c.well` | C. shasta target DNA copies per 20 microliter ddPCR reaction (C.shasta.c.uL times 20). | copies per reaction | 0 to 36746 |
| `C.shasta.c.well.lower` | Lower bound of the 95 percent Poisson confidence interval of C.shasta.c.well. | copies per reaction | 0 to 35960 |
| `C.shasta.c.well.upper` | Upper bound of the 95 percent Poisson confidence interval of C.shasta.c.well. | copies per reaction | 3.56 to 37560 |
| `Accepted.Droplets` | Number of droplets read by the Bio-Rad QX200 droplet reader for the reaction; reactions with fewer than 10,000 droplets were rerun. | droplets | 8993 to 21002 |
| `Velocity.0.8` | Water velocity at 80 percent of total depth on the vertical, Price AA current meter. | feet per second | 0 to 4.4274 |
| `Velocity.0.6` | Water velocity at 60 percent of total depth on the vertical, Price AA current meter. | feet per second | 0 to 5.82377 |
| `Velocity.0.2` | Water velocity at 20 percent of total depth on the vertical, Price AA current meter. | feet per second | 0 to 6.92617 |
| `Point.Depth` | Depth at which the water sample was collected (0 for surface grabs). | feet | 0 to 5.4 |
| `Volume..ml.` | Volume of water filtered for the sample. | milliliters | 450 to 1300 |
| `Distance.Left.Bank.ft.` | Distance of the sampling vertical from the left bank (looking downstream). | feet | 0 to 144.357 |
| `Distance.From.Shore` | Distance of the sampling vertical from the nearer bank. | feet | 0 to 78.7402 |
| `Site` | Cross section (site) code. |  | KLM = Site 1, just downstream of the Klamathon Bridge; I5 = Site 2, just downstream of the Interstate 5 bridge; TOH = Site 3, just downstream of the Tree of Heaven Campground; BVR = Site 4, just upstream of the Beaver Creek confluence; KMN = Site 5, just downstream of the Kinsman fish trapping site; SV = Site 6, in Seiad Valley just upstream of the Grider Creek confluence (river left) |
| `Type` | Sample type. |  | G = Surface grab sample collected in a sterile Whirl-Pak bag upstream of the raft; P = Point sample collected at 60 percent of total depth with a US DH-59 isokinetic sampler |
| `Q.cfs` | Total discharge of the cross section on the sampling day, midsection method. | cubic feet per second | 1238.74 to 2358.01 |
| `Q.L.s` | Total discharge of the cross section on the sampling day. | liters per second | 35077 to 66771.3 |
| `Q.` | Share of the total discharge carried by the subsection around the vertical. | percent | 0 to 13.9789 |
| `SubQ..L.s.` | Discharge of the subsection around the vertical. | liters per second | 0 to 6285.6 |
| `JulDay` | Day of year of sampling in 2021 (106 = 16 April, 112 = 22 April). | day of year | 106 to 112 |
| `Elution.volum..uL.` | Volume of buffer in which the extracted DNA was eluted. | microliters | 100 to 150 |
| `DistanceFromShore.m` | Distance of the sampling vertical from the nearer bank. | meters | 0 to 24 |
| `Distance.Left.Bank.m` | Distance of the sampling vertical from the left bank (looking downstream). | meters | 0 to 44 |
| `Width.Interval` | Width of the discharge subsection represented by the vertical. | meters | 0 to 8 |
| `Width.Interval.Round.m` | Width of the discharge subsection represented by the vertical, rounded. | meters | 0 to 8 |
| `Chinook.vol.cor.conc` | O. tshawytscha eDNA concentration in the water sample, rounded to the nearest copy: Chinook.c.well x Elution.volum..uL. / 15 x Volume..ml. / 1000. | copies per liter | 0 to 3733 |
| `C.shasta.vol.cor.conc` | C. shasta eDNA concentration in the water sample, rounded to the nearest copy: C.shasta.c.well x Elution.volum..uL. / 15 x Volume..ml. / 1000. | copies per liter | 0 to 257222 |
| `SiteN` | Site number, 1 (upstream) to 6 (downstream). | site number | 1 to 6 |
| `Meters` | River distance of the cross section from the mouth of the Klamath River (river kilometer x 1000), one value per site, used as the along-river coordinate of the spatial model. | meters | 212602 to 300533 |
| `V.M.` | Water velocity at 60 percent of total depth on the vertical (Velocity.0.6 converted to metric). | meters per second | 0 to 1.77509 |
| `D.M.` | Depth at which the water sample was collected (Point.Depth converted to metric). | meters | 0 to 1.64592 |
| `DistanceFromShore.m.z` | Distance from the nearer bank, identical to DistanceFromShore.m; the analysis script standardizes this column for the models. | meters | 0 to 24 |
| `D.M..z` | Sample depth standardized to mean zero and unit standard deviation across the 192 samples; this column is the one used by the models. | standard deviations | -0.843549 to 2.61381 |
| `V.M..z` | Water velocity standardized to mean zero and unit standard deviation across the 192 samples. | standard deviations | -2.08069 to 2.1651 |
| `D.thal.m` | Distance of the sampling vertical from the thalweg, the vertical with the greatest discharge. | meters | 0 to 30 |
| `Type01` | Sample type as an indicator. |  | 0 = Surface grab (Type G); 1 = Point sample at depth (Type P) |

## 2021FIELD_DATA_Klamath.2.csv

Field sheet: one row per water sample collected (253 rows), including the surface and depth samples used in the analysis, the depth-integrated samples, extra samples, field blanks, and the hydraulic measurements at each vertical.

| Column | Definition | Units | Range or codes |
|---|---|---|---|
| `Sample ID` | Water sample identifier: site code, sample type, and vertical number. |  | Unique text identifier |
| `Sample description 1` | Sample label as written on the ddPCR plate map. |  | Text |
| `Chin.Conc.Eye` | O. tshawytscha ddPCR concentration from the first run, copies per microliter of reaction, with the positive/negative droplet threshold set by eye. | copies per microliter of reaction | 0 to 27.4 |
| `Cs.Conc.Eye` | C. shasta ddPCR concentration from the first run, copies per microliter of reaction, with the threshold set by eye. | copies per microliter of reaction | 0 to 1837 |
| `Chin.Conc.Eye.R` | O. tshawytscha ddPCR concentration from the replicate run of a subset of extracts (measurement error check). | copies per microliter of reaction | 0 to 9.92 |
| `Cs.Conc.Eye.R` | C. shasta ddPCR concentration from the replicate run of a subset of extracts (measurement error check). | copies per microliter of reaction | 0 to 1904 |
| `Used (1/2)` | Which ddPCR run was carried into the analysis for samples run twice. |  | 1 = First run; 2 = Replicate run |
| `Transect` | Vertical (transect station) number on the cross section, from the left bank. | vertical number | 0 to 20 |
| `Velocity 0.8` | Water velocity at 80 percent of total depth on the vertical. | feet per second | 0 to 4.4274 |
| `Velocity 0.6` | Water velocity at 60 percent of total depth on the vertical. | feet per second | 0 to 5.82377 |
| `Velocity 0.2` | Water velocity at 20 percent of total depth on the vertical. | feet per second | 0 to 6.92617 |
| `Point Depth` | Depth at which the sample was collected (0 for surface grabs). | feet | 0 to 6 |
| `Volume (ml)` | Volume of water filtered for the sample. | milliliters | 450 to 1840 |
| `DI total` | Depth-integrated samples: total depth of the vertical. | feet | 0 to 5.4 |
| `DI min` | Depth-integrated samples: shallowest depth sampled. | feet | 0 to 6 |
| `DI max` | Depth-integrated samples: deepest depth sampled. | feet | 0 to 7.2 |
| `Distance Left Bank(ft)` | Distance of the vertical from the left bank (looking downstream). | feet | 0 to 136.155 |
| `Distance From Shore` | Distance of the vertical from the nearer bank. | feet | 0 to 111.877 |
| `Site` | Cross section (site) code. |  | KLM = Site 1, just downstream of the Klamathon Bridge; I5 = Site 2, just downstream of the Interstate 5 bridge; TOH = Site 3, just downstream of the Tree of Heaven Campground; BVR = Site 4, just upstream of the Beaver Creek confluence; KMN = Site 5, just downstream of the Kinsman fish trapping site; SV = Site 6, in Seiad Valley just upstream of the Grider Creek confluence (river left) |
| `Type` | Sample type. |  | G = Surface grab, Whirl-Pak bag; P = Point sample at 60 percent of depth, US DH-59 sampler; DI = Depth-integrated sample, US DH-59 sampler; FB = Field blank, store-bought drinking water poured into a Whirl-Pak bag in the field; EX = Extra sample collected beyond the planned design |
| `Q cfs` | Total discharge of the cross section on the sampling day. | cubic feet per second | 1238.74 to 2358.01 |
| `Q L/s` | Total discharge of the cross section on the sampling day. | liters per second | 35077 to 66771.3 |
| `Q%` | Share of total discharge carried by the subsection around the vertical. | percent | 0 to 13.9789 |
| `SubQ (L/s)` | Discharge of the subsection around the vertical. | liters per second | 0 to 628560 |
| `JulDay` | Day of year of sampling in 2021. | day of year | 106 to 112 |
| `Elution volum (uL)` | Volume of buffer in which the extracted DNA was eluted. | microliters | 100 to 150 |

## 2021FIELD_DATA_Klamath_Discharge.csv

Discharge measurement by the midsection method: one row per subsection (vertical) of each cross section (107 rows).

| Column | Definition | Units | Range or codes |
|---|---|---|---|
| `subsection` | Subsection (vertical) number from the left bank. | subsection number | 0 to 20 |
| `Distance Left Bank meters` | Distance of the vertical from the left bank (looking downstream). | meters | 0 to 41.5 |
| `feet` | Distance of the vertical from the left bank. | feet | 0 to 136.155 |
| `Depth` | Total depth at the vertical, Ellsworth A-55 sounding reel. | feet | 0 to 8.55 |
| `clicks .8` | Price AA current meter revolutions counted at 80 percent of depth. | revolutions | 0 to 60 |
| `seconds .8` | Timing interval for the count at 80 percent of depth. | seconds | 0 to 45 |
| `fps.8` | Velocity at 80 percent of depth from the meter rating. | feet per second | 0 to 4.4274 |
| `clicks .6` | Price AA current meter revolutions counted at 60 percent of depth. | revolutions | 0 to 79 |
| `seconds .6` | Timing interval for the count at 60 percent of depth. | seconds | 0 to 50 |
| `fps.6` | Velocity at 60 percent of depth from the meter rating. | feet per second | 0 to 5.82377 |
| `clicks.2` | Price AA current meter revolutions counted at 20 percent of depth. | revolutions | 0 to 94 |
| `seconds.2` | Timing interval for the count at 20 percent of depth. | seconds | 0 to 40 |
| `fps.2` | Velocity at 20 percent of depth from the meter rating; blank where no count was made at that depth. | feet per second | 0 to 6.92617 |
| `Width Interval` | Width of the subsection represented by the vertical (midsection method). | feet | 0 to 26.2467 |
| `Area` | Cross-sectional area of the subsection. | square feet | 0 to 84.1535 |
| `3 point` | Mean velocity of the vertical from the three-point method. | feet per second | 0 to 5.34607 |
| `Q internal` | Discharge of the subsection. | cubic feet per second | 0 to 221.974 |
| `% total` | Share of the total cross-section discharge carried by the subsection. | percent | 0 to 13.9789 |
| `Site` | Cross section (site) code. |  | KLM = Site 1, just downstream of the Klamathon Bridge; I5 = Site 2, just downstream of the Interstate 5 bridge; TOH = Site 3, just downstream of the Tree of Heaven Campground; BVR = Site 4, just upstream of the Beaver Creek confluence; KMN = Site 5, just downstream of the Kinsman fish trapping site; SV = Site 6, in Seiad Valley just upstream of the Grider Creek confluence (river left) |

## MASTER.FINAL.ddPCRdata.Cshasta2021.csv

Droplet digital PCR output: one row per reaction (259 rows), including field samples, field blanks, replicate runs, and controls.

| Column | Definition | Units | Range or codes |
|---|---|---|---|
| `Sample.ID` | Sample identifier: site code, sample type (G, P, DI, FB, EX) and vertical number; no-template controls and replicate runs included. |  | Text identifier; a few identifiers repeat where a sample was run twice |
| `Chinook.c.ul` | O. tshawytscha target DNA concentration in the ddPCR reaction, Poisson-corrected, reported by Bio-Rad QX Manager. | copies per microliter of reaction | 0 to 27.4 |
| `Chinook.c.uL.lower` | Lower bound of the 95 percent Poisson confidence interval of Chinook.c.ul. | copies per microliter of reaction | 0 to 24.9 |
| `Chinook.c.uL.upper` | Upper bound of the 95 percent Poisson confidence interval of Chinook.c.ul. | copies per microliter of reaction | 0.168 to 30 |
| `Chinook.c.well` | O. tshawytscha target DNA copies per 20 microliter ddPCR reaction (Chinook.c.ul times 20). | copies per reaction | 0 to 549 |
| `Chinook.c.well.lower` | Lower bound of the 95 percent Poisson confidence interval of Chinook.c.well. | copies per reaction | 0 to 498 |
| `Chinook.c.well.upper` | Upper bound of the 95 percent Poisson confidence interval of Chinook.c.well. | copies per reaction | 3.36 to 600 |
| `C.shasta.c.uL` | C. shasta target DNA concentration in the ddPCR reaction, Poisson-corrected, reported by Bio-Rad QX Manager. | copies per microliter of reaction | 0 to 1837 |
| `C.shasta.c.uL.lower` | Lower bound of the 95 percent Poisson confidence interval of C.shasta.c.uL. | copies per microliter of reaction | 0 to 1798 |
| `C.shasta.c.uL.upper` | Upper bound of the 95 percent Poisson confidence interval of C.shasta.c.uL. | copies per microliter of reaction | 0 to 1878 |
| `C.shasta.c.well` | C. shasta target DNA copies per 20 microliter ddPCR reaction (C.shasta.c.uL times 20). | copies per reaction | 0 to 36746 |
| `C.shasta.c.well.lower` | Lower bound of the 95 percent Poisson confidence interval of C.shasta.c.well. | copies per reaction | 0 to 35960 |
| `C.shasta.c.well.upper` | Upper bound of the 95 percent Poisson confidence interval of C.shasta.c.well. | copies per reaction | 0 to 37560 |
| `Accepted Droplets` | Number of droplets read by the Bio-Rad QX200 droplet reader for the reaction; reactions with fewer than 10,000 droplets were rerun. | droplets | 8993 to 21032 |
