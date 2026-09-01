# iSBio - Soil Virome Biogeography
Soil metagenomics (PRJNA1045969) from the International Soil Biogeography Consortium (iSBio)

Analysis pipeline for "Biogeography of soil viral communities across major natural terrestrial ecosystems" across forest, grassland, shrubland, and tundra ecosystems (314 soil metagenomes, 33 countries, six continents), using consistent sampling, sequencing, and analytical procedures throughout.

Across 314 soil metagenomes, we recovered ~19,000 viral operational taxonomic units (vOTUs), substantially expanding known soil viral diversity. This repository contains the R-based analysis workflow used to:

Characterize vOTU diversity — alpha/beta diversity (Shannon index, PCoA/NMDS)
Identify environmental drivers of viral biogeography:dbRDA and random forest variable importance across environmental and soil predictors
Disentangle environmental filtering from spatial structuring: Mantel and partial Mantel tests separating geographic distance, elevational distance, and environmental distance effects on viral community
Link viruses to predicted hosts: using metagenome-assembled genomes (MAGs), including virus-to-host abundance ratios (VHR) comparisons across ecosystems

Findings show pronounced ecosystem-level differentiation in viral richness, composition, predicted lifestyle, and virus-host relationships. While soil pH is the strongest overall environmental predictor, its relative importance varies by ecosystem — forest viral communities are predominantly pH-associated, while elevation and soil texture play stronger roles in grasslands.

# System requirements

The analysis workflow was tested using:

- R version 4.5.0
- Operating system: Windows 11
- No non-standard hardware is required. Some analyses, particularly random forest analyses and community-distance calculations, may benefit from a computer with sufficient memory.

# R package dependencies
- phyloseq
- vegan
- ggplot2
- dplyr
- tidyr
- readr
- randomForest  

Repository structure:
- 01_alpha_beta_diversity.R — Shannon/richness, PCoA, NMDS, PERMANOVA
- 02_dbRDA_drivers.R — distance-based redundancy analysis
- 03_random_forest.R — variable importance for community structure
- 04_mantel_env_correlogram.R -testing vOTU/vFunctions composition against environmental variables
- 05_mantel_distance_decay.R — geographic/elevational/environmental Mantel tests

The processed vOTU abundance table (`count.txt`) and sample metadata (`isbiometa.txt`) are provided in this repository. These materials reproduce the principal community-diversity, environmental-driver, and spatial analyses reported in the manuscript. Analyses of viral lifestyle, virus–host relationships, and auxiliary metabolic genes are not included in this repository.

License

The analysis code is released under the MIT License. See the LICENSE file
for details.
