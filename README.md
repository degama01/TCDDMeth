# Dioxin Disrupts Dynamic DNA Methylation Patterns in Genes that Govern Cardiomyocyte Maturation
***


de Gannes, M., Ko, CI., Koch, S., Zhang, X., Biesiada, J.,
Medvedovic, M., Rubinstein, J.,  & Puga A.

## Data
***

All data files required to run the code used to perform the analysis is included in the ```./Data``` folder and will be downloaded there when this repository is cloned.

## Reference Files
***

All raw and processed data from our study used to derive the files for the analysis can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/) and saved to ```./Reference Files/```. They are not required to run the code included in this repository.

## Analyses
***

Scripts for all analyses are saved to ```./Code/```. All scripts should be run with this directory set as the working directory (```./```).

All analyses were performed in [R](https://www.r-project.org/) version 4.0.0 on Windows 10, using the following packages:

* tidyverse
* circlize
* RColorBrewer
* GenomicFeatures
* ComplexHeatmap
* edgeR
* org.Mm.eg.db
* methylKit

These packages can be installed in R by running the following:

```
install.packages("tidyverse")
install.packages("circlize")
install.packages("RColorBrewer")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicFeatures")
BiocManager::install("ComplexHeatmap")
BiocManager::install("edgeR")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("methylKit")

```

### To run the analyses:

  1. Clone this repository.
  2. Install the packages listed above.
  3. Run the following in R:
  ```
  setwd("/path/to/TCDDMeth")
  source("./Code/EntirePipeline.r")
  ```
  
## Output
***
Output from all folders is saved to the ```./Output/``` directory, with figures saved to ```./Output/Figures/```. Note that the output generated from ```RNA-seq.R``` is required to run ```Methylation_Analysis.R```.