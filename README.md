# Dioxin Disrupts Dynamic DNA Methylation Patterns in Genes that Govern Cardiomyocyte Maturation
***


de Gannes, M., Ko, CI., Zhang, X., Biesiada, J.,
Medvedovic, M.,  & Puga A.

## Data
***

Here, we use three datasets generated from our study: 1) raw.data.rda, 2) res.list.rda, and 3) myobj.rds. See ```./Data/README.md``` for details.

* Raw and processed data were downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/) and saved to ```./Data/```.
* Methylation counts data that were processed in our study are included in ```./Data/``` and will be downloaded by cloning this repository.

## Reference Files
***

All reference files are included in this repository, and will be downloaded by cloning it.

## Analyses
***

Scripts for all analyses are saved to ```./Code/```. All scripts should be run with this directory set as the working directory (```./```).

All analyses were performed in [RStudio](https://rstudio.com/) version 1.3 on Windows 10, using the following packages:

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

  1. Clone this repository
  2. Download all the required data as described here and in ```./Data/README.md```
  3. Install the packages listed above
  4. Run the following in RStudio:
  ```
  setwd("/path/to/TCDDMeth")
  source("./Code/EntirePipeline.Rmd")
  ```
  
## Output
***
Output from all folders is saved to the ```./Output/``` directory, with figures saved to ```./Output/Figures/```. Note that the output generated from ```RNA-seq.R``` is required to run ```Composite_Heatmap.R```.