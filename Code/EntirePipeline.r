# 1) Clone the repository from https://github.com/degama01/TCDDMeth
# 2) Install the packages listed in ".README.md"
# 3) Set working directory to "/path/to/TCDDMeth/"
# 4) Run the following scripts:
# Note that the scripts must be run sequentially.
# RNA-seq Analysis
source("./Code/RNA-seq.R") # Prepare required data

# DNA Methylation Analysis
## Note that the output files from the RNA-seq Analysis are needed to run the DNA Methylation Analysis code.
source("./Code/Methylation_Analysis.R")
