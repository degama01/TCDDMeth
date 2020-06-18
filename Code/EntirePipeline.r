# 1) Clone the repository from https://github.com/degama01/TCDDMeth
# 2) Download the required data as described in ".README.md" and in "./Data/README.md"
# 3) Install the packages listed in ".README.md"
# 4) Set working directory to "/path/to/TCDDMeth/"
# 5) Run the following scripts:
# Note that the scripts must be run sequentially.
# RNA-seq Analysis
source("./Code/RNA-seq.R") # Prepare required data

# DNA Methylation Analysis
## Note that the output files from the RNA-seq Analysis are needed to run the DNA Methylation Analysis code.
source("./Code/Methylation_Analysis.R")
