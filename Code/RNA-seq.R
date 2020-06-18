#' # Code Input
#' 
#' ## Preparation of Data for Usage
## ------------------------------------------------------------------------------------------------------------
# Loading Relevant Libraries
library(tidyverse)
library(GenomicFeatures)
library(ComplexHeatmap)
library(circlize)
library(edgeR)

#' 
## ------------------------------------------------------------------------------------------------------------
# Loading the raw and processed data from the Data folder
file = "./Data/table_S2_revised.txt"
raw_data = load("./Data/raw.data.rda")
res_list = load("./Data/res.list.rda")

#' 
#' ### Need to match gene symbols with ID values in raw data Counts and rpkm Counts Table
#' 
## ------------------------------------------------------------------------------------------------------------
# First, make the counts and rpkm counts into two dataframes
CountsTable.df <- data.frame(id = row.names(raw.data$countsTable), raw.data$countsTable)
ScaledCounts.df <- data.frame(id = row.names(raw.data$rpkmCountsTable), raw.data$rpkmCountsTable) # 24571 genes

# Then, make them into tbl_df format for plyr and Tidyverse
CountsTable.df <- tbl_df(CountsTable.df)
ScaledCounts.df <- tbl_df(ScaledCounts.df)

# Now take a list of symbols from res.list along with matching IDs by combining the dataframe comparisons, keeping only id and symbol columns
Symbols_1.df <- data.frame(res.list$B1DMSO24_vs_ESC$gene.result)
Symbols_1.df <- Symbols_1.df %>%
  dplyr::select(id, symbol) %>%
  print()
Symbols_2.df <- data.frame(res.list$B1TCDD24_vs_B1DMSO24$gene.result)
Symbols_2.df <- Symbols_2.df %>%
  dplyr::select(id, symbol) %>%
  print()
Symbols_3.df <- data.frame(res.list$B1TCDD72_vs_B1DMSO24$gene.result)
Symbols_3.df <- Symbols_3.df %>%
  dplyr::select(id, symbol) %>%
  print()
Symbols_4.df <- data.frame(res.list$B1TCDD96_vs_B1DMSO24$gene.result)
Symbols_4.df <- Symbols_4.df %>%
  dplyr::select(id, symbol) %>%
  print()

# Combine the dataframes
Gene_Symbols.df <- do.call("rbind", list(Symbols_1.df, Symbols_2.df, Symbols_3.df, Symbols_4.df))
Gene_Symbols.df <- tbl_df(Gene_Symbols.df)

#' 
## ------------------------------------------------------------------------------------------------------------
# Create table of counts with rownames column name being "Symbol" as in edgeR page 10
CountsTable.df$symbol <- Gene_Symbols.df$symbol[match(CountsTable.df$id, Gene_Symbols.df$id)]
CountsTable.df <- CountsTable.df %>%
  drop_na() %>%
  dplyr::select(id, symbol, APM9, APM10, APM1, APM2, APM3, APM4, APM5, APM6, APM7, APM8) %>%
  print()
CountsTable.df <- tbl_df(CountsTable.df)
Genes <- CountsTable.df$symbol
x <- CountsTable.df %>%
  dplyr::select(-symbol, -id) %>%
  print()
# Make x into a matrix and add symbols as rownames
x <- as.matrix(x)
row.names(x) <- Genes

#' 
## ------------------------------------------------------------------------------------------------------------
# Now create the groups and relevant factors
group <- factor(c(1,1,2,2,3,3,4,4,5,5))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

#' 
## ------------------------------------------------------------------------------------------------------------
# Perform quasi-likelihood F-tests
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef= 2:5) # performs comparisons for each group  in comparison to group 1 (DMSO24h)
topTags(qlf)

# Determine number of genes with FDR value < 0.05
#temp = p.adjust(qlf$table$PValue, method = "BH")
#sum(temp < 0.05) #4821 genes


#' 
#' # Prepare a Table to Put Into Excel of All Genes and Fold Changes to Use for IPA Analysis
#' 
## ------------------------------------------------------------------------------------------------------------
# Create a new dataframe called FC.df from the qlf table that was made.
FC.df <- as.data.frame(qlf$table)
FC.df<- cbind(symbol = rownames(FC.df), FC.df)
FC.df$symbol <- as.character(FC.df$symbol)
FC.df <- tbl_df(FC.df)

# Add column for all FDR values
FC.df$FDR <- p.adjust(FC.df$PValue, method = "BH")

# Add in column for gene IDs
FC.df$id <- Gene_Symbols.df$id[match(FC.df$symbol, Gene_Symbols.df$symbol)]

FC.df <- FC.df %>%
  dplyr::select(id, symbol, logFC.group2, logFC.group3, logFC.group4, logFC.group5, F, PValue, FDR) %>%
  print()

# Now select for all genes with FDR < 0.05, then order by lowest FDR values
FC_filtered.df <- FC.df %>%
  filter(FDR < 0.05) %>% #4821 genes
  arrange(FDR) %>%
  print()

# Write to csv file
write.csv(FC_filtered.df, "./Output/TCDD_B1_FDR_Filtered.csv")

# Write csv file with columns labelled properly
Supp_table_1.df <- FC_filtered.df %>%
  rename("logFC.group2"="logFC_ES", "logFC.group3"="logFC_TCDD24", "logFC.group4"="logFC_TCDD72", "logFC.group5"="logFC_TCDD96") %>%
  print()

write.csv(Supp_table_1.df, "./Output/Supp_table.csv")


#' 
#' 
## ------------------------------------------------------------------------------------------------------------
# Now, take the original FC.df with ALL the (filtered and unfiltered) then add FDR columns for both Raw Counts and Scaled Counts. 
CountsTable.df$FDR <- FC.df$FDR[match(CountsTable.df$symbol, FC.df$symbol)]
ScaledCounts.df$symbol <- Gene_Symbols.df$symbol[match(ScaledCounts.df$id, Gene_Symbols.df$id)]
ScaledCounts.df <- ScaledCounts.df %>%
  drop_na() %>%
  dplyr::select(id, symbol, APM1, APM2, APM9, APM10, APM3, APM4, APM5, APM6, APM7, APM8) %>%
  print()
ScaledCounts.df$FDR <- FC.df$FDR[match(ScaledCounts.df$symbol, FC.df$symbol)]

# Confirm  FDR values match with symbols
CountsTable.df
ScaledCounts.df

#' 
## ------------------------------------------------------------------------------------------------------------
CountsTable.df <- CountsTable.df %>%
  drop_na() %>%
  dplyr::select(symbol, FDR, APM1, APM2, APM9, APM10, APM3, APM4, APM5, APM6, APM7, APM8) %>%
  filter(FDR < 0.05) %>%
  print()
ScaledCounts.df <- ScaledCounts.df %>%
  drop_na() %>%
  dplyr::select(symbol, FDR, APM1, APM2, APM9, APM10, APM3, APM4, APM5, APM6, APM7, APM8) %>%
  filter (FDR< 0.05) %>%
  print()
ScaledCounts.df <- tbl_df(ScaledCounts.df)
# Create vector of genes
Genes <- as.vector(CountsTable.df$symbol)

# Remove symbol and FDR column from each dataframe then add in rownames as genes.
CountsTable.df <- CountsTable.df %>%
  dplyr::select(-symbol, -FDR) %>%
  print()
row.names(CountsTable.df) <- Genes

ScaledCounts.df <- ScaledCounts.df %>%
  dplyr::select(-symbol, -FDR) %>%
  print()
row.names(ScaledCounts.df) <- Genes

#' 
#' 
## ------------------------------------------------------------------------------------------------------------
# Set up for gene expression matrix
mat = as.matrix(CountsTable.df) #4821 genes
base_mean = rowMeans(mat)
log_rpkm <- log10(ScaledCounts.df + 1) # Log transformed original rpkms
mat_scaled.df = log_rpkm - rowMeans(log_rpkm) 
# Write the scaled matrix to csv file for use in making the composite heatmap separately in this same folder
write.csv(mat_scaled.df, "./Output/exp_scaled.csv")

#' 
## ------------------------------------------------------------------------------------------------------------

