#' # Set up relevant libraries
## -----------------------------------------------------------------------------------------------------------
# Loading Libraries
library(tidyverse)
library(ComplexHeatmap)
library(GenomicFeatures)
library(circlize)
library(edgeR)
library(org.Mm.eg.db)
library(methylKit)
library(RColorBrewer)
library(ggplot2)

#' 
#' # Importing Datasets from Jacek and Scaled Gene Expression Table
#' 
## -----------------------------------------------------------------------------------------------------------
# Loading the data files
#load(file="./Data/myobj.rds")
#raw.df=unite(myobj, destrand=FALSE)
expr.df <- read.csv("./Output/exp_scaled.csv")

samples.names =
list("APM11B","APM12B","APM13B","APM14B","APM15B","APM16B","APM17","APM18","APM19B","APM20B")
# tr =
list("ES_None_0","ES_None_0","B1_TCDD_24","B1_TCDD_24","B1_TCDD_72","B1_TCDD_72","B1_TCDD_96","B1_TCDD_9
6","B1_DMSO_24","B1_DMSO_24")
tr = c(0,0,1,1,2,2,3,3,4,4)


#' # Start setting up the large dataframe as the original comparisons were done, this time including the B1TCDD24vsB1DMSO24 comparison. Use the EdgeR methylation as an example to guide the set up, annotations, and labelling of columns.
#' 
## -----------------------------------------------------------------------------------------------------------
#working.df <- raw.df
#working.df <- tbl_df(working.df)
#working.df$chr = str_remove_all(working.df$chr, "[chr]")


#' 
#' 
## -----------------------------------------------------------------------------------------------------------
#working.df <- working.df %>%
  #rename(start = "Locus") %>%
  #dplyr::select(-end, -strand, -starts_with("numTs")) %>%
  #rename(coverage1="APM11_coverage", coverage2="APM12_coverage", coverage3="APM13_coverage", coverage4="APM14_coverage", coverage5="APM15_coverage", coverage6="APM16_coverage", coverage7="APM17_coverage", coverage8="APM18_coverage", coverage9="APM19_coverage", coverage10="APM20_coverage") %>%
  #rename(numCs1="APM11_numCs", numCs2="APM12_numCs", numCs3="APM13_numCs", numCs4="APM14_numCs", numCs5="APM15_numCs", numCs6="APM16_numCs", numCs7="APM17_numCs", numCs8="APM18_numCs", numCs9="APM19_numCs", numCs10="APM20_numCs") %>%
  #mutate(ratioAPM11=APM11_numCs/APM11_coverage, ratioAPM12=APM12_numCs/APM12_coverage, ratioAPM13=APM13_numCs/APM13_coverage, ratioAPM14=APM14_numCs/APM14_coverage, ratioAPM15=APM15_numCs/APM15_coverage, ratioAPM16=APM16_numCs/APM16_coverage, ratioAPM17=APM17_numCs/APM17_coverage, ratioAPM18=APM18_numCs/APM18_coverage, ratioAPM19=APM19_numCs/APM19_coverage, ratioAPM20=APM20_numCs/APM20_coverage) %>%
  #rename(chr = "Chr") %>%
  #drop_na() %>%
  #print()


#' 
#' # For convenience, we sort the DGEList so that all loci are in genomic order, from chromosome 1 to chromosome Y.
#' 
## -----------------------------------------------------------------------------------------------------------
#ChrNames <- c(1:19,"X","Y")
#working.df$Chr <- factor(working.df$Chr, levels=ChrNames)
#o <- order(working.df$Chr, working.df$Locus)
#working.df <- working.df[o,]

#' 
#' 
#' # Obtain dataframes for each of the comparisons for ALL the genes.
#' 
## -----------------------------------------------------------------------------------------------------------
ES_all.df <- read.csv("./Data/ES_all.csv") %>%
  #dplyr::select(Chr, Locus, APM11_numCs, APM12_numCs, APM19_numCs, APM20_numCs, APM11_coverage, APM12_coverage, APM19_coverage, APM20_coverage, ratioAPM11, ratioAPM12, ratioAPM19, ratioAPM20) %>%
  print()

DMSO_24_all.df <- read.csv("./Data/DMSO_24_all.csv") %>%
  #dplyr::select(Chr, Locus, APM19_numCs, APM20_numCs, APM19_coverage, APM20_coverage, ratioAPM19, ratioAPM20) %>%
  print()

TCDD_24_all.df <- read.csv("./Data/TCDD_24_all.csv") %>%
  #dplyr::select(Chr, Locus, APM13_numCs, APM14_numCs, APM19_numCs, APM20_numCs, APM13_coverage, APM14_coverage, APM19_coverage, APM20_coverage, ratioAPM13, ratioAPM14, ratioAPM19, ratioAPM20) %>%
  print()

TCDD_72_all.df <- read.csv("./Data/TCDD_72_all.csv") %>%
  #dplyr::select(Chr, Locus, APM15_numCs, APM16_numCs, APM19_numCs, APM20_numCs, APM15_coverage, APM16_coverage, APM19_coverage, APM20_coverage, ratioAPM15, ratioAPM16, ratioAPM19, ratioAPM20) %>%
  print()

TCDD_96_all.df <- read.csv("./Data/TCDD_96_all.csv") %>%
  #dplyr::select(Chr, Locus, APM17_numCs, APM18_numCs, APM19_numCs, APM20_numCs, APM17_coverage, APM18_coverage, APM19_coverage, APM20_coverage, ratioAPM17, ratioAPM18, ratioAPM19, ratioAPM20) %>%
  print()
  

#' 
#' # Now add in the columns for the means of your tested and reference group (B1DMSO24), beta, and direction of change
## -----------------------------------------------------------------------------------------------------------
ES_all.df <- ES_all.df %>%
  mutate(meanES = rowMeans(dplyr::select(., ratioAPM11:ratioAPM12))) %>%
  mutate(meanB1DMSO24h = rowMeans(dplyr::select(., ratioAPM19:ratioAPM20))) %>%
  mutate(beta = meanB1DMSO24h - meanES) %>%
  mutate(directionES = ifelse(beta<0, "hyper", ifelse(beta>0, "hypo", "no change"))) %>%
  print()

DMSO_24_all.df <- DMSO_24_all.df %>%
  mutate(meanB1DMSO24h = rowMeans(dplyr::select(., ratioAPM19:ratioAPM20))) %>%
  mutate(beta = meanB1DMSO24h - meanB1DMSO24h) %>%
  mutate(directionB1DMSO24h = ifelse(beta<0, "hyper", ifelse(beta>0, "hypo", "no change"))) %>%
  print()

TCDD_24_all.df <- TCDD_24_all.df %>%
  mutate(meanB1TCDD24h = rowMeans(dplyr::select(., ratioAPM13:ratioAPM14))) %>%
  mutate(meanB1DMSO24h = rowMeans(dplyr::select(., ratioAPM19:ratioAPM20))) %>%
  mutate(beta = meanB1DMSO24h - meanB1TCDD24h) %>%
  mutate(directionB1TCDD24h = ifelse(beta<0, "hyper", ifelse(beta>0, "hypo", "no change"))) %>%
  print()

TCDD_72_all.df <- TCDD_72_all.df %>%
  mutate(meanB1TCDD72h = rowMeans(dplyr::select(., ratioAPM15:ratioAPM16))) %>%
  mutate(meanB1DMSO24h = rowMeans(dplyr::select(., ratioAPM19:ratioAPM20))) %>%
  mutate(beta = meanB1DMSO24h - meanB1TCDD72h) %>%
  mutate(directionB1TCDD72h = ifelse(beta<0, "hyper", ifelse(beta>0, "hypo", "no change"))) %>%
  print()

TCDD_96_all.df <- TCDD_96_all.df %>%
  mutate(meanB1TCDD96h = rowMeans(dplyr::select(., ratioAPM17:ratioAPM18))) %>%
  mutate(meanB1DMSO24h = rowMeans(dplyr::select(., ratioAPM19:ratioAPM20))) %>%
  mutate(beta = meanB1DMSO24h - meanB1TCDD96h) %>%
  mutate(directionB1TCDD96h = ifelse(beta<0, "hyper", ifelse(beta>0, "hypo", "no change"))) %>%
  print()

#' # Now load all the dataframes with the significant genes for each comparison
#' 
## -----------------------------------------------------------------------------------------------------------
# Load the dataframes
# Loading the data from Jacek
raw_ES.df <- read.csv("./Data/ESvsB1DMSO24.csv")
raw_72.df <- read.csv("./Data/B1TCDD72vsB1DMSO24.csv")
raw_96.df <- read.csv("./Data/B1TCDD96vsB1DMSO24.csv")


#' 
#' # For each of the comparison files, count number of total significant changes, and number of hypo versus hyper methylated sites
#' 
## -----------------------------------------------------------------------------------------------------------

# Remove NAs from beta column
raw_ES.df <- raw_ES.df %>%
  drop_na(beta) %>%
  print()

raw_72.df <- raw_72.df %>%
  drop_na(beta) %>%
  print()

raw_96.df <- raw_96.df %>%
  drop_na(beta) %>%
  print()


#' # Create supplementary tables
## -----------------------------------------------------------------------------------------------------------
Supp_ES.df <- raw_ES.df %>%
  dplyr:: select(names, LOCATION, GENEID, symbol, beta, p.value) %>%
  rename("names" = "Names", "LOCATION"="Location", "GENEID"="Gene ID", "symbol"="Symbol", "p.value"="p value") %>%
  print()
write.csv(Supp_ES.df, "./Output/Supp_ES.csv")

Supp_72.df <- raw_72.df %>%
  dplyr:: select(names, LOCATION, GENEID, symbol, beta, p.value) %>%
  rename("names" = "Names", "LOCATION"="Location", "GENEID"="Gene ID", "symbol"="Symbol", "p.value"="p value") %>%
  print()
write.csv(Supp_72.df, "./Output/Supp_72.csv")

Supp_96.df <- raw_96.df %>%
  dplyr:: select(names, LOCATION, GENEID, symbol, beta, p.value) %>%
  rename("names" = "Names", "LOCATION"="Location", "GENEID"="Gene ID", "symbol"="Symbol", "p.value"="p value") %>%
  print()
write.csv(Supp_96.df, "./Output/Supp_96.csv")


#' 
#' 
## -----------------------------------------------------------------------------------------------------------

# ESvsB1DMSO24
hypo_ES <- sum(raw_ES.df$beta > 0 & -log10(raw_ES.df$p.value) >= 5) #2616 sites

hyper_ES <- sum(raw_ES.df$beta < 0 & -log10(raw_ES.df$p.value) >= 5) #3 sites

total_ES <- sum(hypo_ES+hyper_ES) #2619 sites

#24vsB1DMSO24
## no significant changes

# 72vsB1DMSO24
hypo_72 <- sum(raw_72.df$beta > 0) # 59 sites

hyper_72 <- sum(raw_72.df$beta < 0) #514 sites

total_72 <- sum(hypo_72+hyper_72) #573 sites

# 96vsB1DMSO24

hypo_96 <- sum(raw_96.df$beta > 0) # 0 sites

hyper_96 <- sum(raw_96.df$beta < 0) #3026 sites

total_96 <- sum(hypo_96+hyper_96) #3026 sites


#' 
#' 
#' 
#' 
## -----------------------------------------------------------------------------------------------------------
# Rename the labels for the samples in the expression matrix so that they are the same as those for the methylation datasets
expr.df <- expr.df %>%
  dplyr::rename(APM11="APM1", APM12="APM2", APM19="APM9", APM20="APM10", APM13="APM3", APM14="APM4", APM15="APM5", APM16="APM6", APM17="APM7", APM18="APM8") %>%
  print()


#' 
#' # Create a vector of labels for the samples called "type" as in the Composite Heatmap example
## -----------------------------------------------------------------------------------------------------------
type <- c("ES", "ES", "B1DMSO24", "B1DMSO24", "B1TCDD24", "B1TCDD24", "B1TCDD72", "B1TCDD72", "B1TCDD96", "B1TCDD96")

#' 
#' 
#' # Remove irrelevant columns from the methylation datasets, filter out genomic targets with no symbol, and also add columns for Chromosome Number, Locus, and the direction of the methylation change
#' 
## -----------------------------------------------------------------------------------------------------------
ES.df <- raw_ES.df %>%
  dplyr::select(names, seqnames, start, end, LOCATION, GENEID, symbol, meanES, beta, p.value) %>%
  drop_na() %>%
  mutate(directionES = ifelse(beta<0, "hyper", "hypo")) %>%
  rename(seqnames = "Chr") %>%
  rename(start = "Locus") %>%
  dplyr::select (-end) %>%
  print()
ES.df$Chr = gsub("[^[:digit:]]", "", ES.df$Chr)
ES.df <- tbl_df(ES.df)
ES.df<- ES.df %>%
  mutate(Chr = case_when(str_sub(names, 4, 4) == 'X' ~ 'X',str_sub(names, 4, 4) == 'Y' ~ 'Y',TRUE ~ Chr)) %>%
  print()

DMSO_24.df <- raw_72.df %>%
  dplyr::select(names, seqnames, start, end, LOCATION, GENEID, symbol, meanB1DMSO24h, beta, p.value) %>%
  drop_na() %>%
  mutate(directionB1DMSO24h = "no change") %>%
  rename(seqnames = "Chr") %>%
  rename(start = "Locus") %>%
  dplyr::select (-end) %>%
  print()
DMSO_24.df$Chr = gsub("[^[:digit:]]", "", DMSO_24.df$Chr)
DMSO_24.df <- tbl_df(DMSO_24.df)
DMSO_24.df<- DMSO_24.df %>%
  mutate(Chr = case_when(str_sub(names, 4, 4) == 'X' ~ 'X',str_sub(names, 4, 4) == 'Y' ~ 'Y',TRUE ~ Chr)) %>%
  print()

#TCDD_24.df <- raw_24.df %>%
 # dplyr::select(names, seqnames, start, end, LOCATION, GENEID, symbol, meanB1TCDD24h, beta) %>%
  #drop_na() %>%
  #mutate(directionB1TCDD24h = ifelse(beta<0, "hyper", "hypo")) %>%
  #rename(seqnames = "Chr") %>%
  #rename(start = "Locus") %>%
  #dplyr::select (-end) %>%
  #print()
#TCDD_24.df$Chr = gsub("[^[:digit:]]", "", TCDD_24.df$Chr)
# TCDD_24.df <- tbl_df(TCDD_24.df)
#TCDD_24.df<- TCDD_24.df %>%
  #mutate(Chr = case_when(str_sub(names, 4, 4) == 'X' ~ 'X',str_sub(names, 4, 4) == 'Y' ~ 'Y',TRUE ~ Chr)) %>%
  #print()

TCDD_72.df <- raw_72.df %>%
  dplyr::select(names, seqnames, start, end, LOCATION, GENEID, symbol, meanB1TCDD72h, beta, p.value) %>%
  drop_na() %>%
  mutate(directionB1TCDD72h = ifelse(beta<0, "hyper", "hypo")) %>%
  rename(seqnames = "Chr") %>%
  rename(start = "Locus") %>%
  dplyr::select (-end) %>%
  print()
TCDD_72.df$Chr = gsub("[^[:digit:]]", "", TCDD_72.df$Chr)
TCDD_72.df <- tbl_df(TCDD_72.df)
TCDD_72.df<- TCDD_72.df %>%
  mutate(Chr = case_when(str_sub(names, 4, 4) == 'X' ~ 'X',str_sub(names, 4, 4) == 'Y' ~ 'Y',TRUE ~ Chr)) %>%
  print()

TCDD_96.df <- raw_96.df %>%
  dplyr::select(names, seqnames, start, end, LOCATION, GENEID, symbol, meanB1TCDD96h, beta, p.value) %>%
  drop_na() %>%
  mutate(directionB1TCDD96h = ifelse(beta<0, "hyper", "hypo")) %>%
  rename(seqnames = "Chr") %>%
  rename(start = "Locus") %>%
  dplyr::select (-end) %>%
  print()

TCDD_96.df$Chr = gsub("[^[:digit:]]", "", TCDD_96.df$Chr)
TCDD_96.df <- tbl_df(TCDD_96.df)
TCDD_96.df<- TCDD_96.df %>%
  mutate(Chr = case_when(str_sub(names, 4, 4) == 'X' ~ 'X',str_sub(names, 4, 4) == 'Y' ~ 'Y',TRUE ~ Chr)) %>%
  print()

#' 
#' # Add in a column for the genomic distance from the CpG to the TSS of the associated gene
#' 
## -----------------------------------------------------------------------------------------------------------
TSS <- nearestTSS(ES.df$Chr, ES.df$Locus, species="Mm")
ES.df$Strand <- TSS$strand
ES.df$Distance <- TSS$distance
ES.df$Width <- TSS$width

TSS <- nearestTSS(DMSO_24.df$Chr, DMSO_24.df$Locus, species="Mm")
DMSO_24.df$Strand <- TSS$strand
DMSO_24.df$Distance <- TSS$distance
DMSO_24.df$Width <- TSS$width

TSS <- nearestTSS(TCDD_72.df$Chr, TCDD_72.df$Locus, species="Mm")
TCDD_72.df$Strand <- TSS$strand
TCDD_72.df$Distance <- TSS$distance
TCDD_72.df$Width <- TSS$width

TSS <- nearestTSS(TCDD_96.df$Chr, TCDD_96.df$Locus, species="Mm")
TCDD_96.df$Strand <- TSS$strand
TCDD_96.df$Distance <- TSS$distance
TCDD_96.df$Width <- TSS$width


#' 
#' # Create new dataframes for all significant targets for each comparison and combine them.
## -----------------------------------------------------------------------------------------------------------

# Set up initial dataframes
ES_select.df <- ES.df %>%
  dplyr::select(Chr, Locus, LOCATION, GENEID, symbol, Strand, Distance, Width, p.value) %>%
  group_by(symbol) %>% arrange(p.value) %>% 
  filter(row_number() == 1) %>%
  print()
ES_select.df <- tbl_df(ES_select.df)

DMSO_24_select.df <- DMSO_24.df %>%
  dplyr::select(Chr, Locus, LOCATION, GENEID, symbol, Strand, Distance, Width, p.value) %>%
  group_by(symbol) %>% arrange(p.value) %>% 
  filter(row_number() == 1) %>%
  print()
DMSO_24_select.df <- tbl_df(DMSO_24_select.df)

TCDD_72_select.df <- TCDD_72.df %>%
  dplyr::select(Chr, Locus, LOCATION, GENEID, symbol, Strand, Distance, Width, p.value) %>%
  group_by(symbol) %>% arrange(p.value) %>% 
  filter(row_number() == 1) %>%
  print()
TCDD_72_select.df <- tbl_df(TCDD_72_select.df)

TCDD_96_select.df <- TCDD_96.df %>%
  dplyr::select(Chr, Locus, LOCATION, GENEID, symbol, Strand, Distance, Width, p.value) %>%
  group_by(symbol) %>% arrange(p.value) %>% 
  filter(row_number() == 1) %>%
  print()
TCDD_96_select.df <- tbl_df(TCDD_96_select.df)



#' 
#' 
#' 
## -----------------------------------------------------------------------------------------------------------
# Combine the dataframes
Merged_select.df <- rbind(ES_select.df, DMSO_24_select.df, TCDD_72_select.df, TCDD_96_select.df)
Merged_select.df <- Merged_select.df %>%
  ungroup(symbol) %>%
  arrange(symbol) %>%
  print()

# Remove duplicates by Locus
Merged_select.df <- distinct(Merged_select.df, Locus, .keep_all = TRUE)

# Keep only loci with smallest p-value remaining
Merged_select.df <- Merged_select.df %>%
  group_by(symbol) %>% arrange(p.value) %>% 
  filter(row_number() == 1) %>%
  ungroup(symbol) %>%
  arrange(symbol) %>%
  print()

#' 
#' 
#' # Now determine the list of significant DMRs overlapping with your expression genes and filter for only those genes in a new dataset for methylation, and a new dataset for gene expression
#' 
## -----------------------------------------------------------------------------------------------------------
Gene_List <- intersect(Merged_select.df$symbol, expr.df$X)

meth_filtered.df <- Merged_select.df %>%
  filter(symbol %in% Gene_List)%>%
  print()

expr_filtered.df <- expr.df %>%
  rename(X = "Symbol") %>%
  filter(Symbol %in% Gene_List) %>%
  arrange(Symbol) %>%
  print()


#' 
#' # Now merge each of the dataframe comparisons containing ALL genes with the filtered methylated by locus, removing any symbols nonexistant (aka removing NAs) as these would be represent all non-significant genes
#' 
## -----------------------------------------------------------------------------------------------------------

ES_merged.df <- merge(ES_all.df,meth_filtered.df,by="Locus")
ES_merged.df <- ES_merged.df %>%
  rename(Chr.x = "Chr") %>%
  rename(LOCATION = "Location") %>%
  rename(symbol = "Symbol") %>%
  dplyr::select(Chr, Locus, Location, GENEID, Symbol, Strand, Distance, Width, ratioAPM11, ratioAPM12, meanES, directionES, p.value) %>%
  filter(Symbol %in% Gene_List)%>%
  drop_na() %>%
  print()
ES_merged.df <- tbl_df(ES_merged.df)

DMSO_24_merged.df <- merge(DMSO_24_all.df,meth_filtered.df,by="Locus")
DMSO_24_merged.df <- DMSO_24_merged.df %>%
  rename(Chr.x = "Chr") %>%
  rename(LOCATION = "Location") %>%
  rename(symbol = "Symbol") %>%
  dplyr::select(Chr, Locus, Location, GENEID, Symbol, Strand, Distance, Width, ratioAPM19, ratioAPM20, meanB1DMSO24h, directionB1DMSO24h, p.value) %>%
  filter(Symbol %in% Gene_List)%>%
  drop_na() %>%
  print()
DMSO_24_merged.df <- tbl_df(DMSO_24_merged.df)

TCDD_24_merged.df <- merge(TCDD_24_all.df,meth_filtered.df,by="Locus")
TCDD_24_merged.df <- TCDD_24_merged.df %>%
  rename(Chr.x = "Chr") %>%
  rename(LOCATION = "Location") %>%
  rename(symbol = "Symbol") %>%
  dplyr::select(Chr, Locus, Location, GENEID, Symbol, Strand, Distance, Width, ratioAPM13, ratioAPM14, meanB1TCDD24h, directionB1TCDD24h, p.value) %>%
  filter(Symbol %in% Gene_List)%>%
  drop_na() %>%
  print()
TCDD_24_merged.df <- tbl_df(TCDD_24_merged.df)

TCDD_72_merged.df <- merge(TCDD_72_all.df,meth_filtered.df,by="Locus")
TCDD_72_merged.df <- TCDD_72_merged.df %>%
  rename(Chr.x = "Chr") %>%
  rename(LOCATION = "Location") %>%
  rename(symbol = "Symbol") %>%
  dplyr::select(Chr, Locus, Location, GENEID, Symbol, Strand, Distance, Width, ratioAPM15, ratioAPM16, meanB1TCDD72h, directionB1TCDD72h, p.value) %>%
  filter(Symbol %in% Gene_List)%>%
  drop_na() %>%
  print()
TCDD_72_merged.df <- tbl_df(TCDD_72_merged.df)

TCDD_96_merged.df <- merge(TCDD_96_all.df,meth_filtered.df,by="Locus")
TCDD_96_merged.df <- TCDD_96_merged.df %>%
  rename(Chr.x = "Chr") %>%
  rename(LOCATION = "Location") %>%
  rename(symbol = "Symbol") %>%
  dplyr::select(Chr, Locus, Location, GENEID, Symbol, Strand, Distance, Width, ratioAPM17, ratioAPM18,  meanB1TCDD96h, directionB1TCDD96h, p.value) %>%
  filter(Symbol %in% Gene_List)%>%
  drop_na() %>%
  print()
TCDD_96_merged.df <- tbl_df(TCDD_96_merged.df)


#' 
#' 
#' # Finish adding relevant variables for Complex Heatmap
## -----------------------------------------------------------------------------------------------------------
# Matrix of methylation values
meth.df <- Reduce(merge, list(ES_merged.df,DMSO_24_merged.df, TCDD_24_merged.df, TCDD_72_merged.df, TCDD_96_merged.df))
meth.df <- tbl_df(meth.df)
meth.df <- meth.df %>%
  group_by(Symbol) %>% 
  filter(row_number() == 1) %>% #Slc2a17 had same p-value in duplicates but I kept only the one with greatest change in methylation in the TCDD76h timepoint.
  ungroup() %>%
  arrange(Symbol) %>%
  print()

df_meth <- meth.df %>%
  filter(Symbol %in% expr_filtered.df$Symbol) %>%
  dplyr::select(ratioAPM11, ratioAPM12, ratioAPM19, ratioAPM20, ratioAPM13, ratioAPM14, ratioAPM15, ratioAPM16, ratioAPM17, ratioAPM18) %>%
  rename(ratioAPM11="APM11", ratioAPM12="APM12", ratioAPM19="APM19", ratioAPM20="APM20", ratioAPM13="APM13", ratioAPM14="APM14", ratioAPM15="APM15", ratioAPM16="APM16", ratioAPM17="APM17", ratioAPM18="APM18") %>%
  print()
mat_meth <- as.matrix(df_meth)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------
# Matrix of expression values
mat_expr <- expr_filtered.df %>%
  ungroup() %>%
  arrange(Symbol) %>%
  print()

#Ensure that order of genes are the same as in the methylation dataframe
mat_expr <- mat_expr[ order(match(mat_expr$Symbol, meth.df$Symbol)), ]
mat_expr <- tbl_df(mat_expr)

Genes <- mat_expr$Symbol
write.csv(Genes, "./Output/Gene_List.csv")

mat_expr <- expr_filtered.df %>%
  dplyr::select(-Symbol) %>%
  print()

mat_expr <- as.matrix(mat_expr)


#' 
#' # Make separate heatmaps for comparing methylation and RNA-seq for each timepoint
#' 
## -----------------------------------------------------------------------------------------------------------
# Set up different matrices, then order the columns by DMSO first, then the tested cells

mat_ES_meth <- mat_meth[, c(3:4, 1:2)]
mat_24_meth <- mat_meth[, c(3:4, 5:6)]
mat_72_meth <- mat_meth[, c(3:4, 7:8)]
mat_96_meth <- mat_meth[, c(3:4, 9:10)]

mat_ES_expr <- mat_expr[, c(3:4, 1:2)]
mat_24_expr <- mat_expr[, c(3:4, 5:6)]
mat_72_expr <- mat_expr[, c(3:4, 7:8)]
mat_96_expr <- mat_expr[, c(3:4, 9:10)]

type_ES <- c("B1DMSO24", "B1DMSO24", "ES", "ES")
type_24 <- c("B1DMSO24", "B1DMSO24", "B1TCDD24", "B1TCDD24")
type_72 <- c("B1DMSO24", "B1DMSO24", "B1TCDD72", "B1TCDD72")
type_96 <- c("B1DMSO24", "B1DMSO24", "B1TCDD96", "B1TCDD96")


#' 
## -----------------------------------------------------------------------------------------------------------
#Reset the column names of each matrix to the identity of the actual sample and the replicate
colnames(mat_ES_meth) <- c("DMSO_24_1", "DMSO_24_2", "ES_1", "ES_2")
colnames(mat_ES_expr) <- c("DMSO_24_1", "DMSO_24_2", "ES_1", "ES_2")

colnames(mat_24_meth) <- c("DMSO_24_1", "DMSO_24_2", "TCDD_24_1", "TCDD_24_2")
colnames(mat_24_expr) <- c("DMSO_24_1", "DMSO_24_2", "TCDD_24_1", "TCDD_24_2")

colnames(mat_72_meth) <- c("DMSO_24_1", "DMSO_24_2", "TCDD_72_1", "TCDD_72_2")
colnames(mat_72_expr) <- c("DMSO_24_1", "DMSO_24_2", "TCDD_72_1", "TCDD_72_2")

colnames(mat_96_meth) <- c("DMSO_24_1", "DMSO_24_2", "TCDD_96_1", "TCDD_96_2")
colnames(mat_96_expr) <- c("DMSO_24_1", "DMSO_24_2", "TCDD_96_1", "TCDD_96_2")

#' 
#' 
#' 
#' 
## -----------------------------------------------------------------------------------------------------------
# generate directions for methylation
direction = rowMeans(mat_meth[, c(1:2,5:10)]) - rowMeans(mat_meth[, 3:4])
direction = ifelse(direction > 0, "hyper", "hypo") # average direction across all timepoints (samples)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------
# generate directions for methylation
direction_96 = rowMeans(mat_96_meth[, 3:4]) - rowMeans(mat_96_meth[, 1:2])
direction_96 = ifelse(direction_96 > 0, "hyper", "hypo") # for direction of methylation for TCDD96


#' 
## -----------------------------------------------------------------------------------------------------------

dist <- c(meth.df$Distance) # distance from DMRs to TSS of the assiciated genes. Positive values means the TSS is downstream of the CpG and negative values means the TSS is upstream.


#' 
## -----------------------------------------------------------------------------------------------------------
# matrix for correlation between methylation and expression
cor_pvalue = -log10(sapply(seq_len(nrow(mat_meth)), function(i) {
    cor.test(mat_meth[i, ], mat_expr[i, ])$p.value
}))

#' 
## -----------------------------------------------------------------------------------------------------------
# annotation to genes
location = as.character(meth.df$Location) # annotation to the gene models/ Location

#' 
#' 
#' 
#' # Build the Complex Heatmap Using the Variables Created
#' 
#' The different sources of information and corresponding variables are:
#' 
#' type: the label which shows whether the sample is TCDD or Control.
#' mat_meth: a matrix in which rows correspond to differetially methylated regions (DMRs). The value in the matrix is the mean methylation level in the DMR in every sample.
#' mat_expr: a matrix in which rows correspond to genes which are associated to the DMRs (i.e. the nearest gene to the DMR). The value in the matrix is the expression level for each gene in each sample. Expression is scaled for every gene across samples.
#' direction: direction of the methylation change (hyper meaning higher methylation in treatment samples, hypo means lower methylation in DMSO samples).
#' cor_pvalue: p-value for the correlation test between methylation and expression of the associated gene.
#' anno_gene: annotation to the gene models (location).
#' dist: distance from DMRs to TSS of the assiciated genes.
#' 
#' The clustering of columns for the methylation matrix are calculated first so that columns in the expression matrix can be adjusted to have the same column order as in the methylation matrix.
#' 
## -----------------------------------------------------------------------------------------------------------
#Reset the column names of each matrix to the identity of the actual sample and the replicate
colnames(mat_meth) <- c("ES_1", "ES_2", "DMSO_24_1", "DMSO_24_2", "TCDD_24_1", "TCDD_24_2", "TCDD_72_1", "TCDD_72_2", "TCDD_96_1", "TCDD_96_2")
colnames(mat_expr) <- c("ES_1", "ES_2", "DMSO_24_1", "DMSO_24_2", "TCDD_24_1", "TCDD_24_2", "TCDD_72_1", "TCDD_72_2", "TCDD_96_1", "TCDD_96_2")

#' 
#' 
#' 
## -----------------------------------------------------------------------------------------------------------
library(RColorBrewer)
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
direction_col = c("hyper" = "red", "hypo" = "blue")
expr_col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")) # -1, 0, 1 for original
pvalue_col_fun = colorRamp2(c(0, 2, 4), c("white", "white", "red"))
location_col = structure(brewer.pal(length(unique(location)), "Set1"), 
    names = unique(location))
dist_col_fun = colorRamp2(c(-550000, 200000), c("black", "white"))

#' 
#' We first define two column annotations and then make the complex heatmaps.
#' 
## -----------------------------------------------------------------------------------------------------------
ht_opt(
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 8), 
    heatmap_column_names_gp = gpar(fontsize = 8),
    heatmap_column_title_gp = gpar(fontsize = 10),
    heatmap_row_title_gp = gpar(fontsize = 8)
)

ha = HeatmapAnnotation(type = type, 
    col = list(type = c("ES" = "pink", "B1DMSO24" = "royalblue", "B1TCDD24" = "purple", "B1TCDD72" = "orange", "B1TCDD96" = "green")),
    annotation_name_side = "left")
ha2 = HeatmapAnnotation(type = type, 
    col = list(type = c("ES" = "pink", "B1DMSO24" = "royalblue", "B1TCDD24" = "purple", "B1TCDD72" = "orange", "B1TCDD96" = "green")), 
    show_legend = FALSE)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------
pdf("./Output/Figures/Composite_Heatmap_All.pdf", width = 15, height = 15)
ht_list = Heatmap(mat_meth, name = "methylation", col = meth_col_fun, cluster_columns = FALSE, show_column_names = TRUE,
    top_annotation = ha, column_title = "Methylation") +
    #rowAnnotation(Gene = anno_text(Genes), gp = gpar(fontsize = 0.1)) #+
    Heatmap(direction_96, name = "direction", col = direction_col) +
    Heatmap(mat_expr, name = "expression", 
        col = expr_col_fun, 
        cluster_columns = FALSE, show_column_names = TRUE,
        top_annotation = ha2, column_title = "Expression") +
    Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun) +
    Heatmap(location, name = "location", col = location_col) +
    Heatmap(dist, name = "dist_tss", col = dist_col_fun)

My_Heatmap <- draw(ht_list, row_km = 2, row_split = direction_96,
    #column_title = "Comprehensive correspondence between methylation, expression and other genomic features", 
    #column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
    merge_legends = TRUE, heatmap_legend_side = "bottom")

dev.off()

#' 
## -----------------------------------------------------------------------------------------------------------
ht_opt(RESET = TRUE)

#' 
#' 
#' ## Make the separate heatmaps comparing each timepoint
#' 
#' ### ES
#' 
## -----------------------------------------------------------------------------------------------------------
# generate directions for methylation
direction_ES = rowMeans(mat_ES_meth[, 3:4]) - rowMeans(mat_ES_meth[, 1:2])
direction_ES = ifelse(direction_ES > 0, "hyper", "hypo") # for direction of methylation for ES

ht_opt(
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 8), 
    heatmap_column_names_gp = gpar(fontsize = 8),
    heatmap_column_title_gp = gpar(fontsize = 10),
    heatmap_row_title_gp = gpar(fontsize = 8)
)

ha = HeatmapAnnotation(type = type_ES, 
    col = list(type = c("B1DMSO24" = "royalblue", "ES" = "pink")),
    annotation_name_side = "left")
ha2 = HeatmapAnnotation(type = type_ES, 
    col = list(type = c("B1DMSO24" = "royalblue", "ES" = "pink")), 
    show_legend = FALSE)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------
pdf("./Output/Figures/Heatmap_ESvsDMSO24.pdf", width = 15, height = 15)
ht_list = Heatmap(mat_ES_meth, name = "methylation", col = meth_col_fun, cluster_columns = FALSE, show_column_names = TRUE,
    top_annotation = ha, column_title = "Methylation") +
    #rowAnnotation(Gene = anno_text(Genes), gp = gpar(fontsize = 0.1)) #+
    Heatmap(direction_ES, name = "direction", col = direction_col) +
    Heatmap(mat_ES_expr, name = "expression", 
        col = expr_col_fun, 
        cluster_columns = FALSE, show_column_names = TRUE,
        top_annotation = ha2, column_title = "Expression") +
    Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun) +
    Heatmap(location, name = "location", col = location_col) +
    Heatmap(dist, name = "dist_tss", col = dist_col_fun)

My_Heatmap <- draw(ht_list, row_km = 2, row_split = direction_ES,
    #column_title = "Comprehensive correspondence between methylation, expression and other genomic features", 
    #column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
    merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

#' 
## -----------------------------------------------------------------------------------------------------------
ht_opt(RESET = TRUE)

#' 
#' ### TCDD24
#' 
## -----------------------------------------------------------------------------------------------------------
# generate directions for methylation
direction_24 = rowMeans(mat_24_meth[, 3:4]) - rowMeans(mat_24_meth[, 1:2])
direction_24 = ifelse(direction_24 > 0, "hyper", "hypo") # for direction of methylation for TCDD24

ht_opt(
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 8), 
    heatmap_column_names_gp = gpar(fontsize = 8),
    heatmap_column_title_gp = gpar(fontsize = 10),
    heatmap_row_title_gp = gpar(fontsize = 8)
)

ha = HeatmapAnnotation(type = type_24, 
    col = list(type = c("B1DMSO24" = "royalblue", "B1TCDD24" = "purple")),
    annotation_name_side = "left")
ha2 = HeatmapAnnotation(type = type_24, 
    col = list(type = c("B1DMSO24" = "royalblue", "B1TCDD24" = "purple")), 
    show_legend = FALSE)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------
pdf("./Output/Figures/Heatmap_TCDD24vsDMSO24.pdf", width = 15, height = 15)
ht_list = Heatmap(mat_24_meth, name = "methylation", col = meth_col_fun, cluster_columns = FALSE, show_column_names = TRUE,
    top_annotation = ha, column_title = "Methylation") +
    #rowAnnotation(Gene = anno_text(Genes), gp = gpar(fontsize = 0.1)) #+
    Heatmap(direction_24, name = "direction", col = direction_col) +
    Heatmap(mat_24_expr, name = "expression", 
        col = expr_col_fun, 
        cluster_columns = FALSE, show_column_names = TRUE,
        top_annotation = ha2, column_title = "Expression") +
    Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun) +
    Heatmap(location, name = "location", col = location_col) +
    Heatmap(dist, name = "dist_tss", col = dist_col_fun)

My_Heatmap <- draw(ht_list, row_km = 2, row_split = direction_24,
    #column_title = "Comprehensive correspondence between methylation, expression and other genomic features", 
    #column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
    merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

#' 
## -----------------------------------------------------------------------------------------------------------
ht_opt(RESET = TRUE)

#' 
#' ### TCDD72
#' 
## -----------------------------------------------------------------------------------------------------------
# generate directions for methylation
direction_72 = rowMeans(mat_24_meth[, 3:4]) - rowMeans(mat_72_meth[, 1:2])
direction_72 = ifelse(direction_72 > 0, "hyper", "hypo") # for direction of methylation for TCDD24

ht_opt(
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 8), 
    heatmap_column_names_gp = gpar(fontsize = 8),
    heatmap_column_title_gp = gpar(fontsize = 10),
    heatmap_row_title_gp = gpar(fontsize = 8)
)

ha = HeatmapAnnotation(type = type_72, 
    col = list(type = c("B1DMSO24" = "royalblue","B1TCDD72" = "orange")),
    annotation_name_side = "left")
ha2 = HeatmapAnnotation(type = type_72, 
    col = list(type = c("B1DMSO24" = "royalblue", "B1TCDD72" = "orange")), 
    show_legend = FALSE)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------
pdf("./Output/Figures/Heatmap_TCDD72vsDMSO24.pdf", width = 15, height = 15)
ht_list = Heatmap(mat_72_meth, name = "methylation", col = meth_col_fun, cluster_columns = FALSE, show_column_names = TRUE,
    top_annotation = ha, column_title = "Methylation") +
    #rowAnnotation(Gene = anno_text(Genes), gp = gpar(fontsize = 0.1)) #+
    Heatmap(direction_72, name = "direction", col = direction_col) +
    Heatmap(mat_72_expr, name = "expression", 
        col = expr_col_fun, 
        cluster_columns = FALSE, show_column_names = TRUE,
        top_annotation = ha2, column_title = "Expression") +
    Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun) +
    Heatmap(location, name = "location", col = location_col) +
    Heatmap(dist, name = "dist_tss", col = dist_col_fun)

My_Heatmap <- draw(ht_list, row_km = 2, row_split = direction_72,
    #column_title = "Comprehensive correspondence between methylation, expression and other genomic features", 
    #column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
    merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

#' 
## -----------------------------------------------------------------------------------------------------------
ht_opt(RESET = TRUE)

#' 
#' ### TCDD96
#' 
## -----------------------------------------------------------------------------------------------------------
# generate directions for methylation
direction_96 = rowMeans(mat_96_meth[, 3:4]) - rowMeans(mat_96_meth[, 1:2])
direction_96 = ifelse(direction_96 > 0, "hyper", "hypo") # for direction of methylation for TCDD96

ht_opt(
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 8), 
    heatmap_column_names_gp = gpar(fontsize = 8),
    heatmap_column_title_gp = gpar(fontsize = 10),
    heatmap_row_title_gp = gpar(fontsize = 8)
)

ha = HeatmapAnnotation(type = type_96, 
    col = list(type = c("B1DMSO24" = "royalblue", "B1TCDD96" = "green")),
    annotation_name_side = "left")
ha2 = HeatmapAnnotation(type = type_96, 
    col = list(type = c("B1DMSO24" = "royalblue", "B1TCDD96" = "green")), 
    show_legend = FALSE)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------
pdf("./Output/Figures/Heatmap_TCDD96vsDMSO24.pdf", width = 15, height = 15)
ht_list = Heatmap(mat_96_meth, name = "methylation", col = meth_col_fun, cluster_columns = FALSE, show_column_names = TRUE,
    top_annotation = ha, column_title = "Methylation") +
    #rowAnnotation(Gene = anno_text(Genes), gp = gpar(fontsize = 0.1)) #+
    Heatmap(direction_96, name = "direction", col = direction_col) +
    Heatmap(mat_96_expr, name = "expression", 
        col = expr_col_fun, 
        cluster_columns = FALSE, show_column_names = TRUE,
        top_annotation = ha2, column_title = "Expression") +
    Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun) +
    Heatmap(location, name = "location", col = location_col) +
    Heatmap(dist, name = "dist_tss", col = dist_col_fun)

My_Heatmap <- draw(ht_list, row_km = 2, row_split = direction_96,
    #column_title = "Comprehensive correspondence between methylation, expression and other genomic features", 
    #column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
    merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

#' 
## -----------------------------------------------------------------------------------------------------------
ht_opt(RESET = TRUE)

#' 
#' # Obtain List of Genes with Significant Correlation between Methylation and Expression
#' 
## -----------------------------------------------------------------------------------------------------------
# Create Correlation table with correlation p-values
Expression_ES <- rowMeans(mat_expr[, 1:2])
Expression_DMSO <- rowMeans(mat_expr[, 3:4])
Expression_24 <- rowMeans(mat_expr[, 5:6])
Expression_72 <- rowMeans(mat_expr[, 7:8])
Expression_96 <- rowMeans(mat_expr[, 9:10])

Correlation_ES.df <- meth.df %>%
  dplyr::select(Symbol, meanES, meanB1DMSO24h) %>%
  mutate(Beta = meanES - meanB1DMSO24h) %>%
  mutate(Expression = Expression_ES) %>%
  mutate(cor_pvalue = cor_pvalue) %>%
  rename(Symbol = "Gene") %>%
  dplyr::select(Gene, Beta, Expression, cor_pvalue)

Corr_Filtered_ES <- Correlation_ES.df %>%
  filter(cor_pvalue > 1.301) %>% #111 genes
  arrange(desc(cor_pvalue)) %>%
  print()

Correlation_24.df <- meth.df %>%
  dplyr::select(Symbol, meanB1TCDD24h, meanB1DMSO24h) %>%
  mutate(Beta = meanB1TCDD24h - meanB1DMSO24h) %>%
  mutate(Expression = Expression_DMSO) %>%
  mutate(cor_pvalue = cor_pvalue) %>%
  rename(Symbol = "Gene") %>%
  dplyr::select(Gene, Beta, Expression, cor_pvalue)

Corr_Filtered_24 <- Correlation_24.df %>%
  filter(cor_pvalue > 1.301) %>% #111 genes
  arrange(desc(cor_pvalue)) %>%
  print()

Correlation_72.df <- meth.df %>%
  dplyr::select(Symbol, meanB1TCDD72h, meanB1DMSO24h) %>%
  mutate(Beta = meanB1TCDD72h - meanB1DMSO24h) %>%
  mutate(Expression = Expression_72) %>%
  mutate(cor_pvalue = cor_pvalue) %>%
  rename(Symbol = "Gene") %>%
  dplyr::select(Gene, Beta, Expression, cor_pvalue)

Corr_Filtered_72 <- Correlation_72.df %>% #111 genes
  filter(cor_pvalue > 1.301) %>%
  arrange(desc(cor_pvalue)) %>%
  print()

Correlation_96.df <- meth.df %>%
  dplyr::select(Symbol, meanB1TCDD96h, meanB1DMSO24h) %>%
  mutate(Beta = meanB1TCDD96h - meanB1DMSO24h) %>%
  mutate(Expression = Expression_96) %>%
  mutate(cor_pvalue = cor_pvalue) %>% #111 genes
  rename(Symbol = "Gene") %>%
  dplyr::select(Gene, Beta, Expression, cor_pvalue)

Corr_Filtered_96 <- Correlation_96.df %>%
  filter(cor_pvalue > 1.301) %>% # 111 genes
  arrange(desc(cor_pvalue)) %>%
  print()


#' # Create new filtered dataframes for each of the groups, adding in the location column
#' 
## -----------------------------------------------------------------------------------------------------------

# Set up dataframes
Methylation_Location_ES.df <- Corr_Filtered_ES
Methylation_Location_24.df <- Corr_Filtered_24
Methylation_Location_72.df <- Corr_Filtered_72
Methylation_Location_96.df <- Corr_Filtered_96

#' 
## -----------------------------------------------------------------------------------------------------------
# Create matching "key" of genes with the location column

Location_key.df <- meth.df %>%
  rename(Symbol = "Gene") %>%
  dplyr::select(Gene, Location, Distance) %>%
  rename("Distance" = "Distance (bp)") %>%
  print()

## -----------------------------------------------------------------------------------------------------------
# Filter genes in key by only the genes with significant correlation
Location_key_filtered.df <- Location_key.df %>%
  filter(Gene %in% Corr_Filtered_96$Gene) %>%
  print()

#' 
## -----------------------------------------------------------------------------------------------------------
# Match each of the groups with the key to add in a location column
Filtered_ES_Location.df <- merge(Methylation_Location_ES.df,Location_key_filtered.df,by="Gene")
write.csv(Filtered_ES_Location.df, "./Output/Corr_ES.csv")
Filtered_24_Location.df <- merge(Methylation_Location_24.df, Location_key_filtered.df, by="Gene")
write.csv(Filtered_24_Location.df, "./Output/Corr_24.csv")
Filtered_72_Location.df <- merge(Methylation_Location_72.df,Location_key_filtered.df,by="Gene")
write.csv(Filtered_72_Location.df, "./Output/Corr_72.csv")
Filtered_96_Location.df <- merge(Methylation_Location_96.df, Location_key_filtered.df, by="Gene")
write.csv(Filtered_96_Location.df, "./Output/Corr_96.csv")


#' 
## -----------------------------------------------------------------------------------------------------------
# Classify each dataframe by location
## ES
Filtered_ES_Promoter.df <- Filtered_ES_Location.df %>%
  filter(Location == "promoter") %>%
  print()

Filtered_ES_Intron.df <- Filtered_ES_Location.df %>%
  filter(Location == "intron") %>%
  print()

Filtered_ES_coding.df <- Filtered_ES_Location.df %>%
  filter(Location == "coding") %>%
  print()

Filtered_ES_threeUTR.df <- Filtered_ES_Location.df %>%
  filter(Location == "threeUTR") %>%
  print()
## 24
Filtered_24_Promoter.df <- Filtered_24_Location.df %>%
  filter(Location == "promoter") %>%
  print()

Filtered_24_Intron.df <- Filtered_24_Location.df %>%
  filter(Location == "intron") %>%
  print()

Filtered_24_coding.df <- Filtered_24_Location.df %>%
  filter(Location == "coding") %>%
  print()

Filtered_24_threeUTR.df <- Filtered_24_Location.df %>%
  filter(Location == "threeUTR") %>%
  print()
## 72
Filtered_72_Promoter.df <- Filtered_72_Location.df %>%
  filter(Location == "promoter") %>%
  print()

Filtered_72_Intron.df <- Filtered_72_Location.df %>%
  filter(Location == "intron") %>%
  print()

Filtered_72_coding.df <- Filtered_72_Location.df %>%
  filter(Location == "coding") %>%
  print()

Filtered_72_threeUTR.df <- Filtered_72_Location.df %>%
  filter(Location == "threeUTR") %>%
  print()

##96
Filtered_96_Promoter.df <- Filtered_96_Location.df %>%
  filter(Location == "promoter") %>%
  print()

Filtered_96_Intron.df <- Filtered_96_Location.df %>%
  filter(Location == "intron") %>%
  print()

Filtered_96_coding.df <- Filtered_96_Location.df %>%
  filter(Location == "coding") %>%
  print()

Filtered_96_threeUTR.df <- Filtered_96_Location.df %>%
  filter(Location == "threeUTR") %>%
  print()


#' 
#' 
#' 
#' # Quantify percentage of hypo- and hyper-methylated, and downregulated and upregulated genes for each treatment and create bar graphs for this using ggplot2
#' 
## -----------------------------------------------------------------------------------------------------------
# Create summary dataframe

Methylation_Direction_ES.df <- data.frame("Number_hypo" = sum(Corr_Filtered_ES$Beta < 0), "Number_hyper" = sum(Corr_Filtered_ES$Beta > 0), "Number_up" = sum(Corr_Filtered_ES$Expression > 0), "Number_down" = sum(Corr_Filtered_ES$Expression < 0))
  

#' 
## -----------------------------------------------------------------------------------------------------------
# Add columns for percentages
Methylation_Direction_ES.df <- tbl_df(Methylation_Direction_ES.df)
Methylation_Direction_ES.df <- Methylation_Direction_ES.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

## -----------------------------------------------------------------------------------------------------------
# Repeat for other groups

## 24 hours
Methylation_Direction_24.df <- data.frame("Number_hypo" = sum(Corr_Filtered_24$Beta < 0), "Number_hyper" = sum(Corr_Filtered_24$Beta > 0), "Number_up" = sum(Corr_Filtered_24$Expression > 0), "Number_down" = sum(Corr_Filtered_24$Expression < 0))

Methylation_Direction_24.df <- tbl_df(Methylation_Direction_24.df)
Methylation_Direction_24.df <- Methylation_Direction_24.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

## 72 hours
Methylation_Direction_72.df <- data.frame("Number_hypo" = sum(Corr_Filtered_72$Beta < 0), "Number_hyper" = sum(Corr_Filtered_72$Beta > 0), "Number_up" = sum(Corr_Filtered_72$Expression > 0), "Number_down" = sum(Corr_Filtered_72$Expression < 0))

Methylation_Direction_72.df <- tbl_df(Methylation_Direction_72.df)
Methylation_Direction_72.df <- Methylation_Direction_72.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()
## 96 hours
Methylation_Direction_96.df <- data.frame("Number_hypo" = sum(Corr_Filtered_96$Beta < 0), "Number_hyper" = sum(Corr_Filtered_96$Beta > 0), "Number_up" = sum(Corr_Filtered_96$Expression > 0), "Number_down" = sum(Corr_Filtered_96$Expression < 0))

Methylation_Direction_96.df <- tbl_df(Methylation_Direction_96.df)
Methylation_Direction_96.df <- Methylation_Direction_96.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()


#' 
## -----------------------------------------------------------------------------------------------------------
# Now do the same for each location

## ES
Methylation_Promoter_ES.df <- data.frame("Number_hypo" = sum(Filtered_ES_Promoter.df$Beta < 0), "Number_hyper" = sum(Filtered_ES_Promoter.df$Beta > 0), "Number_up" = sum(Filtered_ES_Promoter.df$Expression > 0), "Number_down" = sum(Filtered_ES_Promoter.df$Expression < 0))

Methylation_Promoter_ES.df <- tbl_df(Methylation_Promoter_ES.df)
Methylation_Promoter_ES.df <- Methylation_Promoter_ES.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_Intron_ES.df <- data.frame("Number_hypo" = sum(Filtered_ES_Intron.df$Beta < 0), "Number_hyper" = sum(Filtered_ES_Intron.df$Beta > 0), "Number_up" = sum(Filtered_ES_Intron.df$Expression > 0), "Number_down" = sum(Filtered_ES_Intron.df$Expression < 0))

Methylation_Intron_ES.df <- tbl_df(Methylation_Intron_ES.df)
Methylation_Intron_ES.df <- Methylation_Intron_ES.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_coding_ES.df <- data.frame("Number_hypo" = sum(Filtered_ES_coding.df$Beta < 0), "Number_hyper" = sum(Filtered_ES_coding.df$Beta > 0), "Number_up" = sum(Filtered_ES_coding.df$Expression > 0), "Number_down" = sum(Filtered_ES_coding.df$Expression < 0))

Methylation_coding_ES.df <- tbl_df(Methylation_coding_ES.df)
Methylation_coding_ES.df <- Methylation_coding_ES.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_threeUTR_ES.df <- data.frame("Number_hypo" = sum(Filtered_ES_threeUTR.df$Beta < 0), "Number_hyper" = sum(Filtered_ES_threeUTR.df$Beta > 0), "Number_up" = sum(Filtered_ES_threeUTR.df$Expression > 0), "Number_down" = sum(Filtered_ES_threeUTR.df$Expression < 0))

Methylation_threeUTR_ES.df <- tbl_df(Methylation_threeUTR_ES.df)
Methylation_threeUTR_ES.df <- Methylation_threeUTR_ES.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

##24
Methylation_Promoter_24.df <- data.frame("Number_hypo" = sum(Filtered_24_Promoter.df$Beta < 0), "Number_hyper" = sum(Filtered_24_Promoter.df$Beta > 0), "Number_up" = sum(Filtered_24_Promoter.df$Expression > 0), "Number_down" = sum(Filtered_24_Promoter.df$Expression < 0))

Methylation_Promoter_24.df <- tbl_df(Methylation_Promoter_24.df)
Methylation_Promoter_24.df <- Methylation_Promoter_24.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_Intron_24.df <- data.frame("Number_hypo" = sum(Filtered_24_Intron.df$Beta < 0), "Number_hyper" = sum(Filtered_24_Intron.df$Beta > 0), "Number_up" = sum(Filtered_24_Intron.df$Expression > 0), "Number_down" = sum(Filtered_24_Intron.df$Expression < 0))

Methylation_Intron_24.df <- tbl_df(Methylation_Intron_24.df)
Methylation_Intron_24.df <- Methylation_Intron_24.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_coding_24.df <- data.frame("Number_hypo" = sum(Filtered_24_coding.df$Beta < 0), "Number_hyper" = sum(Filtered_24_coding.df$Beta > 0), "Number_up" = sum(Filtered_24_coding.df$Expression > 0), "Number_down" = sum(Filtered_24_coding.df$Expression < 0))

Methylation_coding_24.df <- tbl_df(Methylation_coding_24.df)
Methylation_coding_24.df <- Methylation_coding_24.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_threeUTR_24.df <- data.frame("Number_hypo" = sum(Filtered_24_threeUTR.df$Beta < 0), "Number_hyper" = sum(Filtered_24_threeUTR.df$Beta > 0), "Number_up" = sum(Filtered_24_threeUTR.df$Expression > 0), "Number_down" = sum(Filtered_24_threeUTR.df$Expression < 0))

Methylation_threeUTR_24.df <- tbl_df(Methylation_threeUTR_24.df)
Methylation_threeUTR_24.df <- Methylation_threeUTR_24.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

##72
Methylation_Promoter_72.df <- data.frame("Number_hypo" = sum(Filtered_72_Promoter.df$Beta < 0), "Number_hyper" = sum(Filtered_72_Promoter.df$Beta > 0), "Number_up" = sum(Filtered_72_Promoter.df$Expression > 0), "Number_down" = sum(Filtered_72_Promoter.df$Expression < 0))

Methylation_Promoter_72.df <- tbl_df(Methylation_Promoter_72.df)
Methylation_Promoter_72.df <- Methylation_Promoter_72.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_Intron_72.df <- data.frame("Number_hypo" = sum(Filtered_72_Intron.df$Beta < 0), "Number_hyper" = sum(Filtered_72_Intron.df$Beta > 0), "Number_up" = sum(Filtered_72_Intron.df$Expression > 0), "Number_down" = sum(Filtered_72_Intron.df$Expression < 0))

Methylation_Intron_72.df <- tbl_df(Methylation_Intron_72.df)
Methylation_Intron_72.df <- Methylation_Intron_72.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_coding_72.df <- data.frame("Number_hypo" = sum(Filtered_72_coding.df$Beta < 0), "Number_hyper" = sum(Filtered_72_coding.df$Beta > 0), "Number_up" = sum(Filtered_72_coding.df$Expression > 0), "Number_down" = sum(Filtered_72_coding.df$Expression < 0))

Methylation_coding_72.df <- tbl_df(Methylation_coding_72.df)
Methylation_coding_72.df <- Methylation_coding_72.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_threeUTR_72.df <- data.frame("Number_hypo" = sum(Filtered_72_threeUTR.df$Beta < 0), "Number_hyper" = sum(Filtered_72_threeUTR.df$Beta > 0), "Number_up" = sum(Filtered_72_threeUTR.df$Expression > 0), "Number_down" = sum(Filtered_72_threeUTR.df$Expression < 0))

Methylation_threeUTR_72.df <- tbl_df(Methylation_threeUTR_72.df)
Methylation_threeUTR_72.df <- Methylation_threeUTR_72.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

##96
Methylation_Promoter_96.df <- data.frame("Number_hypo" = sum(Filtered_96_Promoter.df$Beta < 0), "Number_hyper" = sum(Filtered_96_Promoter.df$Beta > 0), "Number_up" = sum(Filtered_96_Promoter.df$Expression > 0), "Number_down" = sum(Filtered_96_Promoter.df$Expression < 0))

Methylation_Promoter_96.df <- tbl_df(Methylation_Promoter_96.df)
Methylation_Promoter_96.df <- Methylation_Promoter_96.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_Intron_96.df <- data.frame("Number_hypo" = sum(Filtered_96_Intron.df$Beta < 0), "Number_hyper" = sum(Filtered_96_Intron.df$Beta > 0), "Number_up" = sum(Filtered_96_Intron.df$Expression > 0), "Number_down" = sum(Filtered_96_Intron.df$Expression < 0))

Methylation_Intron_96.df <- tbl_df(Methylation_Intron_96.df)
Methylation_Intron_96.df <- Methylation_Intron_96.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_coding_96.df <- data.frame("Number_hypo" = sum(Filtered_96_coding.df$Beta < 0), "Number_hyper" = sum(Filtered_96_coding.df$Beta > 0), "Number_up" = sum(Filtered_96_coding.df$Expression > 0), "Number_down" = sum(Filtered_96_coding.df$Expression < 0))

Methylation_coding_96.df <- tbl_df(Methylation_coding_96.df)
Methylation_coding_96.df <- Methylation_coding_96.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()

Methylation_threeUTR_96.df <- data.frame("Number_hypo" = sum(Filtered_96_threeUTR.df$Beta < 0), "Number_hyper" = sum(Filtered_96_threeUTR.df$Beta > 0), "Number_up" = sum(Filtered_96_threeUTR.df$Expression > 0), "Number_down" = sum(Filtered_96_threeUTR.df$Expression < 0))

Methylation_threeUTR_96.df <- tbl_df(Methylation_threeUTR_96.df)
Methylation_threeUTR_96.df <- Methylation_threeUTR_96.df %>%
  mutate(Percent_hypo = Number_hypo/(Number_hypo+Number_hyper)*100, Percent_hyper = Number_hyper/(Number_hypo+Number_hyper)*100, Percent_up = Number_up/(Number_up+Number_down)*100, Percent_down = Number_down/(Number_up+Number_down)*100) %>%
  print()


#' 
#' 
## -----------------------------------------------------------------------------------------------------------
# Create vector of groups
Group <- c("ES", "24 hrs", "72 hrs", "96 hrs")

# Merge dataframes for each comparison
## All groups and locations
Methylation_Direction.df <- rbind(Methylation_Direction_ES.df, Methylation_Direction_24.df, Methylation_Direction_72.df, Methylation_Direction_96.df)

Methylation_Direction.df <- tbl_df(Methylation_Direction.df)

## Promoter
Methylation_Promoter.df <- rbind(Methylation_Promoter_ES.df, Methylation_Promoter_24.df, Methylation_Promoter_72.df, Methylation_Promoter_96.df)

Methylation_Promoter.df <- tbl_df(Methylation_Promoter.df)
## Intron
Methylation_Intron.df <- rbind(Methylation_Intron_ES.df, Methylation_Intron_24.df, Methylation_Intron_72.df, Methylation_Intron_96.df)

Methylation_Intron.df <- tbl_df(Methylation_Intron.df)

## Coding
Methylation_coding.df <- rbind(Methylation_coding_ES.df, Methylation_coding_24.df, Methylation_coding_72.df, Methylation_coding_96.df)

Methylation_coding.df <- tbl_df(Methylation_coding.df)
## threeUTR
Methylation_threeUTR.df <- rbind(Methylation_threeUTR_ES.df, Methylation_threeUTR_24.df, Methylation_threeUTR_72.df, Methylation_threeUTR_96.df)

Methylation_threeUTR.df <- tbl_df(Methylation_threeUTR.df)

# Add in Group column
## All groups
Methylation_Direction.df <- Methylation_Direction.df %>%
  mutate(Group = Group) %>%
  dplyr::select(Group, Percent_hypo, Percent_hyper, Percent_up, Percent_down) %>%
  print()

## Promoter
Methylation_Promoter.df <- Methylation_Promoter.df %>%
  mutate(Group = Group) %>%
  dplyr::select(Group, Percent_hypo, Percent_hyper, Percent_up, Percent_down) %>%
  print()
## Intron
Methylation_Intron.df <- Methylation_Intron.df %>%
  mutate(Group = Group) %>%
  dplyr::select(Group, Percent_hypo, Percent_hyper, Percent_up, Percent_down) %>%
  print()
## Coding
Methylation_coding.df <- Methylation_coding.df %>%
  mutate(Group = Group) %>%
  dplyr::select(Group, Percent_hypo, Percent_hyper, Percent_up, Percent_down) %>%
  print()
## threeUTR
Methylation_threeUTR.df <- Methylation_threeUTR.df %>%
  mutate(Group = Group) %>%
  dplyr::select(Group, Percent_hypo, Percent_hyper, Percent_up, Percent_down) %>%
  print()
# Organize data "tidier"
## All groups
Methylation_Direction.df <- Methylation_Direction.df %>%
  dplyr::rename(Hypo = "Percent_hypo", Hyper = "Percent_hyper", Up = "Percent_up", Down = "Percent_down") %>%
  gather(Direction, Percent, Hypo:Down, na.rm = FALSE, convert = FALSE) %>%
  print()

## Promoter
Methylation_Promoter.df <- Methylation_Promoter.df %>%
  dplyr::rename(Hypo = "Percent_hypo", Hyper = "Percent_hyper", Up = "Percent_up", Down = "Percent_down") %>%
  gather(Direction, Percent, Hypo:Down, na.rm = FALSE, convert = FALSE) %>%
  print()
## Intron
Methylation_Intron.df <- Methylation_Intron.df %>%
  dplyr::rename(Hypo = "Percent_hypo", Hyper = "Percent_hyper", Up = "Percent_up", Down = "Percent_down") %>%
  gather(Direction, Percent, Hypo:Down, na.rm = FALSE, convert = FALSE) %>%
  print()
## Coding
Methylation_coding.df <- Methylation_coding.df %>%
  dplyr::rename(Hypo = "Percent_hypo", Hyper = "Percent_hyper", Up = "Percent_up", Down = "Percent_down") %>%
  gather(Direction, Percent, Hypo:Down, na.rm = FALSE, convert = FALSE) %>%
  print()

## threeUTR
Methylation_threeUTR.df <- Methylation_threeUTR.df %>%
  dplyr::rename(Hypo = "Percent_hypo", Hyper = "Percent_hyper", Up = "Percent_up", Down = "Percent_down") %>%
  gather(Direction, Percent, Hypo:Down, na.rm = FALSE, convert = FALSE) %>%
  print()
#Turn your 'treatment' column into a character vector
## All groups
Methylation_Direction.df$Group <- as.character(Methylation_Direction.df$Group)

## Promoter
Methylation_Promoter.df$Group <- as.character(Methylation_Promoter.df$Group)
## Intron
Methylation_Intron.df$Group <- as.character(Methylation_Intron.df$Group)
## Coding
Methylation_coding.df$Group <- as.character(Methylation_coding.df$Group)
## threeUTR
Methylation_threeUTR.df$Group <- as.character(Methylation_threeUTR.df$Group)
#Then turn it back into a factor with the levels in the correct order
## All groups
Methylation_Direction.df$Group <- factor(Methylation_Direction.df$Group, levels=unique(Methylation_Direction.df$Group))

## Promoter
Methylation_Promoter.df$Group <- factor(Methylation_Promoter.df$Group, levels=unique(Methylation_Promoter.df$Group))
## Intron
Methylation_Intron.df$Group <- factor(Methylation_Intron.df$Group, levels=unique(Methylation_Intron.df$Group))
## Coding
Methylation_coding.df$Group <- factor(Methylation_coding.df$Group, levels=unique(Methylation_coding.df$Group))
## threeUTR
Methylation_threeUTR.df$Group <- factor(Methylation_threeUTR.df$Group, levels=unique(Methylation_threeUTR.df$Group))
#Turn your 'treatment' column into a character vector
## All groups
Methylation_Direction.df$Direction <- as.character(Methylation_Direction.df$Direction)

## Promoter
Methylation_Promoter.df$Direction <- as.character(Methylation_Promoter.df$Direction)
## Intron
Methylation_Intron.df$Direction <- as.character(Methylation_Intron.df$Direction)
## Coding
Methylation_coding.df$Direction <- as.character(Methylation_coding.df$Direction)
## threeUTR
Methylation_threeUTR.df$Direction <- as.character(Methylation_threeUTR.df$Direction)
#Then turn it back into a factor with the levels in the correct order
## All groups
Methylation_Direction.df$Direction <- factor(Methylation_Direction.df$Direction, levels=unique(Methylation_Direction.df$Direction))

## Promoter
Methylation_Promoter.df$Direction <- factor(Methylation_Promoter.df$Direction, levels=unique(Methylation_Promoter.df$Direction))
## Intron
Methylation_Intron.df$Direction <- factor(Methylation_Intron.df$Direction, levels=unique(Methylation_Intron.df$Direction))

## Coding
Methylation_coding.df$Direction <- factor(Methylation_coding.df$Direction, levels=unique(Methylation_coding.df$Direction))
## threeUTR
Methylation_threeUTR.df$Direction <- factor(Methylation_threeUTR.df$Direction, levels=unique(Methylation_threeUTR.df$Direction))


#' 
#' 
#' # Create bar graphs using ggplot
## -----------------------------------------------------------------------------------------------------------
# All locations
plot.Methylation_Direction <- ggplot(Methylation_Direction.df, aes(x = Group, y = Percent, fill = Direction))

plot.Methylation_Direction +
  geom_col(position = "dodge", color = "black") +
  ylab(label = "Percentage of genes")+
  xlab(label= "Time of TCDD treatment")+
  ggtitle(label = "Direction of change among all genomic regions")+
theme(plot.title = element_text(hjust = 0.5, face="bold", size=18)) +
  theme(axis.text.x = element_text(face="bold", angle = 0,hjust = 1,vjust = 0.5, size = 18, colour = "black"), axis.text.y = element_text(face="bold", size=18, colour = "black"), axis.title.y=element_text(face="bold", size=18), axis.title.x=element_text(face="bold", size=18), legend.text=element_text(face="bold", size=18), legend.title=element_text(face="bold", size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(size = 2, colour = "black"))
ggsave("./Output/Figures/Methylation_Expression_Percents.png", dpi = 300)

## -----------------------------------------------------------------------------------------------------------
# Promoter
plot.Methylation_Promoter <- ggplot(Methylation_Promoter.df, aes(x = Group, y = Percent, fill = Direction))

plot.Methylation_Promoter +
  geom_col(position = "dodge", color = "black") +
  ylab(label = "Percentage of genes")+
  xlab(label= "Time of TCDD treatment")+
  ggtitle(label = "Direction of change among promoters")+
theme(plot.title = element_text(hjust = 0.5, face="bold", size=18)) +
  theme(axis.text.x = element_text(face="bold", angle = 0,hjust = 1,vjust = 0.5, size = 18, colour = "black"), axis.text.y = element_text(face="bold", size=18, colour = "black"), axis.title.y=element_text(face="bold", size=18), axis.title.x=element_text(face="bold", size=18), legend.text=element_text(face="bold", size=18), legend.title=element_text(face="bold", size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(size = 2, colour = "black"))
ggsave("./Output/Figures/Methylation_Promoter_Percents.png", dpi = 300)

#Intron
plot.Methylation_Intron <- ggplot(Methylation_Intron.df, aes(x = Group, y = Percent, fill = Direction))

plot.Methylation_Intron +
  geom_col(position = "dodge", color = "black") +
  ylab(label = "Percentage of genes")+
  xlab(label= "Time of TCDD treatment")+
  ggtitle(label = "Direction of change among introns")+
theme(plot.title = element_text(hjust = 0.5, face="bold", size=18)) +
  theme(axis.text.x = element_text(face="bold", angle = 0,hjust = 1,vjust = 0.5, size = 18, colour = "black"), axis.text.y = element_text(face="bold", size=18, colour = "black"), axis.title.y=element_text(face="bold", size=18), axis.title.x=element_text(face="bold", size=18), legend.text=element_text(face="bold", size=18), legend.title=element_text(face="bold", size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(size = 2, colour = "black"))
ggsave("./Output/Figures/Methylation_Intron_Percents.png", dpi = 300)

# Coding
plot.Methylation_coding <- ggplot(Methylation_coding.df, aes(x = Group, y = Percent, fill = Direction))

plot.Methylation_coding +
  geom_col(position = "dodge", color = "black") +
  ylab(label = "Percentage of genes")+
  xlab(label= "Time of TCDD treatment")+
  ggtitle(label = "Direction of change among coding regions")+
theme(plot.title = element_text(hjust = 0.5, face="bold", size=18)) +
  theme(axis.text.x = element_text(face="bold", angle = 0,hjust = 1,vjust = 0.5, size = 18, colour = "black"), axis.text.y = element_text(face="bold", size=18, colour = "black"), axis.title.y=element_text(face="bold", size=18), axis.title.x=element_text(face="bold", size=18), legend.text=element_text(face="bold", size=18), legend.title=element_text(face="bold", size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(size = 2, colour = "black"))
ggsave("./Output/Figures/Methylation_Coding_Percents.png", dpi = 300)

#threeUTR
plot.Methylation_threeUTR <- ggplot(Methylation_threeUTR.df, aes(x = Group, y = Percent, fill = Direction))
plot.Methylation_threeUTR +
  geom_col(position = "dodge", color = "black") +
  ylab(label = "Percentage of genes")+
  xlab(label= "Time of TCDD treatment")+
  ggtitle(label = "Direction of change among 3'UTR regions")+
theme(plot.title = element_text(hjust = 0.5, face="bold", size=18)) +
  theme(axis.text.x = element_text(face="bold", angle = 0,hjust = 1,vjust = 0.5, size = 18, colour = "black"), axis.text.y = element_text(face="bold", size=18, colour = "black"), axis.title.y=element_text(face="bold", size=18), axis.title.x=element_text(face="bold", size=18), legend.text=element_text(face="bold", size=18), legend.title=element_text(face="bold", size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(size = 2, colour = "black"))
ggsave("./Output/Figures/Methylation_threeUTR_Percents.png", dpi = 300)

#' 
#' 
#' 
#' # Make list of all genes affeced by TCDD treatment
#' 
## -----------------------------------------------------------------------------------------------------------
# Make list of all genes affected at 24, 72, and 96 hours of TCDD treatment for input into IPA analysis.
TCDD_genes <- c(Corr_Filtered_24$Gene, Corr_Filtered_72$Gene, Corr_Filtered_96$Gene)
# Remove duplicate genes
TCDD_genes <- unique(TCDD_genes) #112 genes total

# Also check for common genes affected
TCDD_Common <- intersect(Corr_Filtered_24$Gene, Corr_Filtered_72$Gene)
TCDD_Common <- intersect(Corr_Filtered_24$Gene, Corr_Filtered_96$Gene) # same 112 genes affected at all 3 timepoints!


#' # Finally, make gene table for IPA Analysis
#' 
## -----------------------------------------------------------------------------------------------------------
# Obtain the original filtered expression table
Expression_Genes.df <- read.csv("./Output/TCDD_B1_FDR_Filtered.csv")
rownames(Expression_Genes.df) <- NULL
tbl_df(Expression_Genes.df)

# Filter only for your genes in your newly acquired Gene List
Final_Genes.df <- Expression_Genes.df %>%
  filter(symbol %in% TCDD_Common) %>%
  arrange(FDR) %>%
  print()
write.csv(Final_Genes.df, "./Output/Final_Gene_List.csv")


#' 
#' 
#' # Make Correlation Plot
## -----------------------------------------------------------------------------------------------------------
# Make equation to get R squared
fit <- lm(Expression ~ Beta, data = Corr_Filtered_96)
summary(fit)


#' 
## -----------------------------------------------------------------------------------------------------------
equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(coef(x)[2], digits = 2),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));                 
}
p <- ggplot(Corr_Filtered_96, aes(Beta, Expression))
p + geom_point(aes(colour = Gene), size = 4)+
  theme(legend.position = "none") +
  geom_smooth(method = 'lm')+
  xlim(-0, 0.75) +
  ylim(-0.5, 0.4) +
  ylab(label = "Expression") +
  xlab(label = "Difference in Methylation After 96 Hours TCDD Treatment")+
  annotate("rect", xmin = 0., xmax = 0.6, ymin = 0.25, ymax = 0.4, fill="white", colour="red") +
  annotate("text", x = 0.3, y = 0.35, label = equation(fit), parse = TRUE)+
  annotate("text", x=0.3, y=0.3, label = "p = 0.03") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(filename = "./Output/Correlation.png", dpi = 300)
# Calculating Pearson's product-moment correlation
 cor.test(Corr_Filtered_96$Beta, Corr_Filtered_96$Expression, method = "pearson", conf.level = 0.95)

#' 
#' 
#' 
#' # Create Supplementary Table for the 239 genes both differentially expressed and methylated
#' 
## -----------------------------------------------------------------------------------------------------------
DM_DE.df <- meth_filtered.df %>%
  rename("LOCATION" = "Location", "GENEID" = "Gene ID", "symbol" = "Symbol", "Distance" = "Distance (bp)", "Width" = "Width (bp)", "p.value" = "P Value") %>%
  dplyr::select(-"P Value") %>%
  print()

write.csv(DM_DE.df, "./Output/DM_DE.csv")
  

#' 
## -----------------------------------------------------------------------------------------------------------
