## Script for BioSpyder TempO-seq analysis and figure preparation for PHH manuscript
## Author: Muriel Heldring
## Date: 18 March, 2022

## Requirements:
## - A_BioSpyder analysis pipeline_Lib.R
## - TP53RegulatedGenes.R
## - Raw count data files
## - ncbi_dataset.tsv

rm(list=ls())
gc()

######### ######### ######### ######### ######### 
######### Install packages ## ######### #########
######### ######### ######### ######### ######### 


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# 
# install.packages("rlang")
# install.packages("tidyverse")
# install.packages("colorspace")
# install.packages("gridExtra")
# install.packages("data.table")
# install.packages("pheatmap")
# install.packages("reshape2")
# install.packages("compare")
# install.packages("PoiClaClu")
# install.packages("hexbin")
# install.packages("proj4")
# install.packages("ggalt")
# install.packages("PoiClaClu")
# install.packages("ggfortify")
# install.packages("Rtsne")
# install.packages("viridis")
# install.packages("ggbeeswarm")
# install.packages("PerformanceAnalytics")
# install.packages("psych")
# install.packages("VennDiagram")
# install.packages("gplots")

# BiocManager::install("vsn")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("AnnotationDbi")

###### Load libraries ######
#detach("package:rlang", unload=TRUE)
library(ggbeeswarm)
library(gprofiler2)
library(tidyverse)
library(colorspace)
library(gridExtra)
library(data.table)
library(DESeq2)
library(pheatmap)
library(reshape2)
library(compare)
library(readxl)
library(PoiClaClu)
library(hexbin)
library(ggfortify)
library(Rtsne)
library(viridis)
library(doParallel)
library(foreach)
library(boot)
library(RColorBrewer)
library(PerformanceAnalytics)
library("psych")
library(drc)
library(VennDiagram)
library(gplots)

######### ######### ######### ######### ######### 
######### Settings  ######### ######### #########
######### ######### ######### ######### ######### 

# Working directories #
setwd("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/")

DATE <- "20220411"

# Folder containing only raw count table files (csv or excel files)
CountTableDir <- "/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Input/CountTables/" 
# Store metadata file and raw count files here
inputDir <-  '/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Input/' 
#Give path to output dir
outputDir <- "/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/"
outputDir_GR <- '/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/Figures/'
PATH_TO_PROTDATA <- "/data/muriel/Projects/DDP/DataAnalysis/Dynamics/Output/RData/"

scriptDir <- "/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/RScripts"
source(file.path(scriptDir, "A_BioSpyder analysis pipeline_Lib.R"))

load("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/TP53RegulatedGenes.R")

###### INPUT USER ######

MetaName <- "meta_PHH-HepG2WT.txt"                        #use template and save as txt file
CountThreshold <- 100000                                  #Exclude samples lower than this number for library size
filterForLibSize <- T                                     #TRUE if you want to filter for library size. FALSE is recommended in combination with DESeq2 package
TreatmentFilter <- c("CPT","MEDIUM")               #Treatments in the TREATMENT column of interest. Leave empty vector c() if you want to use all treatments.
VariableQCs <- "TREATMENT"                                #Variable to check distribution Library size and counts (choose meanID, CELL_ID, TREATMENT, SAMPLE_ID, EXP_ID, TIMEPOINT or REPLICATE)
RemoveSamples <- c("NoCells")                             #Give list of SAMPLE_IDs for samples which you would like to remove for further analysis (for example MAQC, NoCells samples, outliers)
ControlSample_design <- c('MEDIUM_conc0_TP8')             #Sample name (meanID, same for replicates) of control for log2FC and p-value calculations
#DESIGN: TREATMENT'_conc'CONCENTRATION'_TP'TIMEPOINT
ControlSample_meanID <- c('HepG2_HepG2_WT_MEDIUM_conc0_TP8') #Sample name (meanID, same for replicates) of control for log2FC and p-value calculations
#meanID: EXP_ID'_'CELL_ID'_'TREATMENT'_conc'CONCENTRATION'_TP'TIMEPOINT

NormMethod <- "DESeq2"                                    #Fill in desired normalization method prior to Deseq2 DEG calculation: "CPM" or "DESeq2"
Filtering <- "padj"                                       #Fill in filtering method ("padj", "log2FC" or "padj&log2FC")
CheckGenes <- c("TP53", "BTG2", "MDM2", "ATM", 
                "CHEK1", "CHEK2", "CDKN1A")               #Check expression of individual genes among samples (Give gene symbols)
#You can list as many genes as you like, however graph may become too full
Geneset <- tp53_downstream                                #Fill in file name of specific gene list of interest (GeneSymbols, txt format). If not used, fill in NA.

COLS_OF_INTEREST <- c("CELL_ID","TREATMENT","CONCENTRATION",
                      "TIMEPOINT","REPLICATE","CMP_ID",
                      "EXP_ID","PLATE_ID","SEX","CANCER",
                      "AGE","SAMPLE_CODE")

## Define interesting genes ##
PROBES_OF_INTEREST <- c("TP53_7287","MDM2_23384", "MDM2_27496", "MDM2_27497", "CDKN1A_1219", "BTG2_13191")
COLS <- c("CELL_ID","REPLICATE")
GENES_OF_INTEREST <- c("TP53", "MDM2","CDKN1A","BTG2")
PROBES_OF_INTEREST_EXT <- c("TP53_7287","MDM2_23384", "MDM2_27496", "MDM2_27497","MDM2_sum","MDM2_mean","CDKN1A_1219","BTG2_13191") # "ATM_16849","CHEK1_1288","CHEK2_1289",
PROTS_OF_INTEREST <- c("TP53","MDM2","CDKN1A","BTG2")

# DGE analysis parameters
custom <- TRUE
Threshold_padj <- 0.01                                     #Filtering of genes with p-adj lower than threshold
Threshold_log2FC <- 0.1                                   #Filtering of genes with absolute log2FC higher than threshold (genes lower than -0.5, and genes higher than 0.5)


######### ######### ######### ######### ######### 
######### Functions ######### ######### #########
######### ######### ######### ######### ######### 


# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)}

# Do correlation test and ouput p.value
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j])))
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

cor.test.p.xy <- function(x,y){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(y), 
    Vectorize(function(i,j) FUN(x[,i], y[,j])))
  dimnames(z) <- list(colnames(x), colnames(y))
  z
}

plotCor <- function(cormat, cormat.pval, title = "",labelSize = 12, textSize = 3, full = F){
  # Make a full or half full matrix
  if (full){
    melted_cormat <- reshape2::melt(cormat)
    melted_cormat.pval <- reshape2::melt(cormat.pval)}
  else {
    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- reshape2::melt(upper_tri)
    upper_tri.pval <- get_upper_tri(cormat.pval)
    melted_cormat.pval <- reshape2::melt(upper_tri.pval)}
  
  significant <- c()
  for (i in seq(1,length(melted_cormat$value))) {
    x <- paste0(as.character(melted_cormat$value[i]),"\n","p=",as.character(melted_cormat.pval$value[i]),
                ifelse(melted_cormat.pval$value[i] < 0.001,"***",
                       ifelse(melted_cormat.pval$value[i] < 0.01,"**",
                              ifelse(melted_cormat.pval$value[i] < 0.1,"*",
                                     ""))))
    #print(x)
    significant <- c(significant,x)}
  significant[significant=="NA\np=NANA"] <- NA
  significant[significant=="1\np=0***"] <- "1"
  melted_cormat$significant <- significant
  
  gplot <- ggplot()+
    geom_tile(data = melted_cormat, aes(Var2, Var1, fill = value),color = "white") +
    scale_fill_gradient2(low = plasma(8)[1], high = plasma(8)[5], mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name=expression(rho)) +
    theme_minimal()+ 
    coord_fixed() + 
    geom_text(data = melted_cormat, aes(Var2, Var1, label = significant), 
              color = "black", size = textSize) +
    ggtitle(title) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, 
                                 size = labelSize, hjust = 1),
      axis.text.y = element_text(size = labelSize),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5,size = (labelSize+4)),
      #title = element_text(size = (labelSize+4)),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(0, 1),
      legend.direction = "vertical")+
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                 title.position = "top", title.hjust = 0.5))
  gplot
}

jaccard <- function(a,b){
  numerator <- length(intersect(a,b))
  denominator <- length(union(a,b))
  return(numerator/denominator)
}

# Bootstrap function
correlation <- function(data, indices, probes, decims = 2,  
                        rows_select = "TP53_7287",
                        cols_select = c("TP53_7287","MDM2_sum","CDKN1A_1219","BTG2_13191")) {
  d <- data[indices,probes] # allows boot to select sample
  cormat <- round(cor(d),decims)
  return(cormat[rows_select,cols_select])
}

plotBootstrap <- function(results, columns, colnames){
  df <- data.frame(results$t[,columns])
  colnames(df) <- colnames
  df <- pivot_longer(df, cols = colnames(df))
  df <- df %>% group_by(name) %>% summarise(mean = mean(value), sd = sd(value))
  df$name <- factor(df$name, levels = colnames)
  return(df)}

removeProbeID <- function(gene_probe) {
  sapply(gene_probe, function(probe) {
    geneName <- strsplit(probe,"_")[[1]][1]})
}

save_pheatmap_pdf <- function(x, filename, width=6, height=9) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

######### ######### ######### ######### ######### 
######### Create output directory  #### #########
######### ######### ######### ######### ######### 


# Create dirs
dir.create(paste0(outputDir_GR,DATE,"_TempOSeq_Figs"))


######### ######### ######### ######### ######### 
######### Prefiltering data # ######### #########
######### ######### ######### ######### ######### 


######  Load data & library size filter ###### 
meta <- read.delim(paste0(inputDir, MetaName), stringsAsFactors = FALSE)
count <- fread(paste0(CountTableDir, "count_PHH-HepGWT.txt"))
count <- as.data.frame(count)
rownames(count) <- count[,1]
count$V1 <- NULL

## Excluding specific samples for further analysis

meta <- meta[meta$ERROR_CODE == 0,]

meta <- meta[which(!grepl(paste0(RemoveSamples, collapse = "|"), meta$CELL_ID)),]

count <- count[, rownames(meta)]

# Make extra columns for the design
meta$meanID <- paste0(meta$EXP_ID, "_",
                      meta$CELL_ID, "_",
                      meta$TREATMENT, "_conc",
                      meta$CONCENTRATION,"_TP", 
                      meta$TIMEPOINT)

meta$SAMPLENAME <- paste0(meta$EXP_ID, "_",
                          meta$CELL_ID, "_",
                          meta$TREATMENT, "_TP", 
                          meta$TIMEPOINT, "_conc",
                          meta$CONCENTRATION, "_rep",
                          meta$REPLICATE)

meta$DESIGN <- paste0(meta$TREATMENT, "_conc",
                      meta$CONCENTRATION,"_TP", 
                      meta$TIMEPOINT)

## Check ##

if(!all(rownames(meta) %in% colnames(count)) |
   !all(colnames(count) %in% rownames(meta))){
  warning("No identical sample names (meta & counts)")
  print(paste0("Unmatched samples in count table: ", paste0(c(colnames(count)[which(!colnames(count) %in% rownames(meta))]), collapse = ", ")))
  print(paste0("Unmatched samples in meta: ", paste0(c(rownames(meta)[which(!rownames(meta) %in% colnames(count))]), collapse = ", ")))
} else {print("All sample names are identical")}

# Set the same order of the samples
count <- count[, rownames(meta)]

## Filter out treatments of interest

if (!is.null(TreatmentFilter)) {
  meta <- meta[meta$TREATMENT %in% TreatmentFilter,]
  
  count <- count[, rownames(meta)]
}

## Check whether columns of count and rows of meta are indeed in the right order

if (all(rownames(meta) == colnames(count))) {
  print("Columns of count and rows of meta are in the right order.")
} else {
  warning("Columns of count and rows of meta are ***NOT*** the right order!")
}

# Make appropriate labelling
labelingDF <- tibble(DESIGN = unique(meta$DESIGN),
                     LABEL = factor(c("0.1 uM CDDP, 8hr",
                                      "   1 uM CDDP, 8hr",
                                      "3.3 uM CDDP, 8hr",
                                      " 10 uM CDDP, 8hr",
                                      " 33 uM CDDP, 8hr",
                                      "100 uM CDDP, 8hr",
                                      "   0 uM CDDP, 8hr",
                                      "0.1 uM CDDP, 24hr",
                                      "   1 uM CDDP, 24hr",
                                      "3.3 uM CDDP, 24hr",
                                      " 10 uM CDDP, 24hr",
                                      " 33 uM CDDP, 24hr",
                                      "100 uM CDDP, 24hr",
                                      "   0 uM CDDP, 24hr"),
                                    levels = c("   0 uM CDDP, 8hr",
                                               "   0 uM CDDP, 24hr",
                                               "0.1 uM CDDP, 8hr",
                                               "0.1 uM CDDP, 24hr",
                                               "   1 uM CDDP, 8hr",
                                               "   1 uM CDDP, 24hr",
                                               "3.3 uM CDDP, 8hr",
                                               "3.3 uM CDDP, 24hr",
                                               " 10 uM CDDP, 8hr",
                                               " 10 uM CDDP, 24hr",
                                               " 33 uM CDDP, 8hr",
                                               " 33 uM CDDP, 24hr",
                                               "100 uM CDDP, 8hr",
                                               "100 uM CDDP, 24hr")))
meta_labelled <- left_join(meta,labelingDF, by = "DESIGN")

## Library size filter ##

meta_labelled$LIB_SIZE <- colSums(count, na.rm = TRUE)

# Get info on libsize
summary(meta_labelled$LIB_SIZE)
df <- data.frame(x = rownames(meta_labelled), y = meta_labelled[,"LIB_SIZE"])
dft <- transpose(df, make.names = "x")
barplot(height = sort(as.numeric(as.vector(dft[1,]))))
abline(h = CountThreshold, col = "red")

## Filter out samples with small library sizes ##

meta_labelled$CELL_ID <- relevel(as.factor(meta_labelled$CELL_ID), ref = "HepG2_WT")
ggplot(data = meta_labelled) + 
  geom_hline(yintercept = CountThreshold) + 
  geom_point(aes(x = CELL_ID, y = LIB_SIZE, color = LABEL), stat = "identity") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis(discrete = T, name = "Condition") +
  theme_classic() + 
  xlab("Sample ID") + ylab("Library size") +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 18))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS1C_LibrarySizeFiltering.pdf"),width = 10, height = 5)

if (filterForLibSize) {
  meta_LowLibSize <- meta[which(meta$LIB_SIZE < CountThreshold), ]
  CELLID_LowLibSize <- as.data.frame(table(meta_LowLibSize$CELL_ID))
  meta_Filtered <- meta[which(meta$LIB_SIZE > CountThreshold), ] 
  count_Filtered <- count[, rownames(meta_Filtered)]
  print("Filtered for lib size")
} else {
  meta_LowLibSize <- meta[which(meta$LIB_SIZE < CountThreshold), ]
  CELLID_LowLibSize <- as.data.frame(table(meta_LowLibSize$CELL_ID))
  meta_Filtered <- meta
  count_Filtered <- count
  print("Not filtering for lib size")
}

ggplot(data = meta_Filtered) + 
  geom_hline(yintercept = CountThreshold) + 
  geom_point(aes(x = CELL_ID, y = LIB_SIZE, color = DESIGN), stat = "identity") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis(discrete = T, name = "Condition") +
  theme_classic() + 
  xlab("Sample ID") + ylab("Library size") +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 18))

## Check if CountTable contains NAs ##

count_NA <- count_Filtered[rowSums(is.na(count_Filtered)) > 0, colSums(is.na(count_Filtered)) > 0]
if(nrow(count_NA) == 0){
  print("No NAs in counts table")
} else {
  warning("NAs present in counts table")
  print(table(count_NA)); print(dim(count_NA))
}

count_Filtered <- na.omit(count_Filtered)

## Number of genes with at least one zero count among samples
NrGenes_zero <- sum(apply(count_Filtered, MARGIN = 1, function(x) sum(x == 0)) > 0)
message(paste0(NrGenes_zero, "/", nrow(count_Filtered), " genes have at least one zero count among samples and will not be used for size factor calculation"))

## Number of genes that have zeros in more than 10% of the samples
#head(count_Filtered,c(100,10))
NrGenes_zero <- sum(apply(count_Filtered, MARGIN = 1, function(x) sum(x == 0)) > 0.1*dim(count_Filtered)[2])
NrGenes_zero_names <- rownames(count_Filtered)[!apply(count_Filtered, MARGIN = 1, function(x) sum(x == 0)) > 0.1*dim(count_Filtered)[2]]
count_Filtered <- count_Filtered[NrGenes_zero_names,]

## Number of genes with at least one zero count among samples
NrGenes_zero <- sum(apply(count_Filtered, MARGIN = 1, function(x) sum(x == 0)) > 0)
message(paste0(NrGenes_zero, "/", nrow(count_Filtered), " genes have at least one zero count among samples and will not be used for size factor calculation"))

## Remove samples that have not at least one measurement within each treatment group (design)
# Get the 50 PHH samples and HepG2 names
meta_Filtered$CELL_ID <- as.factor(meta_Filtered$CELL_ID)
samples <- levels(droplevels(meta_Filtered$CELL_ID))

## Make the CELL_ID a factor with HepG2 as reference
meta_Filtered$CELL_ID <- relevel(meta_Filtered$CELL_ID, ref = "HepG2_WT")
meta_Filtered$SAMPLE_CODE <- as.numeric(meta_Filtered$CELL_ID)

## Make the colors and names for the tSNE plots
scdf_tmp <- meta_Filtered[,c("DESIGN","SAMPLE_CODE")]
scdf_tmp <- unique(scdf_tmp)
scdf_tmp <- scdf_tmp[order(scdf_tmp$SAMPLE_CODE),]
scdf_tmp$VIRCOLOR <- viridis(length(unique(scdf_tmp$DESIGN)))

# Get a named vector with the codes
sampleCodes <- scdf_tmp$SAMPLE_CODE
names(sampleCodes) <- scdf_tmp$CELL_ID

# Get a named vector with the colors
colorsCodes <- scdf_tmp$VIRCOLOR
names(colorsCodes) <- scdf_tmp$DESIGN

# Make a code and color column
meta_Filtered$COLOR <- colorsCodes[meta_Filtered$DESIGN]

# Get PHH samples
phhSamples <- rownames(meta_Filtered)[!meta_Filtered$EXP_ID == "HepG2"]


# Remove low count genes
# Here, I choose to remove genes that have a total sum lower than the number of samples
omit <- rowSums(count_Filtered) == 0 #< dim(count_Filtered)[2]
omitGenes <- rownames(count_Filtered)[omit]
omitGenes
removed <- count_Filtered[omit,]
count_Filtered <- count_Filtered[!omit,]

## Number of genes with at least one zero count among samples

NrGenes_zero <- sum(apply(count_Filtered, MARGIN = 1, function(x) sum(x == 0)) > 0)
message(paste0(NrGenes_zero, "/", nrow(count_Filtered), " genes have at least one zero count among samples and will not be used for size factor calculation"))


######### ######### ######### ######### ######### 
######### Run DESeq ######### ######### #########
######### ######### ######### ######### ######### 

## Make the control sample the first rows/columns ##

meta_Filtered$DESIGN <- as.factor(meta_Filtered$DESIGN)
meta_Filtered$DESIGN <- relevel(meta_Filtered$DESIGN, ref = ControlSample_design)
count_Filtered <- count_Filtered[, rownames(meta_Filtered)]

## Check whether design starts with controlSample and if columns of count and rows of meta are indeed in the right order

if (all(rownames(meta_Filtered) == colnames(count_Filtered)) &
    (levels(meta_Filtered$DESIGN)[1] == ControlSample_design)) {
  print("Columns of count and rows of meta are in the right order.")
} else {
  warning("Columns of count and rows of meta are ***NOT*** the right order!")
}

## Normalization & DEG calculation ##

dds <- DESeqDataSetFromMatrix(countData = count_Filtered,
                              colData = meta_Filtered,
                              design = ~ DESIGN)


if (NormMethod == "DESeq2"){
  dds <- estimateSizeFactors(dds)
} else if (NormMethod == "CPM"){
  sizeFactors(dds) <- colSums(count_Filtered_design)/1000000
} else {
  warning("Incorrect input for 'NormMethod', please double check above!")
}

## Run DESeq

Norm <- as.data.frame(counts(dds, normalized = TRUE)) #Normalized by size factors

Norm <- log2(Norm + 1) #Log2 normalization of size factor corrected counts

# Do deseq for PHH
meta_Filtered_phh <- meta_Filtered[which(meta_Filtered$EXP_ID == "PHH"), ] 
count_Filtered_phh <- count_Filtered[, rownames(meta_Filtered_phh)]
dds_phh <- DESeqDataSetFromMatrix(countData = count_Filtered_phh,
                                  colData = meta_Filtered_phh,
                                  design = ~ DESIGN)


if (NormMethod == "DESeq2"){
  dds_phh <- estimateSizeFactors(dds_phh)
} else if (NormMethod == "CPM"){
  sizeFactors(dds) <- colSums(count_Filtered_design)/1000000
} else {
  warning("Incorrect input for 'NormMethod', please double check above!")
}

## Run DESeq

Norm_phh <- as.data.frame(counts(dds_phh, normalized = TRUE)) #Normalized by size factors

Norm_phh <- log2(Norm_phh + 1) #Log2 normalization of size factor corrected counts

vsd_phh <- assay(vst(dds_phh, blind=FALSE))

# Optional: run DESeq
dds_deseq <- DESeq(dds_phh)

# Do deseq for HepG2
meta_Filtered_hepg2 <- meta_Filtered[which(meta_Filtered$EXP_ID == "HepG2"), ] 
count_Filtered_hepg2 <- count_Filtered[, rownames(meta_Filtered_hepg2)]
dds_hepg2 <- DESeqDataSetFromMatrix(countData = count_Filtered_hepg2,
                                    colData = meta_Filtered_hepg2,
                                    design = ~ DESIGN)


if (NormMethod == "DESeq2"){
  dds_hepg2 <- estimateSizeFactors(dds_hepg2)
} else if (NormMethod == "CPM"){
  sizeFactors(dds) <- colSums(count_Filtered_design)/1000000
} else {
  warning("Incorrect input for 'NormMethod', please double check above!")
}

## Run DESeq

Norm_hepg2 <- as.data.frame(counts(dds_hepg2, normalized = TRUE)) #Normalized by size factors

Norm_hepg2 <- log2(Norm_hepg2 + 1) #Log2 normalization of size factor corrected counts

vsd_hepg2 <- assay(vst(dds_hepg2, blind=FALSE))

# Optional: run DESeq
dds_deseq_hepg2 <- DESeq(dds_hepg2)

######### ######### ######### ######### ######### 
######### Make enrichment figures ##### #########
######### ######### ######### ######### ######### 

# Make custom background
custom_bg <- unique(removeProbeID(rownames(count_Filtered_phh)))

# Get the results per treatment group
res_0.1 <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc0.1_TP8",ControlSample_design),alpha = 0.05))
res_1 <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc1_TP8",ControlSample_design),alpha = 0.05))
res_3.3 <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc3.3_TP8",ControlSample_design),alpha = 0.05))
res_10 <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc10_TP8",ControlSample_design),alpha = 0.05))
res_33 <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc33_TP8",ControlSample_design),alpha = 0.05))
res_100 <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc100_TP8",ControlSample_design),alpha = 0.05))

# Remove the NAs that originated from independent filtering
res_0.1 <- res_0.1[!is.na(res_0.1$padj),]
res_1 <- res_1[!is.na(res_1$padj),]
res_3.3 <- res_3.3[!is.na(res_3.3$padj),]
res_10 <- res_10[!is.na(res_10$padj),]
res_33 <- res_33[!is.na(res_33$padj),]
res_100 <- res_100[!is.na(res_100$padj),]

# Save the gene lists
deg_0.1 <- rownames(res_0.1[res_0.1$padj < Threshold_padj & (res_0.1$log2FoldChange) > Threshold_log2FC,])
deg_1 <- rownames(res_1[res_1$padj < Threshold_padj & (res_1$log2FoldChange) > Threshold_log2FC,])
deg_3.3 <- rownames(res_3.3[res_3.3$padj < Threshold_padj & (res_3.3$log2FoldChange) > Threshold_log2FC,])
deg_10 <- rownames(res_10[res_10$padj < Threshold_padj & (res_10$log2FoldChange) > Threshold_log2FC,])
deg_33 <- rownames(res_33[res_33$padj < Threshold_padj & (res_33$log2FoldChange) > Threshold_log2FC,])
deg_100 <- rownames(res_100[res_100$padj < Threshold_padj & (res_100$log2FoldChange) > Threshold_log2FC,])

# Remove probe IDs
deg_0.1 <- unique(removeProbeID(deg_0.1))
deg_1 <- unique(removeProbeID(deg_1))
deg_3.3 <- unique(removeProbeID(deg_3.3))
deg_10 <- unique(removeProbeID(deg_10))
deg_33 <- unique(removeProbeID(deg_33))
deg_100 <- unique(removeProbeID(deg_100))

# Run g:Profiler on all gene sets to search for cell death enrichment
if (custom){
  GEA <- gost(query = list("0.1" = deg_0.1,"1" = deg_1,
                           "3.3" = deg_3.3,"10" = deg_10,
                           "33" = deg_33,"100" = deg_100), 
              organism = "hsapiens", ordered_query = FALSE, 
              multi_query = FALSE, significant = F, exclude_iea = FALSE, 
              measure_underrepresentation = FALSE, evcodes = FALSE, 
              user_threshold = 0.05, correction_method = "fdr", 
              domain_scope = "annotated", custom_bg = custom_bg, 
              numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)
} else {
  GEA <- gost(query = list("0.1" = deg_0.1,"1" = deg_1,
                           "3.3" = deg_3.3,"10" = deg_10,
                           "33" = deg_33,"100" = deg_100), 
              organism = "hsapiens", ordered_query = FALSE, 
              multi_query = FALSE, significant = F, exclude_iea = FALSE, 
              measure_underrepresentation = FALSE, evcodes = FALSE, 
              user_threshold = 0.05, correction_method = "g_SCS", 
              domain_scope = "annotated", custom_bg = NULL, 
              numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)}

GEA_res <- GEA$result
GEA_res <- GEA_res %>% mutate(significant = factor(significant, levels = c("TRUE","FALSE")))
GEA_res_sign <- GEA_res %>% filter(significant == TRUE)

# Make supervised data frame for enrichment
terms_of_interest <- factor(c("cell death",
                              "regulation of cell death",
                              "regulation of programmed cell death",
                              "programmed cell death",
                              "necrotic cell death",
                              "apoptotic process",
                              "regulation of apoptotic process"))

sup <- GEA_res[GEA_res$term_name %in% terms_of_interest,]
sup <- sup %>% mutate(query = factor(query, levels = c("0.1","1","3.3","10","33","100")))
ggplot(sup) + 
  geom_bar(aes(x = query, y = recall, fill = term_name), stat = "identity", position = "dodge") +
  scale_fill_viridis(discrete = T, name = "GO term") +
  xlab(expression("Cisplatin concentration ("*mu*"M)")) + ylab("Recall") + 
  theme_classic()
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS4A_supervised_GEA_CellDeath.pdf"),width = 6, height = 3)

terms_of_interest <- factor(c("cellular response to DNA damage stimulus",
                              "DNA repair",
                              "DNA integrity checkpoint",
                              "DNA damage response",
                              "regulation of cell cycle",
                              "cell cycle checkpoint signaling"))
sup <- GEA_res[GEA_res$term_name %in% terms_of_interest,]
sup <- sup %>% mutate(query = factor(query, levels = c("0.1","1","3.3","10","33","100")))
ggplot(sup) +
  geom_bar(aes(x = query, y = recall, fill = term_name), stat = "identity", position = "dodge") +
  scale_fill_viridis(discrete = T, name = "GO term") +
  xlab(expression("Cisplatin concentration ("*mu*"M)")) + ylab("Recall") +
  theme_classic()
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS4B_supervised_GEA_CellSurvival.pdf"),width = 6, height = 3)

# Do unsupervised
# Make unsupervised figures for top 10 enrichment terms
# 0.1
unsup_0.1 <- GEA_res[GEA_res$query == "0.1",][1:10,]# & GEA_res$p_value < 0.1,]
unsup_0.1 <- unsup_0.1 %>% mutate(term_name = factor(term_name, levels = unique(term_name)),
                                  term_nameadj = gsub("(.{47})","\\1\n",term_name),
                                  term_nameadj = factor(term_nameadj,levels = term_nameadj))
ggplot(unsup_0.1) + 
  geom_bar(aes(y = term_nameadj, x = recall, fill = significant), stat = "identity") +
  scale_fill_manual("Significant\n(p < 0.05)", 
                    labels = c("TRUE" = "Yes", "FALSE" = "No"), 
                    values = c("FALSE" = "firebrick3","TRUE" = "darkgreen")) +
  xlab("Recall") + ylab("") + 
  ggtitle(expression("0.1 "*mu*"M cisplatin")) +
  theme_classic() +
  theme(text = element_text(size=18),
        title = element_text(size=18))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS4C_unsupervised_GEA_CellSurvival.pdf"),width = 8, height = 7)
# 1
unsup_1 <- GEA_res[GEA_res$query == "1",][1:10,]
unsup_1 <- unsup_1 %>% mutate(term_name = factor(term_name, levels = unique(term_name)),
                              term_nameadj = gsub("(.{50})","\\1\n",term_name),
                              term_nameadj = factor(term_nameadj,levels = term_nameadj))
unsup_1 <- na.omit(unsup_1)
ggplot(unsup_1) + 
  geom_bar(aes(y = term_nameadj, x = recall, fill = significant), stat = "identity") +
  scale_fill_manual("Significant\n(p < 0.05)", labels = c("TRUE" = "Yes", "FALSE" = "No"), values = c("TRUE" = "darkgreen","FALSE" = "firebrick3")) +
  xlab("Recall") + ylab("") + 
  ggtitle(expression("1 "*mu*"M cisplatin")) +
  theme_classic() +
  theme(text = element_text(size=14),
        title = element_text(size=18))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_unsupervised_GEA_CellSurvival.pdf"),width = 6, height = 7)
# 3.3
unsup_3.3 <- GEA_res[GEA_res$query == "3.3",][1:10,]
unsup_3.3 <- unsup_3.3 %>% mutate(term_name = factor(term_name, levels = unique(term_name)),
                                  term_nameadj = gsub("(.{35})","\\1\n",term_name),
                                  term_nameadj = factor(term_nameadj,levels = term_nameadj))
ggplot(unsup_3.3) + 
  geom_bar(aes(y = term_nameadj, x = recall, fill = significant), stat = "identity") +
  scale_fill_manual("Significant\n(p < 0.05)", labels = c("TRUE" = "Yes", "FALSE" = "No"), values = c("FALSE" = "firebrick3","TRUE" = "darkgreen")) +
  xlab("Recall") + ylab("") + 
  ggtitle(expression("3.3 "*mu*"M cisplatin")) +
  theme_classic() +
  theme(text = element_text(size=18),
        title = element_text(size=18))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS4D_unsupervised_GEA_CellSurvival.pdf"),width = 8, height = 7)
# 10
unsup_10 <- GEA_res[GEA_res$query == "10",][1:10,]
unsup_10 <- unsup_10 %>% mutate(term_name = factor(term_name, levels = unique(term_name)),
                                term_nameadj = gsub("(.{35})","\\1\n",term_name),
                                term_nameadj = factor(term_nameadj,levels = term_nameadj))
ggplot(unsup_10) + 
  geom_bar(aes(y = term_nameadj, x = recall, fill = significant), stat = "identity") +
  scale_fill_manual("Significant\n(p < 0.05)", labels = c("TRUE" = "Yes", "FALSE" = "No"), values = c("FALSE" = "firebrick3","TRUE" = "darkgreen")) +
  xlab("Recall") + ylab("") + 
  ggtitle(expression("10 "*mu*"M cisplatin")) +
  theme_classic() +
  theme(text = element_text(size=18),
        title = element_text(size=18))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS4E_unsupervised_GEA_CellSurvival.pdf"),width = 8, height = 7)
# 33
unsup_33 <- GEA_res[GEA_res$query == "33",][1:10,]
unsup_33 <- unsup_33 %>% mutate(term_name = factor(term_name, levels = unique(term_name)),
                                term_nameadj = gsub("(.{43})","\\1\n",term_name),
                                term_nameadj = factor(term_nameadj,levels = term_nameadj))
ggplot(unsup_33) + 
  geom_bar(aes(y = term_nameadj, x = recall, fill = significant), stat = "identity") +
  scale_fill_manual("Significant\n(p < 0.05)", labels = c("TRUE" = "Yes", "FALSE" = "No"), values = c("FALSE" = "firebrick3","TRUE" = "darkgreen")) +
  xlab("Recall") + ylab("") + 
  ggtitle(expression("33 "*mu*"M cisplatin")) +
  theme_classic() +
  theme(text = element_text(size=18),
        title = element_text(size=18))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS4F_unsupervised_GEA_CellSurvival.pdf"),width = 8, height = 7)
# 100
unsup_100 <- GEA_res[GEA_res$query == "100",][1:10,]
unsup_100 <- unsup_100 %>% mutate(term_name = factor(term_name, levels = unique(term_name)),
                                  term_nameadj = gsub("(.{40})","\\1\n",term_name),
                                  term_nameadj = factor(term_nameadj,levels = term_nameadj))
ggplot(unsup_100) + 
  geom_bar(aes(y = term_nameadj, x = recall, fill = significant), stat = "identity") +
  scale_fill_manual("Significant\n(p < 0.05)", labels = c("TRUE" = "Yes", "FALSE" = "No"), values = c("FALSE" = "firebrick3","TRUE" = "darkgreen")) +
  xlab("Recall") + ylab("") + 
  ggtitle(expression("100 "*mu*"M cisplatin")) +
  theme_classic() +
  theme(text = element_text(size=18),
        title = element_text(size=18))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS4G_unsupervised_GEA_CellSurvival.pdf"),width = 8, height = 7)

unsup <- bind_rows(unsup_0.1,
                   unsup_3.3,
                   unsup_10,
                   unsup_33,
                   unsup_100) %>% dplyr::select(c(query, significant, recall, term_name))
write_excel_csv(unsup, file = paste0(outputDir,"GEA_unsupervised.csv"), delim = "\t")

######### ######### ######### ######### ######### 
######### Make data for TXG-MAPr  ##### #########
######### ######### ######### ######### ######### 

## TXG-mapr data
# Get entrezIDs
entrez <- read_delim(paste0(inputDir,"ncbi_dataset.tsv"))
entrez <- unique(entrez %>% dplyr::select(c("NCBI GeneID","Symbol"))) %>% rename("Symbol" = "gene",
                                                                                 "NCBI GeneID" = "gene_id")
# Results PHH
res_3_8h_phh <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc3.3_TP8",ControlSample_design),alpha = 0.05))
res_3_8h_phh <- res_3_8h_phh[!is.na(res_3_8h_phh$padj),]
res_3_8h_phh$gene <- sapply(rownames(res_3_8h_phh),function(x){str_split(x,"_")[[1]][1]})
res_3_8h_phh$experiment_id <- "res_3_8h_phh"
res_3_8h_phh <- res_3_8h_phh %>% rename("log2FoldChange" = "log2_fc",
                                        "pvalue" = "p_value",
                                        "padj" = "p_adj")
res_3_8h_phh <- left_join(res_3_8h_phh,entrez, by = "gene") %>% 
  dplyr::select(c(experiment_id,gene_id,log2_fc,p_value,p_adj))

res_3_24h_phh <- data.frame(results(dds_deseq, contrast=c("DESIGN","CPT_conc3.3_TP24",ControlSample_design),alpha = 0.05))
res_3_24h_phh <- res_3_24h_phh[!is.na(res_3_24h_phh$padj),]
res_3_24h_phh$gene <- sapply(rownames(res_3_24h_phh),function(x){str_split(x,"_")[[1]][1]})
res_3_24h_phh$experiment_id <- "res_3_24h_phh"
res_3_24h_phh <- res_3_24h_phh %>% rename("log2FoldChange" = "log2_fc",
                                          "pvalue" = "p_value",
                                          "padj" = "p_adj")
res_3_24h_phh <- left_join(res_3_24h_phh,entrez, by = "gene") %>% 
  dplyr::select(c(experiment_id,gene_id,log2_fc,p_value,p_adj))

# Results HepG2
res_3_8h_hepg2 <- data.frame(results(dds_deseq_hepg2, contrast=c("DESIGN","CPT_conc3.3_TP8",ControlSample_design),alpha = 0.05))
res_3_8h_hepg2 <- res_3_8h_hepg2[!is.na(res_3_8h_hepg2$padj),]
res_3_8h_hepg2$gene <- sapply(rownames(res_3_8h_hepg2),function(x){str_split(x,"_")[[1]][1]})
res_3_8h_hepg2$experiment_id <- "res_3_8h_hepg2"
res_3_8h_hepg2 <- res_3_8h_hepg2 %>% rename("log2FoldChange" = "log2_fc",
                                            "pvalue" = "p_value",
                                            "padj" = "p_adj")
res_3_8h_hepg2 <- left_join(res_3_8h_hepg2,entrez, by = "gene") %>% 
  dplyr::select(c(experiment_id,gene_id,log2_fc,p_value,p_adj))

res_3_24h_hepg2 <- data.frame(results(dds_deseq_hepg2, contrast=c("DESIGN","CPT_conc3.3_TP24",ControlSample_design),alpha = 0.05))
res_3_24h_hepg2 <- res_3_24h_hepg2[!is.na(res_3_24h_hepg2$padj),]
res_3_24h_hepg2$gene <- sapply(rownames(res_3_24h_hepg2),function(x){str_split(x,"_")[[1]][1]})
res_3_24h_hepg2$experiment_id <- "res_3_24h_hepg2"
res_3_24h_hepg2 <- res_3_24h_hepg2 %>% rename("log2FoldChange" = "log2_fc",
                                              "pvalue" = "p_value",
                                              "padj" = "p_adj")
res_3_24h_hepg2 <- left_join(res_3_24h_hepg2,entrez, by = "gene") %>% 
  dplyr::select(c(experiment_id,gene_id,log2_fc,p_value,p_adj))

# Bind data frames
results <- rbind(res_3_8h_phh, res_3_24h_phh,
                 res_3_8h_hepg2, res_3_24h_hepg2)

write_delim(results, file = paste0(outputDir,"TXGmaprData/DESeq_results.txt"), delim = "\t")

######### ######### ######### ######### ######### 
######### Analyse data from TXG-MAPr  # #########
######### ######### ######### ######### ######### 

# Load data
txg_output <- read_csv(paste0(outputDir,"TXGmaprData/WGCNA_PHH_newEGS_all_2022-03-18 16_26_57.csv"))
# Rename column
txg_output <- txg_output %>% dplyr::rename("condition" = "uploaded condition")
# Get module name
txg_output <- txg_output %>% mutate(ModName = sapply(module, function(x){str_split(x, ":")[[1]][2]}),
                                    Experiment = ifelse(condition == "upload.res_3_24h_hepg2", "HepG2\n24 hr",
                                                        ifelse(condition == "upload.res_3_8h_hepg2", "HepG2\n8 hr",
                                                               ifelse(condition == "upload.res_3_24h_phh", "PHH\n24 hr",
                                                                      "PHH\n8 hr"))))

# Select modules
txg_subset <- txg_output %>% group_by(condition) %>% arrange(desc(EGS)) %>%
  slice_head(n = 20)
write_delim(txg_subset, file = paste0(outputDir,"TXGmaprData/High_EGS_Modules.txt"), delim = "\t")

display_venn(
  lapply(txg_subset %>% group_split(), function(x){x %>% pull(ModName)}),
  category.names = sapply(txg_subset %>% group_split(), function(x){unique(x %>% pull(Experiment))}),
  fill = viridis(4),
  width = 3750,
  cex = 2, 
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.25, 0.25, 0.15, 0.15)
)

venn.diagram(x = lapply(txg_subset %>% group_split(), function(x){x %>% pull(ModName)}),
             filename = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS1E_venndiagram.png"),
             category.names = sapply(txg_subset %>% group_split(), function(x){unique(x %>% pull(Experiment))}),
             fill = viridis(4),
             width = 3500, height = 3500,
             main.cex = 2,
             main.fontface = "bold",
             main.fontfamily = "sansserif",
             main = expression("3.3 "*mu*"M Cisplatin"),
             fontface = "bold",
             cex = 2, 
             cat.cex = 1.8,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.23, 0.23, 0.15, 0.15),
             fontfamily = "sansserif",
             cat.fontfamily = "sansserif")
             

######### ######### ######### ######### ######### 
###### Fit Hill equations in expression data ####
######### ######### ######### ######### ######### 

## Start correlation analysis ##

# Transpose Norm
Normt <- data.frame(t(as.matrix(Norm)))
Norm_phht <- data.frame(t(as.matrix(Norm_phh)))
vsd_phht <- data.frame(t(as.matrix(vsd_phh)))

# Make MDM2_sum and MDM2_mean
Normt <- Normt %>% rownames_to_column(var = "PlateWellID") %>% mutate(MDM2_sum = MDM2_23384 + MDM2_27496 + MDM2_27497,
                                                                      MDM2_mean = (MDM2_23384 + MDM2_27496 + MDM2_27497)/3)
Norm_phht <- Norm_phht %>% rownames_to_column(var = "PlateWellID") %>% mutate(MDM2_sum = MDM2_23384 + MDM2_27496 + MDM2_27497,
                                                                              MDM2_mean = (MDM2_23384 + MDM2_27496 + MDM2_27497)/3)
vsd_phht <- vsd_phht %>% rownames_to_column(var = "PlateWellID") %>% mutate(MDM2_sum = MDM2_23384 + MDM2_27496 + MDM2_27497,
                                                                            MDM2_mean = (MDM2_23384 + MDM2_27496 + MDM2_27497)/3)

# Hill EC50 data
meta_Filtered$PlateWellID <- rownames(meta_Filtered)
meta_Filtered <- meta_Filtered %>% mutate(TRTMT_GROUP = paste0(EXP_ID," ",TIMEPOINT," h"))

Hill <- left_join(meta_Filtered %>% dplyr::select(c(PlateWellID, CONCENTRATION, TRTMT_GROUP)),
                  Normt %>% dplyr::select(c(PlateWellID, MDM2_sum,TP53_7287,CDKN1A_1219,BTG2_13191)),
                  by = "PlateWellID") 

# Make colors
cols_hill <- viridis(4)


# TP53
DR.tp53 <- drm(TP53_7287 ~ CONCENTRATION, 
               data= Hill %>% filter(CONCENTRATION <= 10),
               curveid = TRTMT_GROUP,
               robust = 'mean', #non-robust least squares estimation ("mean")
               fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
Hill_output <- data.frame(gene = c("TP53"),
                          Type = c("PHH", "PHH", "HepG2", "HepG2"),
                          Time = c(8,24,8,24),
                          EC50 = c(DR.tp53$coefficients[13],NA,DR.tp53$coefficients[15:16]), 
                          Expression = predict(DR.tp53,data.frame(CONCENTRATION = DR.tp53$coefficients[13:16], 
                                                                  TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h"))))
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_Hill_fit_TP53.pdf"),width = 5, height = 5)
plot(DR.tp53, xlab = expression("Dose ("*mu*"M)"), ylab = "TP53 expression", col = cols_hill, lwd = 3, legend = F, pch = c(0,1,2,5), cex = 1.5,
     cex.lab = 1.5, cex.axis = 1.2)
points(x = DR.tp53$coefficients[13:16], y = predict(DR.tp53,data.frame(CONCENTRATION = DR.tp53$coefficients[13:16], 
                                                                       TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h"))), 
       col = cols_hill, pch = 20, cex = 2)
dev.off()

# MDM2
DR.mdm2 <- drm(MDM2_sum ~ CONCENTRATION, 
               data= Hill %>% filter((TRTMT_GROUP %in% c("PHH 8 h","PHH 24 h","HepG2 24 h") & CONCENTRATION <= 10) |
                                       (TRTMT_GROUP %in% c("HepG2 8 h") & CONCENTRATION <= 3.3)),
               curveid = TRTMT_GROUP,
               robust = 'mean', #non-robust least squares estimation ("mean")
               fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
Hill_output <- rbind(Hill_output,
                     data.frame(gene = c("MDM2"),
                                Type = c("PHH", "PHH", "HepG2", "HepG2"),
                                Time = c(8,24,8,24),
                                EC50 = DR.mdm2$coefficients[13:16], 
                                Expression = predict(DR.mdm2,data.frame(CONCENTRATION = DR.mdm2$coefficients[13:16], 
                                                                        TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h")))))
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3E_Hill_fit_MDM2.pdf"),width = 5, height = 5)
plot(DR.mdm2, xlab = expression("Dose ("*mu*"M)"), ylab = "MDM2 expression", col = cols_hill, lwd = 3, legend = F, pch = c(0,1,2,5), cex = 1.5,
     cex.lab = 1.5, cex.axis = 1.2)
points(x = DR.mdm2$coefficients[13:16], y = predict(DR.mdm2,data.frame(CONCENTRATION = DR.mdm2$coefficients[13:16], 
                                                                       TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h"))), 
       col = cols_hill, pch = 20, cex = 2)
dev.off()

# CDKN1A
DR.cdkn1a <- drm(CDKN1A_1219 ~ CONCENTRATION, 
                 data= Hill %>% filter(CONCENTRATION <= 10),
                 curveid = TRTMT_GROUP,
                 robust = 'mean', #non-robust least squares estimation ("mean")
                 fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
Hill_output <- rbind(Hill_output,
                     data.frame(gene = c("CDKN1A"),
                                Type = c("PHH", "PHH", "HepG2", "HepG2"),
                                Time = c(8,24,8,24),
                                EC50 = DR.cdkn1a$coefficients[13:16], 
                                Expression = predict(DR.cdkn1a,data.frame(CONCENTRATION = DR.cdkn1a$coefficients[13:16], 
                                                                          TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h")))))
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3F_Hill_fit_CDKN1A.pdf"),width = 5, height = 5)
plot(DR.cdkn1a, xlab = expression("Dose ("*mu*"M)"), ylab = "CDKN1A expression", col = cols_hill, lwd = 3, legend = F, pch = c(0,1,2,5), cex = 1.5,
     cex.lab = 1.5, cex.axis = 1.2)
points(x = DR.cdkn1a$coefficients[13:16], y = predict(DR.cdkn1a,data.frame(CONCENTRATION = DR.cdkn1a$coefficients[13:16], 
                                                                           TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h"))), 
       col = cols_hill, pch = 20, cex = 2)
dev.off()

# BTG2
DR.btg2 <- drm(BTG2_13191 ~ CONCENTRATION, 
               data= Hill %>% filter(CONCENTRATION <= 33),
               curveid = TRTMT_GROUP,
               robust = 'mean', #non-robust least squares estimation ("mean")
               fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
Hill_output <- rbind(Hill_output,
                     data.frame(gene = c("BTG2"),
                                Type = c("PHH", "PHH", "HepG2", "HepG2"),
                                Time = c(8,24,8,24),
                                EC50 = DR.btg2$coefficients[13:16], 
                                Expression = predict(DR.btg2,data.frame(CONCENTRATION = DR.btg2$coefficients[13:16], 
                                                                        TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h")))))
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3G_Hill_fit_BTG2.pdf"),width = 5, height = 5)
plot(DR.btg2, xlab = expression("Dose ("*mu*"M)"), ylab = "BTG2 expression", col = cols_hill, lwd = 3, legend = F, pch = c(0,1,2,5), cex = 1.5,
     cex.lab = 1.5, cex.axis = 1.2)
legend("topleft",
       legend = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h","EC50"),
       col = c(cols_hill,"black"),
       lty = c(1,2,3,4,NA), lwd = 3,
       pch = c(0,1,2,5,20))
points(x = DR.btg2$coefficients[13:16], y = predict(DR.btg2,data.frame(CONCENTRATION = DR.btg2$coefficients[13:16], 
                                                                       TRTMT_GROUP = c("PHH 8 h","PHH 24 h","HepG2 8 h","HepG2 24 h"))), 
       col = cols_hill, pch = 20, cex = 2)
dev.off()

write_delim(Hill_output, file = paste0(outputDir,"Table_EC50.txt"), delim = "\t")


# Add a column with the patient IDs
countsData_Norm <- left_join(Normt,meta_Filtered %>% 
                               dplyr::select(c("EXP_ID","CELL_ID","CONCENTRATION","TIMEPOINT","DESIGN","SAMPLE_CODE","PlateWellID")) ,by = "PlateWellID")
countsData_Norm_all <- left_join(Normt,meta_Filtered ,by = "PlateWellID")
countsData_Norm_phh <- left_join(Norm_phht,meta_Filtered_phh %>% rownames_to_column(var = "PlateWellID"),by = "PlateWellID")
countsData_vsd_phh <- left_join(vsd_phht,meta_Filtered_phh %>% rownames_to_column(var = "PlateWellID"),by = "PlateWellID")

# Check technical variation
tmp <- countsData_Norm_all %>% group_by(meanID) %>% mutate(TechRep = seq(1,n())) %>% ungroup()
cn <- colnames(tmp)
tmp_genes <- tmp[,0:(which(cn == "MDM2_sum")-1)]

# Get 10 highest and 10 lowest expressed genes
tmp_mean <- tmp_genes %>% summarise_if(is.numeric,mean) %>% c(., recursive=TRUE)
tmp_mean <- sort(tmp_mean)
tmp_low <- names(tmp_mean[0:10])
tmp_high <- names(tmp_mean[(length(tmp_mean)-10):length(tmp_mean)])

# Mean rho low
rho_low <- c()
for (i in seq(1:10)) {
  tmp_x <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = tmp_low[i]) %>% dplyr::select(-c(meanID, CELL_ID))
  rho <- c(cor.test(tmp_x$`Rep 1`, tmp_x$`Rep 2`)$estimate,
           cor.test(tmp_x$`Rep 1`, tmp_x$`Rep 3`)$estimate,
           cor.test(tmp_x$`Rep 2`, tmp_x$`Rep 3`)$estimate)
  meanrho <- mean(rho)
  rho_low <- c(rho_low,meanrho)
}
mean_rho_low <- mean(rho_low)
sd_rho_low <- sd(rho_low)
# Mean rho high
rho_high <- c()
for (i in seq(1:10)) {
  tmp_x <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = tmp_high[i]) %>% dplyr::select(-c(meanID, CELL_ID))
  rho <- c(cor.test(tmp_x$`Rep 1`, tmp_x$`Rep 2`)$estimate,
           cor.test(tmp_x$`Rep 1`, tmp_x$`Rep 3`)$estimate,
           cor.test(tmp_x$`Rep 2`, tmp_x$`Rep 3`)$estimate)
  meanrho <- mean(rho)
  rho_high <- c(rho_high,meanrho)
}
mean_rho_high <- mean(rho_high)
sd_rho_high <- sd(rho_high)
rho_tibble <- tibble(class = factor(c("Low", "High"),levels = c("Low", "High")), 
                     mean = c(mean_rho_low, mean_rho_high), 
                     sd = c(sd_rho_low, sd_rho_high))

# Make summary plot 
ggplot(rho_tibble) +
  geom_errorbar(aes(x = class, ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_point(aes(x = class, y = mean, color = class), size = 3, show.legend = FALSE) +
  theme(legend.position = "none") +
  scale_color_viridis_d() +
  xlab("Expression group") + ylab("Correlation coefficient") +
  ylim(0,1) + 
  theme_classic()
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3M_SummaryCorrelations.pdf"),width = 3, height = 2)

# Make plots for low and high expressed genes:

for (i in seq(1:10)) {
  # Low expressed genes
  pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_CorrelationMatrix_low_",tmp_low[i],".pdf"),width = 8, height = 8)
  tmp_x <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = tmp_low[i])
  pairs.panels(tmp_x %>% dplyr::select(-c(meanID, CELL_ID)), 
               hist.col="#7d7d7d", 
               method = "pearson",
               lm = T,
               show.points=TRUE, 
               stars=TRUE, 
               gap=0.05, 
               pch=20, 
               ellipses=FALSE, 
               scale=FALSE,
               jiggle=F,
               main=tmp_low[i], 
               cex.main=3,
               cex.axis=2.5,
               col="#c20202", 
               pty="m", 
               font=2)
  dev.off()
  
  # High expressed genes
  pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_CorrelationMatrix_high_",tmp_high[i],".pdf"),width = 8, height = 8)
  tmp_x <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = tmp_high[i])
  pairs.panels(tmp_x %>% dplyr::select(-c(meanID, CELL_ID)), 
               hist.col="#7d7d7d", 
               method = "pearson",
               lm = T,
               show.points=TRUE, 
               stars=TRUE, 
               gap=0.05, 
               pch=20, 
               ellipses=FALSE, 
               scale=FALSE,
               jiggle=F,
               main=tmp_high[i], 
               cex.main=3,
               cex.axis=2.5,
               col="#c20202", 
               pty="m", 
               font=2)
  dev.off()
  
}


# For TP53 
tmp_p53 <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = TP53_7287)
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3I_CorrelationMatrix_TP53.pdf"),width = 8, height = 8)
pairs.panels(tmp_p53 %>% dplyr::select(-c(meanID, CELL_ID)), 
             hist.col="#7d7d7d", 
             method = "pearson",
             lm = T,
             show.points=TRUE, 
             stars=TRUE, 
             gap=0.05, 
             pch=20, 
             ellipses=FALSE, 
             scale=FALSE,
             jiggle=F,
             main="TP53", 
             cex.main=3,
             cex.axis=2.5,
             col="#c20202", 
             pty="m", 
             font=2)
dev.off()

# For MDM2_mean
tmp_mdm2 <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = MDM2_mean)
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3J_CorrelationMatrix_MDM2.pdf"),width = 8, height = 8)
pairs.panels(tmp_mdm2 %>% dplyr::select(-c(meanID, CELL_ID)), 
             hist.col="#7d7d7d", 
             method = "pearson",
             lm = T,
             show.points=TRUE, 
             stars=TRUE, 
             gap=0.05, 
             pch=20, 
             ellipses=FALSE, 
             scale=FALSE,
             jiggle=F,
             main="MDM2", 
             cex.main=3,
             cex.axis=2.5,
             col="#c20202", 
             pty="m", 
             font=2)
dev.off()

# For CDKN1A
tmp_p21 <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = CDKN1A_1219)
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3K_CorrelationMatrix_CDKN1A.pdf"),width = 8, height = 8)
pairs.panels(tmp_p21 %>% dplyr::select(-c(meanID, CELL_ID)), 
             hist.col="#7d7d7d", 
             method = "pearson",
             lm = T,
             show.points=TRUE, 
             stars=TRUE, 
             gap=0.05, 
             pch=20, 
             ellipses=FALSE, 
             scale=FALSE,
             jiggle=F,
             main="CDKN1A", 
             cex.main=3,
             cex.axis=2.5,
             col="#c20202", 
             pty="m", 
             font=2)
dev.off()

# For BTG2 - no nice correlation
tmp_btg2 <- pivot_wider(tmp, id_cols = c(meanID, CELL_ID),names_from = TechRep, names_prefix = "Rep ", values_from = BTG2_13191)
pdf(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3L_CorrelationMatrix_BTG2.pdf"),width = 8, height = 8)
pairs.panels(tmp_btg2 %>% dplyr::select(-c(meanID, CELL_ID)), 
             hist.col="#7d7d7d", 
             method = "pearson",
             lm = T,
             show.points=TRUE, 
             stars=TRUE, 
             gap=0.05, 
             pch=20, 
             ellipses=FALSE, 
             scale=FALSE,
             jiggle=F,
             main="BTG2", 
             cex.main=3,
             cex.axis=2.5,
             col="#c20202", 
             pty="m", 
             font=2)
dev.off()

## Make boxplots, scatterplots and heatmaps ##

# Prepare summary data (summarise replicates)
### hmdata for all genes
VarData <- left_join(Normt,meta_Filtered %>% 
                       dplyr::select(c("EXP_ID","CELL_ID","CONCENTRATION","TIMEPOINT","DESIGN","SAMPLE_CODE","PlateWellID","COLOR","meanID")) ,by = "PlateWellID")
VarData <- VarData %>% dplyr::select(-c("MDM2_mean","MDM2_sum")) %>% as_tibble() # %>% filter(EXP_ID == "PHH") %>% select(-c("MDM2_mean","MDM2_sum")) %>% as_tibble() 
pcaData <- VarData %>% dplyr::select(-c("PlateWellID","CONCENTRATION","TIMEPOINT","SAMPLE_CODE","EXP_ID","CELL_ID","COLOR","DESIGN")) %>% 
  group_by(meanID) %>% summarise_all(list(mean)) 
pcaDataMeta <- VarData %>% dplyr::select(c("EXP_ID","CELL_ID","CONCENTRATION","TIMEPOINT","DESIGN","SAMPLE_CODE","COLOR","meanID")) %>% 
  distinct()

### Heatmap data for gene-wise heatmaps
hmdata <- countsData_Norm[,c("EXP_ID","CELL_ID","CONCENTRATION","TIMEPOINT","DESIGN","SAMPLE_CODE", PROBES_OF_INTEREST_EXT)]
hmdata <- hmdata %>% rename(TP53_7287 = "TP53", 
                            CDKN1A_1219 = "CDKN1A",
                            BTG2_13191 = "BTG2",
                            MDM2_mean = "MDM2")
hmdata$TIMEPOINT <- factor(hmdata$TIMEPOINT, levels = c("8","24"))
hmdata$CONCENTRATION <- factor(hmdata$CONCENTRATION, levels = c("0","0.1","1","3.3","10","33","100"))
hmdata$CONC_TIME <- paste0("c",hmdata$CONCENTRATION,"_tp",hmdata$TIMEPOINT)
hmdata$CONC_TIME <- factor(hmdata$CONC_TIME, levels = c("c0_tp8", "c0.1_tp8","c1_tp8","c3.3_tp8", "c10_tp8","c33_tp8", "c100_tp8",
                                                        "c0_tp24", "c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24"))

hmdataSum <- hmdata %>% ungroup() %>% 
  group_by(EXP_ID,DESIGN, CELL_ID,CONCENTRATION,TIMEPOINT,CONC_TIME,SAMPLE_CODE) %>% summarise_all(list(mean),na.rm = F) %>% 
  ungroup()

hmdataMeta <- unique(countsData_Norm_all[,c("EXP_ID","CELL_ID","SAMPLE_CODE","LIVERPATH","CANCER","SEX","CONFLUENCY","AGE",
                                            "DIABETES","SMOKING","ALCOHOL", "CARDIAC.DISEASE","HYPERTENSION")])
rownames(hmdataMeta) <- as.character(hmdataMeta$CELL_ID)
hmdataMeta <- hmdataMeta[as.character(sort(hmdataMeta$CELL_ID)),]

# Make a heatmap for all samples separately on TP53 expression
hmdataSumTP53 <- hmdataSum[,c("EXP_ID","CELL_ID","CONC_TIME", "TP53")]
hmdataSumTP53 <- hmdataSumTP53 %>% spread(key = CONC_TIME, value = TP53)
hmdataSumTP53 <- data.frame(hmdataSumTP53)
rownames(hmdataSumTP53) <- as.character(relevel(factor(hmdataSumTP53$CELL_ID), ref = "HepG2_WT"))
hmdataSumTP53 <- hmdataSumTP53[,c("CELL_ID","c0_tp8", "c0.1_tp8","c1_tp8","c3.3_tp8", "c10_tp8","c33_tp8", "c100_tp8",
                                  "c0_tp24", "c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

hmdata_tp53_all <- hmdataSumTP53[,c("c0_tp8","c0.1_tp8", "c1_tp8","c3.3_tp8", "c10_tp8", "c33_tp8", "c100_tp8",
                                    "c0_tp24","c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

# Make a heatmap for all samples separately on BTG2 expression
hmdataSumBTG2 <- hmdataSum[,c("EXP_ID","CELL_ID","CONC_TIME", "BTG2")]
hmdataSumBTG2 <- hmdataSumBTG2 %>% spread(key = CONC_TIME, value = BTG2)
hmdataSumBTG2 <- data.frame(hmdataSumBTG2)
rownames(hmdataSumBTG2) <- as.character(relevel(factor(hmdataSumBTG2$CELL_ID), ref = "HepG2_WT"))
hmdataSumBTG2 <- hmdataSumBTG2[,c("CELL_ID","c0_tp8", "c0.1_tp8","c1_tp8","c3.3_tp8", "c10_tp8","c33_tp8", "c100_tp8",
                                  "c0_tp24", "c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

hmdata_BTG2_all <- hmdataSumBTG2[,c("c0_tp8","c0.1_tp8", "c1_tp8","c3.3_tp8", "c10_tp8", "c33_tp8", "c100_tp8",
                                    "c0_tp24","c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

# Make a heatmap for all samples separately on MDM2 expression
hmdataSumMDM2 <- hmdataSum[,c("EXP_ID","CELL_ID","CONC_TIME", "MDM2")]
hmdataSumMDM2 <- hmdataSumMDM2 %>% spread(key = CONC_TIME, value = MDM2)
hmdataSumMDM2 <- data.frame(hmdataSumMDM2)
rownames(hmdataSumMDM2) <- as.character(relevel(factor(hmdataSumMDM2$CELL_ID), ref = "HepG2_WT"))
hmdataSumMDM2 <- hmdataSumMDM2[,c("CELL_ID","c0_tp8", "c0.1_tp8","c1_tp8","c3.3_tp8", "c10_tp8","c33_tp8", "c100_tp8",
                                  "c0_tp24", "c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

hmdata_MDM2_all <- hmdataSumMDM2[,c("c0_tp8","c0.1_tp8", "c1_tp8","c3.3_tp8", "c10_tp8", "c33_tp8", "c100_tp8",
                                    "c0_tp24","c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

# Make a heatmap for all samples separately on CDKN1A expression
hmdataSumCDKN1A <- hmdataSum[,c("EXP_ID","CELL_ID","CONC_TIME", "CDKN1A")]
hmdataSumCDKN1A <- hmdataSumCDKN1A %>% spread(key = CONC_TIME, value = CDKN1A)
hmdataSumCDKN1A <- data.frame(hmdataSumCDKN1A)
rownames(hmdataSumCDKN1A) <- as.character(relevel(factor(hmdataSumCDKN1A$CELL_ID), ref = "HepG2_WT"))
hmdataSumCDKN1A <- hmdataSumCDKN1A[,c("CELL_ID","c0_tp8", "c0.1_tp8","c1_tp8","c3.3_tp8", "c10_tp8","c33_tp8", "c100_tp8",
                                      "c0_tp24", "c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

hmdata_CDKN1A_all <- hmdataSumCDKN1A[,c("c0_tp8","c0.1_tp8", "c1_tp8","c3.3_tp8", "c10_tp8", "c33_tp8", "c100_tp8",
                                        "c0_tp24","c0.1_tp24","c1_tp24", "c3.3_tp24", "c10_tp24", "c33_tp24", "c100_tp24")]

## Make 3 groups with low, intermediate and high ##
# Clustering based on TP53 at all concentrations
d1 <- dist(hmdata_tp53_all)
hc1 <- hclust(d1)
plot(hc1)
groups_tp53_all <- cutree(hc1, k=3)
pheatmap(hmdata_tp53_all, cluster_rows = T, cluster_cols = F,
         labels_col = rep("",14))
groups_tp53_all <- sapply(groups_tp53_all, function(g){ifelse(g == 1, "Low",
                                                              ifelse(g == 2,"High", "Intermediate"))})
hmdataMeta <- cbind(hmdataMeta,groups_tp53_all = as.factor(groups_tp53_all))

## Make 3 groups with low, intermediate and high ##
# Clustering based on BTG2 at all concentrations
d2 <- dist(hmdata_BTG2_all)
hc2 <- hclust(d2)
plot(hc2)
groups_BTG2_all <- cutree(hc2, k=3)
pheatmap(hmdata_BTG2_all, cluster_rows = T, cluster_cols = F,
         labels_col = rep("",14))
groups_BTG2_all <- sapply(groups_BTG2_all, function(g){ifelse(g == 1, "Intermediate",
                                                              ifelse(g == 2, "Low","High"))})
hmdataMeta <- cbind(hmdataMeta,groups_BTG2_all = as.factor(groups_BTG2_all))

## Make 3 groups with low, intermediate and high ##
# Clustering based on CDKN1A at all concentrations
d3 <- dist(hmdata_CDKN1A_all)
hc3 <- hclust(d3)
plot(hc3)
groups_CDKN1A_all <- cutree(hc3, k=3)
pheatmap(hmdata_CDKN1A_all, cluster_rows = T, cluster_cols = F,
         labels_col = rep("",14))
groups_CDKN1A_all <- sapply(groups_CDKN1A_all, function(g){ifelse(g == 1, "Low",
                                                                  ifelse(g == 2, "Intermediate","High"))})
hmdataMeta <- cbind(hmdataMeta,groups_CDKN1A_all = as.factor(groups_CDKN1A_all))

## Make 3 groups with low, intermediate and high ##
# Clustering based on MDM2 at all concentrations
d4 <- dist(hmdata_MDM2_all)
hc4 <- hclust(d4)
plot(hc4)
groups_MDM2_all <- cutree(hc4, k=3)
pheatmap(hmdata_MDM2_all, cluster_rows = T, cluster_cols = F,
         labels_col = rep("",14))
groups_MDM2_all <- sapply(groups_MDM2_all, function(g){ifelse(g == 1, "High",
                                                              ifelse(g == 2, "Intermediate","Low"))})
hmdataMeta <- cbind(hmdataMeta,groups_MDM2_all = as.factor(groups_MDM2_all))

##################################
########## FIGURE 1 A ############
##################################
doses <- c("0","0.1","1","3.3","10","33","100")
colAnot <- data.frame(Time = c(rep("8",length(doses)),rep("24",length(doses))),
                      Concentration = rep(doses,2),
                      row.names = colnames(hmdata_tp53_all))
rowAnot <- data.frame(TP53 = as.character(hmdataMeta$groups_tp53_all),
                      CDKN1A = as.character(hmdataMeta$groups_CDKN1A_all),
                      BTG2 = as.character(hmdataMeta$groups_BTG2_all),
                      MDM2 = as.character(hmdataMeta$groups_MDM2_all),
                      row.names = rownames(hmdata_tp53_all))

concclrs <- colorRampPalette(c("#FFFFFF",viridis(2)[1]))(7)# inferno(8)[5]))(7) #brewer.pal(7, "Reds")
names(concclrs) <- as.character(c(0,0.1,1,3.3,10,33,100))

clustclrs <- inferno(5) #colorRampPalette(c("#FFFFFF",viridis(8)[3]))(4)

timeclrs <- colorRampPalette(c("#FFFFFF",viridis(8)[3]))(4)

hmclrs <- colorRampPalette(c(viridis(5)[5],"#FFFFFF",viridis(5)[4]))(100)

ann_colors <- list(
  Concentration = concclrs,
  Time = c("8" = timeclrs[2], "24" = timeclrs[4]),
  TP53 = c("Low" = clustclrs[5], "Intermediate" = clustclrs[4],"High" = clustclrs[3]),
  CDKN1A = c("Low" = clustclrs[5], "Intermediate" = clustclrs[4],"High" = clustclrs[3]),
  BTG2 = c("Low" = clustclrs[5], "Intermediate" = clustclrs[4],"High" = clustclrs[3]),
  MDM2 = c("Low" = clustclrs[5], "Intermediate" = clustclrs[4],"High" = clustclrs[3]))

# Choose dims 5 x 7, portrait
# Fig1B_Heatmap_TP53
hm_tp53 <- pheatmap(hmdata_tp53_all, cluster_rows = T, cluster_cols = F,
                    color = hmclrs,
                    labels_col = rep("",14), 
                    annotation_row = rowAnot[,"TP53", drop=FALSE],
                    cutree_rows = 3,
                    annotation_col = colAnot,
                    annotation_colors = ann_colors,
                    annotation_names_row = F,
                    main = "TP53")
save_pheatmap_pdf(hm_tp53, paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig1B_Heatmap_TP53.pdf"))

##################################
########## FIGURE 1C-D ###########
##################################
# FigS2_Heatmap_MDM2
hm_mdm2 <- pheatmap(hmdata_MDM2_all, cluster_rows = T, cluster_cols = F,
                    color = hmclrs,
                    labels_col = rep("",14), 
                    cutree_rows = 3,
                    annotation_row = rowAnot[,"MDM2", drop=FALSE],
                    annotation_col = colAnot,
                    annotation_colors = ann_colors,
                    annotation_names_row = F,
                    main = "MDM2")
save_pheatmap_pdf(hm_mdm2, paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2A_Heatmap_MDM2.pdf"))
# FigS2_Heatmap_CDKN1A
hm_p21 <- pheatmap(hmdata_CDKN1A_all, cluster_rows = T, cluster_cols = F,
                   color = hmclrs,
                   labels_col = rep("",14), 
                   cutree_rows = 3,
                   annotation_row = rowAnot[,"CDKN1A", drop=FALSE],
                   annotation_col = colAnot,
                   annotation_colors = ann_colors,
                   annotation_names_row = F,
                   main = "CDKN1A")
save_pheatmap_pdf(hm_p21, paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2B_Heatmap_CDKN1A.pdf"))
# FigS2_Heatmap_BTG2
hm_btg2 <- pheatmap(hmdata_BTG2_all, cluster_rows = T, cluster_cols = F,
                    color = hmclrs,
                    labels_col = rep("",14), 
                    cutree_rows = 3,
                    annotation_row = rowAnot[,"BTG2", drop=FALSE],
                    annotation_col = colAnot,
                    annotation_colors = ann_colors,
                    annotation_names_row = F,
                    main = "BTG2")
save_pheatmap_pdf(hm_btg2, paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2C_Heatmap_BTG2.pdf"))

# FigS2_ClusterConservation
select <- hmdataMeta %>% dplyr::select(CELL_ID, groups_BTG2_all, groups_tp53_all,groups_MDM2_all,groups_CDKN1A_all) %>% 
  rename(groups_BTG2_all = "BTG2", groups_tp53_all = "TP53", groups_CDKN1A_all = "CDKN1A", groups_MDM2_all = "MDM2") %>%
  pivot_longer(cols =  c(TP53,BTG2,MDM2,CDKN1A), names_to = "gene", values_to = "group")

select$group <- factor(select$group, levels = c("Low","Intermediate","High"))
select$gene <- factor(select$gene,levels = c("TP53","MDM2","CDKN1A","BTG2"))
select <- select %>% mutate(CELL_ID = factor(CELL_ID, levels = names(sort(groups_tp53_all))))

ggplot(select) + geom_tile(aes(x = gene, y = CELL_ID, fill = group)) + 
  scale_fill_manual(name = "Group",
                    values = c("Low" = clustclrs[5], "Intermediate" = clustclrs[4],"High" = clustclrs[3])) + theme_classic() +
  xlab("Gene") + ylab("Sample ID")
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2D_ClusterConservation.pdf"),width = 5, height = 7)
dev.off()


# Make t-SNE plot
## Do PCA ##
pcaData_n <- as.data.frame(pcaData)
rownames(pcaData_n) <- pcaData_n$meanID
pcaData_n$meanID <- NULL
pcaDataFull <- left_join(pcaDataMeta, pcaData, by = "meanID")
res.pca <- prcomp(pcaData_n, scale = TRUE)

## Plot output ##

expl_var <- res.pca$sdev^2/sum(res.pca$sdev^2)
barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",
        names.arg=paste0("PC",seq(1:50)), col="darkgreen")

autoplot(res.pca, data = pcaDataFull, colour = 'DESIGN') + theme_classic() + scale_color_manual(values=colorsCodes)
autoplot(res.pca, data = pcaDataFull, colour = 'CELL_ID') + theme_classic() + scale_color_viridis(discrete = T) #+ scale_color_manual(values=colorsCodes)
#ggsave(filename = paste0(outputDir_GR,DATE,"_logNorm/PCA_medium_AllGenes.pdf"),width = 8, height = 5.5)

N_perm <- 10
expl_var_perm <- matrix(NA, ncol = length(res.pca$sdev), nrow = N_perm)
for (k in 1:N_perm) {
  expr_perm <- apply(t(pcaData_n),2,sample)
  PC_perm <- prcomp(t(expr_perm), center=TRUE, scale=FALSE)
  expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
}
plot(expl_var[1:50]~seq(1:50), ylab="EXPLAINED VARIANCE",
     col="darkgreen", type='o', xlab="PRINCIPAL COMPONENTS")
lines(colMeans(expl_var_perm)[1:50]~seq(1:50),col="red")
legend("topright", c("Explained by PCS", "Explained by chance"),
       fill=c("darkgreen","red"), inset=0.02)

pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
plot(pval[1:50]~seq(1:50),col="darkred",type='o',
     xlab="PRINCIPAL COMPONENTS",ylab="PVALUE")
optPC<-head(which(pval>=0.05),1)-1
print(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC))

## Executing the algorithm on curated data ##
# lowest = 0.578858
tsne <- Rtsne(pcaDataFull %>% dplyr::select(-c("EXP_ID","CELL_ID","CONCENTRATION","TIMEPOINT","DESIGN","SAMPLE_CODE","COLOR","meanID")), 
              dims = 2, init_dims = 29, perplexity=22, verbose=TRUE, max_iter = 6000,check_duplicates = FALSE)

## Plotting ##
tsne_plot <- as_tibble(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], SampleID = pcaDataFull$CELL_ID))
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=SampleID), size = 2) + 
  xlab("t-SNE variable 1") + ylab("t-SNE variable 2") +
  theme_classic() +
  scale_color_viridis(discrete = T, name = "Sample ID")
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS1D_tSNE_plot_AllGenens.pdf"), width = 8, height = 5)


##################################
######### FIGURE 1 F-G ###########
##################################
# Make plot Fig 1F; group specific basal TP53 expression
hmdata_tp53_all$CELL_ID <- rownames(hmdata_tp53_all)
groupTP53 <- left_join(hmdata_tp53_all,hmdataMeta, by = "CELL_ID")
groupSum <-  groupTP53 %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_tp53_all) %>% summarise(meanTP53_tp8 = mean(c0_tp8, na.rm = T),
                                                                                                    meanTP53_tp24 = mean(c0_tp24, na.rm = T))
groupInfo <- groupTP53 %>% filter(!CELL_ID=="HepG2_WT") %>% mutate(groups_tp53_all = factor(groups_tp53_all, levels = c("Low","Intermediate","High")))
groupSum <- groupSum %>% mutate(groups_tp53_all = factor(groups_tp53_all, levels = c("Low","Intermediate","High")))
groupHepG2 <- groupTP53[groupTP53$EXP_ID == "HepG2",]

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_tp53_all, y = c0_tp8, 
                                    fill = groups_tp53_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_tp53_all, y = c0_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_tp53_all, y = meanTP53_tp8, color = groups_tp53_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_tp53_all, y = c0_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_tp53_all, y = c0_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("TP53 expression cluster") + ylab("Basal TP53 expression") +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig1C_tp8_basalTP53expression.pdf"),width = 3, height = 2.5)

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_tp53_all, y = c0_tp24, 
                                    fill = groups_tp53_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_tp53_all, y = c0_tp24),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_tp53_all, y = meanTP53_tp24, color = groups_tp53_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_tp53_all, y = c0_tp24), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_tp53_all, y = c0_tp24, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("TP53 expression cluster") + ylab("Basal TP53 expression") +
  ggtitle("24 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_tp24_basalTP53expression.pdf"),width = 3, height = 2.5)

groupSum <-  groupTP53 %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_tp53_all) %>% summarise(meanTP53_tp8 = mean(c3.3_tp8, na.rm = T),
                                                                                                    meanTP53_tp24 = mean(c3.3_tp24, na.rm = T))
groupSum <- groupSum %>% mutate(groups_tp53_all = factor(groups_tp53_all, levels = c("Low","Intermediate","High")))
ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_tp53_all, y = c3.3_tp8, 
                                    fill = groups_tp53_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_tp53_all, y = c3.3_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_tp53_all, y = meanTP53_tp8, color = groups_tp53_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_tp53_all, y = c3.3_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_tp53_all, y = c3.3_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("TP53 expression cluster") + 
  ylab(bquote(atop("TP53 expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2E_tp8_TP53expressionAt3.3.pdf"),width = 3, height = 2.5)

# Make plot Fig 1G; group specific basal TP53 expression
hmdata_MDM2_all$CELL_ID <- rownames(hmdata_MDM2_all)
groupMDM2 <- left_join(hmdata_MDM2_all,hmdataMeta, by = "CELL_ID")
groupSum <-  groupMDM2 %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_MDM2_all) %>% summarise(meanMDM2_tp8 = mean(c0_tp8, na.rm = T),
                                                                                                    meanMDM2_tp24 = mean(c0_tp24, na.rm = T))
groupInfo <- groupMDM2 %>% filter(!CELL_ID=="HepG2_WT") %>% mutate(groups_MDM2_all = factor(groups_MDM2_all, levels = c("Low","Intermediate","High")))
groupSum <- groupSum %>% mutate(groups_MDM2_all = factor(groups_MDM2_all, levels = c("Low","Intermediate","High")))
groupHepG2 <- groupMDM2[groupMDM2$EXP_ID == "HepG2",]

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_MDM2_all, y = c0_tp8, 
                                    fill = groups_MDM2_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_MDM2_all, y = c0_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_MDM2_all, y = meanMDM2_tp8, color = groups_MDM2_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_MDM2_all, y = c0_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_MDM2_all, y = c0_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("MDM2 expression cluster") + ylab("Basal MDM2 expression") +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig1D_tp8_basalMDM2expression.pdf"),width = 3, height = 2.5)

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_MDM2_all, y = c0_tp24, 
                                    fill = groups_MDM2_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_MDM2_all, y = c0_tp24),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_MDM2_all, y = meanMDM2_tp24, color = groups_MDM2_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_MDM2_all, y = c0_tp24), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_MDM2_all, y = c0_tp24, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("MDM2 expression cluster") + ylab("Basal MDM2 expression") +
  ggtitle("24 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_tp24_basalMDM2expression.pdf"),width = 3, height = 2.5)

groupSum <-  groupMDM2 %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_MDM2_all) %>% summarise(meanMDM2_tp8 = mean(c3.3_tp8, na.rm = T),
                                                                                                    meanMDM2_tp24 = mean(c3.3_tp24, na.rm = T))
groupSum <- groupSum %>% mutate(groups_MDM2_all = factor(groups_MDM2_all, levels = c("Low","Intermediate","High")))
ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_MDM2_all, y = c3.3_tp8, 
                                    fill = groups_MDM2_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_MDM2_all, y = c3.3_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_MDM2_all, y = meanMDM2_tp8, color = groups_MDM2_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_MDM2_all, y = c3.3_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_MDM2_all, y = c3.3_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("MDM2 expression cluster") + 
  ylab(bquote(atop("MDM2 expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2F_tp8_MDM2expressionAt3.3.pdf"),width = 3, height = 2.5)

# Make plot Fig 1H; group specific basal TP53 expression
hmdata_CDKN1A_all$CELL_ID <- rownames(hmdata_CDKN1A_all)
groupCDKN1A <- left_join(hmdata_CDKN1A_all,hmdataMeta, by = "CELL_ID")
groupSum <-  groupCDKN1A %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_CDKN1A_all) %>% summarise(meanCDKN1A_tp8 = mean(c0_tp8, na.rm = T),
                                                                                                        meanCDKN1A_tp24 = mean(c0_tp24, na.rm = T))
groupInfo <- groupCDKN1A %>% filter(!CELL_ID=="HepG2_WT") %>% mutate(groups_CDKN1A_all = factor(groups_CDKN1A_all, levels = c("Low","Intermediate","High")))
groupSum <- groupSum %>% mutate(groups_CDKN1A_all = factor(groups_CDKN1A_all, levels = c("Low","Intermediate","High")))
groupHepG2 <- groupCDKN1A[groupCDKN1A$EXP_ID == "HepG2",]

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_CDKN1A_all, y = c0_tp8, 
                                    fill = groups_CDKN1A_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_CDKN1A_all, y = c0_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_CDKN1A_all, y = meanCDKN1A_tp8, color = groups_CDKN1A_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_CDKN1A_all, y = c0_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_CDKN1A_all, y = c0_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("CDKN1A expression cluster") + ylab("Basal CDKN1A expression") +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig1E_tp8_basalCDKN1Aexpression.pdf"),width = 3, height = 2.5)

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_CDKN1A_all, y = c0_tp24, 
                                    fill = groups_CDKN1A_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_CDKN1A_all, y = c0_tp24),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_CDKN1A_all, y = meanCDKN1A_tp24, color = groups_CDKN1A_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_CDKN1A_all, y = c0_tp24), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_CDKN1A_all, y = c0_tp24, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("CDKN1A expression cluster") + ylab("Basal CDKN1A expression") +
  ggtitle("24 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_tp24_basalCDKN1Aexpression.pdf"),width = 3, height = 2.5)

groupSum <-  groupCDKN1A %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_CDKN1A_all) %>% summarise(meanCDKN1A_tp8 = mean(c3.3_tp8, na.rm = T),
                                                                                                        meanCDKN1A_tp24 = mean(c3.3_tp24, na.rm = T))
groupSum <- groupSum %>% mutate(groups_CDKN1A_all = factor(groups_CDKN1A_all, levels = c("Low","Intermediate","High")))
ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_CDKN1A_all, y = c3.3_tp8, 
                                    fill = groups_CDKN1A_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_CDKN1A_all, y = c3.3_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_CDKN1A_all, y = meanCDKN1A_tp8, color = groups_CDKN1A_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_CDKN1A_all, y = c3.3_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_CDKN1A_all, y = c3.3_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("CDKN1A expression cluster") + 
  ylab(bquote(atop("CDKN1A expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2G_tp8_CDKN1AexpressionAt3.3.pdf"),width = 3, height = 2.5)

# Make plot Fig 1I; group specific basal BTG2 expression
hmdata_BTG2_all$CELL_ID <- rownames(hmdata_BTG2_all)
groupBTG2 <- left_join(hmdata_BTG2_all,hmdataMeta, by = "CELL_ID")
groupSum <-  groupBTG2 %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_BTG2_all) %>% summarise(meanBTG2_tp8 = mean(c0_tp8, na.rm = T),
                                                                                                    meanBTG2_tp24 = mean(c0_tp24, na.rm = T))
groupInfo <- groupBTG2 %>% filter(!CELL_ID=="HepG2_WT") %>% mutate(groups_BTG2_all = factor(groups_BTG2_all, levels = c("Low","Intermediate","High")))
groupSum <- groupSum %>% mutate(groups_BTG2_all = factor(groups_BTG2_all, levels = c("Low","Intermediate","High")))
groupHepG2 <- groupBTG2[groupBTG2$EXP_ID == "HepG2",]

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_BTG2_all, y = c0_tp8, 
                                    fill = groups_BTG2_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_BTG2_all, y = c0_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_BTG2_all, y = meanBTG2_tp8, color = groups_BTG2_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_BTG2_all, y = c0_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_BTG2_all, y = c0_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("BTG2 expression cluster") + ylab("Basal BTG2 expression") +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig1F_tp8_basalBTG2expression.pdf"),width = 3, height = 2.5)

ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_BTG2_all, y = c0_tp24, 
                                    fill = groups_BTG2_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_BTG2_all, y = c0_tp24),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_BTG2_all, y = meanBTG2_tp24, color = groups_BTG2_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_BTG2_all, y = c0_tp24), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_BTG2_all, y = c0_tp24, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("BTG2 expression cluster") + ylab("Basal BTG2 expression") +
  ggtitle("24 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_tp24_basalBTG2expression.pdf"),width = 3, height = 2.5)

groupSum <-  groupBTG2 %>% filter(!CELL_ID=="HepG2_WT") %>% group_by(groups_BTG2_all) %>% summarise(meanBTG2_tp8 = mean(c3.3_tp8, na.rm = T),
                                                                                                    meanBTG2_tp24 = mean(c3.3_tp24, na.rm = T))
groupSum <- groupSum %>% mutate(groups_BTG2_all = factor(groups_BTG2_all, levels = c("Low","Intermediate","High")))
ggplot() + 
  geom_violin(data = groupInfo, aes(x = groups_BTG2_all, y = c3.3_tp8, 
                                    fill = groups_BTG2_all), alpha = 0.25) +
  scale_fill_manual(values = rep('lightgrey',3)) +
  geom_jitter(data = groupInfo, aes(x = groups_BTG2_all, y = c3.3_tp8),
              width = 0.1, color = "grey50") +
  geom_point(data = groupSum, aes(x = groups_BTG2_all, y = meanBTG2_tp8, color = groups_BTG2_all),
             size = 4) +
  geom_point(data = groupHepG2, aes(x = groups_BTG2_all, y = c3.3_tp8), color = inferno(3)[2], size = 0.75) +
  geom_text(data = groupHepG2, aes(x = groups_BTG2_all, y = c3.3_tp8, label = EXP_ID), 
            size=3 , color = inferno(3)[2], hjust = -0.15) +
  scale_color_viridis(discrete = T) +
  theme_classic() + xlab("BTG2 expression cluster") + 
  ylab(bquote(atop("BTG2 expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  ggtitle("8 hr") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS2H_tp8_BTG2expressionAt3.3.pdf"),width = 3, height = 2.5)

# Add the expression group to the df
hmdataSum <- left_join(hmdataSum,hmdataMeta %>% dplyr::select("CELL_ID","groups_tp53_all"), by = "CELL_ID")

# Select PHH data
bpdata_phh <- hmdataSum %>% filter(EXP_ID == "PHH")
bpdata_hg2 <- hmdata %>% filter(EXP_ID == "HepG2")

# 8 hr & 24 hr data
bpdata_8hr <- bpdata_phh %>% filter(TIMEPOINT == "8")
bpdata_24hr <- bpdata_phh %>% filter(TIMEPOINT == "24")

##################################
########## FIGURE 2A-D ###########
##################################

# HepG2, 8hr, downstream target & TP53 expression
bpdata_hg2Select <- bpdata_hg2 %>% dplyr::select(-c(MDM2_23384,MDM2_27496,MDM2_27497,MDM2_sum)) 
bpdata_hg2Pivot <- pivot_longer(bpdata_hg2Select, 
                                cols = c(TP53, MDM2, CDKN1A, BTG2)) 
bpdata_hg2Mean <- bpdata_hg2Pivot %>% 
  group_by(DESIGN, SAMPLE_CODE, EXP_ID, CELL_ID, TIMEPOINT, CONCENTRATION, CONC_TIME, name) %>% 
  summarise_all(list(mean = mean, sd = sd)) %>% ungroup 
bpdata_hg2Mean <- bpdata_hg2Mean %>% mutate(relevel(CONCENTRATION, ref = "0"))

hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==8 &
                                     bpdata_hg2Mean$name == "TP53",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==8 &
                                      bpdata_hg2Pivot$name == "TP53",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, 
                                 y=mean, yend=mean, color = CONCENTRATION),
               size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("TP53 normalised count") +
  #xlim((-10),(6)) +
  ylim(2,8)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2A_HepG2_8hr_TP53Expr_MeanSD.pdf"),width = 3, height = 2.5)

hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==8 &
                                     bpdata_hg2Mean$name == "MDM2",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==8 &
                                      bpdata_hg2Pivot$name == "MDM2",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, y=mean, yend=mean, color = CONCENTRATION), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("MDM2 normalised count") +
  ylim(7,12)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2B_HepG2_8hr_MDM2Expr_MeanSD.pdf"),width = 3, height = 2.5)

hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==8 &
                                     bpdata_hg2Mean$name == "CDKN1A",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==8 &
                                      bpdata_hg2Pivot$name == "CDKN1A",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, y=mean, yend=mean, color = CONCENTRATION), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("CDKN1A normalised count") +
  ylim(1,7)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2C_HepG2_8hr_CDKN1AExpr_MeanSD.pdf"),width = 3, height = 2.5)

# BTG2, HepG2, 8 hr
hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==8 &
                                     bpdata_hg2Mean$name == "BTG2",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==8 &
                                      bpdata_hg2Pivot$name == "BTG2",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, y=mean, yend=mean, color = CONCENTRATION), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("BTG2 normalised count") +
  ylim(4,10)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2D_HepG2_8hr_BTG2Expr_MeanSD.pdf"),width = 3, height = 2.5)

# HepG2, 24hr, downstream target & TP53 expression
hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==24 &
                                     bpdata_hg2Mean$name == "TP53",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==24 &
                                      bpdata_hg2Pivot$name == "TP53",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, y=mean, yend=mean, color = CONCENTRATION), 
               size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("TP53 normalised count") +
  ylim(2,8)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3A_HepG2_24hr_TP53Expr_MeanSD.pdf"),width = 3, height = 2.5)

hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==24 &
                                     bpdata_hg2Mean$name == "MDM2",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==24 &
                                      bpdata_hg2Pivot$name == "MDM2",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, y=mean, yend=mean, color = CONCENTRATION), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("MDM2 normalised count") +
  ylim(7,12)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3B_HepG2_24hr_MDM2Expr_MeanSD.pdf"),width = 3, height = 2.5)

hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==24 &
                                     bpdata_hg2Mean$name == "CDKN1A",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==24 &
                                      bpdata_hg2Pivot$name == "CDKN1A",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, y=mean, yend=mean, color = CONCENTRATION), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("CDKN1A normalised count") +
  ylim(1,7)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3C_HepG2_24hr_CDKN1AExpr_MeanSD.pdf"),width = 3, height = 2.5)

hline <- data.frame(bpdata_hg2Mean[bpdata_hg2Mean$TIMEPOINT==24 &
                                     bpdata_hg2Mean$name == "BTG2",c("mean","CONCENTRATION")]) 
hline <- hline[order(hline$CONCENTRATION),]
hline$axisLab <- seq(1,7)
ggplot() + 
  geom_point(data = bpdata_hg2Pivot[bpdata_hg2Pivot$TIMEPOINT==24 &
                                      bpdata_hg2Pivot$name == "BTG2",],
             aes(x = CONCENTRATION, y = value, color = CONCENTRATION),size = 2) +
  geom_segment(data = hline, aes(x=axisLab-0.4, xend=axisLab+0.4, y=mean, yend=mean, color = CONCENTRATION), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = "") + 
  ggtitle("HepG2, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("BTG2 normalised count") +
  ylim(4,10)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3D_HepG2_24hr_BTG2Expr_MeanSD.pdf"),width = 3, height = 2.5)

##################################
########## FIGURE 2A-D ############
##################################
# Make a boxplot
# 8 hr
ggplot(bpdata_8hr, 
       aes(x = CONCENTRATION, y = TP53, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("TP53 normalised count") +
  ggtitle("PHHs, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  ylim(2,8)
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2A_PHH_8hr_TP53Expr_Boxplot.pdf"),width = 3, height = 2.5)

# 8 hr, MDM2
ggplot(bpdata_8hr, 
       aes(x = CONCENTRATION, y = MDM2, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("MDM2 normalised count") +
  ylim(7,12) + 
  ggtitle("PHHs, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2B_PHH_8hr_MDM2Expr_Boxplot.pdf"),width = 3, height = 2.5)

# 8 hr, CDKN1A
ggplot(bpdata_8hr, 
       aes(x = CONCENTRATION, y = CDKN1A, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("CDKN1A normalised count") +
  ylim(1,7) + 
  ggtitle("PHHs, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2C_PHH_8hr_CDKN1AExpr_Boxplot.pdf"),width = 3, height = 2.5)

# 8 hr, BTG2
ggplot(bpdata_8hr, 
       aes(x = CONCENTRATION, y = BTG2, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("BTG2 normalised count") +
  ylim(4,10) + 
  ggtitle("PHHs, 8 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2D_PHH_8hr_BTG2Expr_Boxplot.pdf"),width = 3, height = 2.5)

# 24 hr
ggplot(bpdata_24hr, 
       aes(x = CONCENTRATION, y = TP53, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("TP53 normalised count") +
  ylim(2,8) + 
  ggtitle("PHHs, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3A_PHH_24hr_TP53Expr_Boxplot.pdf"),width = 3, height = 2.5)

# 24 hr, MDM2
ggplot(bpdata_24hr, 
       aes(x = CONCENTRATION, y = MDM2, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("MDM2 normalised count") +
  ylim(7,12) + 
  ggtitle("PHHs, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3B_PHH_24hr_MDM2Expr_Boxplot.pdf"),width = 3, height = 2.5)

# 24 hr, CDKN1A
ggplot(bpdata_24hr, 
       aes(x = CONCENTRATION, y = CDKN1A, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("CDKN1A normalised count") +
  ylim(1,7) + 
  ggtitle("PHHs, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3C_PHH_24hr_CDKN1AExpr_Boxplot.pdf"),width = 3, height = 2.5)

# 24 hr, BTG2
ggplot(bpdata_24hr, 
       aes(x = CONCENTRATION, y = BTG2, fill = CONCENTRATION)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_viridis(discrete = T, name = "") + 
  theme(axis.text.x = element_text(angle = 0)) +
  xlab(expression("Concentration ("*mu*"M)")) + 
  ylab("BTG2 normalised count") +
  ylim(4,10) + 
  ggtitle("PHHs, 24 hr") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS3D_PHH_24hr_BTG2Expr_Boxplot.pdf"),width = 3, height = 2.5)

# ######### ######### ######### ######### #########
# ######### ###### Do for 8 hr ######### #########
# ######### ######### ######### ######### #########


######### ######### ######### ######### ######### ######### ######### ######### 
######### XII. Do correlation analysis  ######### ######### ######### #########
######### ######### ######### ######### ######### ######### ######### ######### 

## Make correlation plot for medium 24hr vs medium 8 hr for all genes
# 8
df_med_8 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 0, 
                                                TIMEPOINT == 8))
rownames(df_med_8) <- df_med_8$CELL_ID

# 24
df_med_24 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 0,
                                                 TIMEPOINT == 24))
rownames(df_med_24) <- df_med_24$CELL_ID

##################################
########## FIGURE 3 C ############
##################################

colors <- c("BTG2" = "#440154FF", "MDM2" = "#21908CFF", "CDKN1A" = "#FDE725FF")
ggplot(data = df_med_8) + 
  theme_classic() + 
  geom_smooth(aes(x = TP53,y = BTG2),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = BTG2, color = "BTG2")) + 
  geom_smooth(aes(x = TP53,y = MDM2),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2, color = "MDM2")) + 
  geom_smooth(aes(x = TP53,y = CDKN1A),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = CDKN1A, color = "CDKN1A")) + 
  labs(x = "Basal TP53 expression",
       y = "Basal downstream gene\nexpression",
       color = "Gene") +
  scale_color_manual(values = colors) + 
  annotate(geom="text", x=6.75, y=9.5, color = "#21908CFF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_8[,"TP53"], y = df_med_8[,"MDM2"])$estimate,2))) + 
  annotate(geom="text", x=6.75, y=7, color = "#440154FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_8[,"TP53"], y = df_med_8[,"BTG2"])$estimate,2))) + 
  annotate(geom="text", x=6.75, y=3, color = "#FDE725FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_8[,"TP53"], y = df_med_8[,"CDKN1A"])$estimate,2))) +
  ggtitle("8 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2E_DownstreamTargetCorrelation_8hrBasalTP53",".pdf"),width = 4, height = 2.75)

colors <- c("BTG2" = "#440154FF", "MDM2" = "#21908CFF", "CDKN1A" = "#FDE725FF")
ggplot(data = df_med_24) + 
  theme_classic() + 
  geom_smooth(aes(x = TP53,y = BTG2),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = BTG2, color = "BTG2")) + 
  geom_smooth(aes(x = TP53,y = MDM2),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2, color = "MDM2")) + 
  geom_smooth(aes(x = TP53,y = CDKN1A),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = CDKN1A, color = "CDKN1A")) + 
  labs(x = "Basal TP53 expression",
       y = "Basal downstream gene\nexpression",
       color = "Gene") +
  scale_color_manual(values = colors) + 
  annotate(geom="text", x=6.75, y=8, color = "#21908CFF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_24[,"TP53"], y = df_med_24[,"MDM2"])$estimate,2))) + 
  annotate(geom="text", x=6.75, y=7, color = "#440154FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_24[,"TP53"], y = df_med_24[,"BTG2"])$estimate,2))) + 
  annotate(geom="text", x=6.75, y=2.5, color = "#FDE725FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_24[,"TP53"], y = df_med_24[,"CDKN1A"])$estimate,2))) +
  ggtitle("24 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_DownstreamTargetCorrelation_24hrBasalTP53",".pdf"),width = 4, height = 2.75)

colors <- c("MDM2_23384" = "#440154FF", "MDM2_27496" = "#21908CFF", "MDM2_27497" = "#FDE725FF")
ggplot(data = df_med_8) + 
  theme_classic() + 
  geom_smooth(aes(x = TP53,y = MDM2_23384),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2_23384, color = "MDM2_23384")) + 
  geom_smooth(aes(x = TP53,y = MDM2_27496),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2_27496, color = "MDM2_27496")) + 
  geom_smooth(aes(x = TP53,y = MDM2_27497),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2_27497, color = "MDM2_27497")) + 
  scale_color_manual(values = colors) + 
  annotate(geom="text", x=4.75, y=11.5, color = "#440154FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_8[,"TP53"], y = df_med_8[,"MDM2_23384"])$estimate,2))) + 
  annotate(geom="text", x=4.75, y=8, color = "#21908CFF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_8[,"TP53"], y = df_med_8[,"MDM2_27496"])$estimate,2))) + 
  annotate(geom="text", x=4.75, y=5.75, color = "#FDE725FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_8[,"TP53"], y = df_med_8[,"MDM2_27497"])$estimate,2))) +
  labs(x = "Basal TP53 expression",
       y = "Basal MDM2 expression",
       color = "MDM2 probe") +
  ggtitle("8 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_MDM2TargetCorrelation_8hrBasalTP53",".pdf"),width = 4, height = 2.75)

colors <- c("MDM2_23384" = "#440154FF", "MDM2_27496" = "#21908CFF", "MDM2_27497" = "#FDE725FF")
ggplot(data = df_med_24) + 
  theme_classic() + 
  geom_smooth(aes(x = TP53,y = MDM2_23384),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2_23384, color = "MDM2_23384")) + 
  geom_smooth(aes(x = TP53,y = MDM2_27496),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2_27496, color = "MDM2_27496")) + 
  geom_smooth(aes(x = TP53,y = MDM2_27497),method='lm', formula= y~x, se = F, colour = "grey") + 
  geom_point(aes(x = TP53,y = MDM2_27497, color = "MDM2_27497")) + 
  scale_color_manual(values = colors) + 
  annotate(geom="text", x=5.25, y=11.5, color = "#440154FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_24[,"TP53"], y = df_med_24[,"MDM2_23384"])$estimate,2))) + 
  annotate(geom="text", x=5.25, y=8.75, color = "#21908CFF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_24[,"TP53"], y = df_med_24[,"MDM2_27496"])$estimate,2))) + 
  annotate(geom="text", x=5.25, y=6.5, color = "#FDE725FF",parse = T,
           label= sprintf("rho == %0.2f", round(cor.test(x = df_med_24[,"TP53"], y = df_med_24[,"MDM2_27497"])$estimate,2))) +
  labs(x = "Basal TP53 expression",
       y = "Basal MDM2 expression",
       color = "MDM2 probe") +
  ggtitle("24 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_MDM2TargetCorrelation_24hrBasalTP53",".pdf"),width = 4, height = 2.75)

## Get the data frames

df_cpt0.1_8 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 0.1, 
                                                   TIMEPOINT == 8))
df_cpt0.1_24 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 0.1, 
                                                    TIMEPOINT == 24))
df_cpt1_8 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 1, 
                                                 TIMEPOINT == 8))
df_cpt1_24 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 1, 
                                                  TIMEPOINT == 24))

# # Cisplatin 3.3 uM
df_cpt3_8 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 3.3, 
                                                 TIMEPOINT == 8))
df_cpt3_24 <- as.data.frame(bpdata_phh %>% filter(CONCENTRATION == 3.3, 
                                                  TIMEPOINT == 24))


##################################
########## FIGURE 3 A ############
##################################

a <- df_cpt0.1_8 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
b <- df_cpt1_8 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
c <- df_cpt3_8 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
d <- df_med_8 %>% dplyr::select(CELL_ID, TP53) %>% rename(TP53 = "TP53_0")
df_tp53_8 <- rbind(merge(a,d,by = "CELL_ID"),
                   merge(b,d,by = "CELL_ID"),
                   merge(c,d,by = "CELL_ID"))

cor.test(x = merge(b,d,by = "CELL_ID")$TP53_0, y = merge(b,d,by = "CELL_ID")$TP53)
ggplot(data = df_tp53_8) + 
  theme_classic() + 
  geom_abline(slope = 1, lty = "dashed", color = "grey") +
  geom_point(aes(x = TP53_0,y = TP53, group = CONCENTRATION, color = CONCENTRATION)) + 
  geom_smooth(aes(x = TP53_0,y = TP53, group = CONCENTRATION, color = CONCENTRATION),
              method='lm', formula= y~x, se = F) + 
  xlab("Basal TP53 expression") + 
  ylab("TP53 expression\nafter cisplatin exposure") +
  scale_color_viridis(discrete = T, name = expression("Conc. ("*mu*"M)")) +
  annotate(geom="text", x=5, y=7.5, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_tp53_8 %>% filter(CONCENTRATION == 0.1) %>% pull(TP53_0), 
                                                         y = df_tp53_8 %>% filter(CONCENTRATION == 0.1) %>% pull(TP53))$estimate,3))) + 
  annotate(geom="text", x=5, y=7.2, color = "#21908CFF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_tp53_8 %>% filter(CONCENTRATION == 1) %>% pull(TP53_0), 
                                                         y = df_tp53_8 %>% filter(CONCENTRATION == 1) %>% pull(TP53))$estimate,3))) + 
  annotate(geom="text", x=5, y=6.9, color = "#FDE725FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_tp53_8 %>% filter(CONCENTRATION == 3.3) %>% pull(TP53_0), 
                                                         y = df_tp53_8 %>% filter(CONCENTRATION == 3.3) %>% pull(TP53))$estimate,3))) + 
  xlim(4,8) + ylim(4,8) +
  ggtitle("8 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS5A_ConcentrationCorrelation_8hrBasalTP53",".pdf"),width = 4, height = 2.75)

e <- df_cpt0.1_24 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
f <- df_cpt1_24 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
g <- df_cpt3_24 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
h <- df_med_24 %>% dplyr::select(CELL_ID, TP53) %>% rename(TP53 = "TP53_0")
df_tp53_24 <- rbind(merge(e,h,by = "CELL_ID"),
                    merge(f,h,by = "CELL_ID"),
                    merge(g,h,by = "CELL_ID"))

cor.test(x = merge(b,d,by = "CELL_ID")$TP53_0, y = merge(b,d,by = "CELL_ID")$TP53)
ggplot(data = df_tp53_24) + 
  theme_classic() + 
  geom_abline(slope = 1, lty = "dashed", color = "grey") +
  geom_point(aes(x = TP53_0,y = TP53, group = CONCENTRATION, color = CONCENTRATION)) + 
  geom_smooth(aes(x = TP53_0,y = TP53, group = CONCENTRATION, color = CONCENTRATION),
              method='lm', formula= y~x, se = F) + 
  xlab("Basal TP53 expression") + 
  ylab("TP53 expression\nafter cisplatin exposure") +
  scale_color_viridis(discrete = T, name = expression("Conc. ("*mu*"M)")) +
  annotate(geom="text", x=5, y=7.5, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_tp53_24 %>% filter(CONCENTRATION == 0.1) %>% pull(TP53_0), 
                                                         y = df_tp53_24 %>% filter(CONCENTRATION == 0.1) %>% pull(TP53))$estimate,3))) + 
  annotate(geom="text", x=5, y=7.2, color = "#21908CFF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_tp53_24 %>% filter(CONCENTRATION == 1) %>% pull(TP53_0), 
                                                         y = df_tp53_24 %>% filter(CONCENTRATION == 1) %>% pull(TP53))$estimate,3))) + 
  annotate(geom="text", x=5, y=6.9, color = "#FDE725FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_tp53_24 %>% filter(CONCENTRATION == 3.3) %>% pull(TP53_0), 
                                                         y = df_tp53_24 %>% filter(CONCENTRATION == 3.3) %>% pull(TP53))$estimate,3))) + 
  xlim(4,8) + ylim(4,8) +
  ggtitle("24 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS5A_ConcentrationCorrelation_24hrBasalTP53",".pdf"),width = 4, height = 2.75)

##################################
########## FIGURE 3 B ############
##################################
d <- df_med_8 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
i <- df_cpt0.1_24 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
j <- df_cpt1_24 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53)
k <- df_cpt3_24 %>% dplyr::select(CELL_ID,CONCENTRATION, TP53) 
l <- df_med_24 %>% dplyr::select(CELL_ID, CONCENTRATION,TP53) 
df_8_24 <- rbind(merge(a,i,by = c("CELL_ID","CONCENTRATION"), suffix = c("_8","_24")),
                 merge(b,j,by = c("CELL_ID","CONCENTRATION"), suffix = c("_8","_24")),
                 merge(c,k,by = c("CELL_ID","CONCENTRATION"), suffix = c("_8","_24")),
                 merge(d,l,by = c("CELL_ID","CONCENTRATION"), suffix = c("_8","_24")))

ggplot() + 
  theme_classic() + 
  geom_abline(slope = 1, lty = "dashed", color = "grey") +
  geom_point(data = df_8_24,aes(x = TP53_8,y = TP53_24, group = CONCENTRATION, color = CONCENTRATION)) + 
  geom_smooth(data = df_8_24,aes(x = TP53_8,y = TP53_24, group = CONCENTRATION, color = CONCENTRATION),
              method='lm', formula= y~x, se = F) + 
  xlab("TP53 expression at 8 hr") + 
  ylab("TP53 expression at 24 hr") +
  scale_color_viridis(discrete = T, name = expression("Conc. ("*mu*"M)")) +
  annotate(geom="text", x=7, y=5.6, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_8_24 %>% filter(CONCENTRATION == 0) %>% pull(TP53_8), 
                                                         y = df_8_24 %>% filter(CONCENTRATION == 0) %>% pull(TP53_24))$estimate,3))) +
  annotate(geom="text", x=7, y=5.2, color = "#31688EFF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_8_24 %>% filter(CONCENTRATION == 0.1) %>% pull(TP53_8), 
                                                         y = df_8_24 %>% filter(CONCENTRATION == 0.1) %>% pull(TP53_24))$estimate,3))) +
  annotate(geom="text", x=7, y=4.8, color = "#35B779FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_8_24 %>% filter(CONCENTRATION == 1) %>% pull(TP53_8), 
                                                         y = df_8_24 %>% filter(CONCENTRATION == 1) %>% pull(TP53_24))$estimate,3))) +
  annotate(geom="text", x=7, y=4.4, color = "#FDE725FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_8_24 %>% filter(CONCENTRATION == 3.3) %>% pull(TP53_8), 
                                                         y = df_8_24 %>% filter(CONCENTRATION == 3.3) %>% pull(TP53_24))$estimate,3))) +
  ggtitle("") +
  xlim(4,8) + ylim(4,8)
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig5B_TimepointCorrelation",".pdf"),width = 4, height = 2.75)


##################################
########## FIGURE 3 D ############
##################################
cormat <- round(cor(x = df_med_8[,PROTS_OF_INTEREST], y = df_med_8[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,PROTS_OF_INTEREST], y = df_med_8[,PROTS_OF_INTEREST]),3)
row1 <- cormat[1,]
rowA <- cormat.pval[1,]
cormat <- round(cor(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt0.1_8[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt0.1_8[,PROTS_OF_INTEREST]),3)
row2 <- cormat[1,]
rowB <- cormat.pval[1,]
cormat <- round(cor(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt1_8[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt1_8[,PROTS_OF_INTEREST]),3)
row3 <- cormat[1,]
rowC <- cormat.pval[1,]
cormat <- round(cor(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt3_8[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt3_8[,PROTS_OF_INTEREST]),3)
row4 <- cormat[1,]
rowD <- cormat.pval[1,]
CorMat1 <- matrix(c(row1,row2,row3,row4), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),PROTS_OF_INTEREST))
CorMat2 <- matrix(c(rowA,rowB,rowC,rowD), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),PROTS_OF_INTEREST))
plotCor(CorMat1,CorMat2, title = "8 hr", textSize = 4, full = T)
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2F_CorrelationMatrix_8hr",".pdf"),width = 4.25, height = 4.25)

cormat <- round(cor(x = df_med_24[,PROTS_OF_INTEREST], y = df_med_24[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,PROTS_OF_INTEREST], y = df_med_24[,PROTS_OF_INTEREST]),3)
row1 <- cormat[1,]
rowA <- cormat.pval[1,]
cormat <- round(cor(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt0.1_24[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt0.1_24[,PROTS_OF_INTEREST]),3)
row2 <- cormat[1,]
rowB <- cormat.pval[1,]
cormat <- round(cor(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt1_24[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt1_24[,PROTS_OF_INTEREST]),3)
row3 <- cormat[1,]
rowC <- cormat.pval[1,]
cormat <- round(cor(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt3_24[,PROTS_OF_INTEREST]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt3_24[,PROTS_OF_INTEREST]),3)
row4 <- cormat[1,]
rowD <- cormat.pval[1,]
CorMat1 <- matrix(c(row1,row2,row3,row4), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),PROTS_OF_INTEREST))
CorMat2 <- matrix(c(rowA,rowB,rowC,rowD), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),PROTS_OF_INTEREST))
plotCor(CorMat1,CorMat2, title = "24 hr", textSize = 4, full = T)
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/Fig2G_CorrelationMatrix_24hr",".pdf"),width = 4.25, height = 4.25)

cormat <- round(cor(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row1 <- cormat[1,]
rowA <- cormat.pval[1,]
cormat <- round(cor(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt0.1_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt0.1_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row2 <- cormat[1,]
rowB <- cormat.pval[1,]
cormat <- round(cor(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt1_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt1_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row3 <- cormat[1,]
rowC <- cormat.pval[1,]
cormat <- round(cor(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt3_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt3_8[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row4 <- cormat[1,]
rowD <- cormat.pval[1,]
CorMat1 <- matrix(c(row1,row2,row3,row4), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")))
CorMat2 <- matrix(c(rowA,rowB,rowC,rowD), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")))
plotCor(CorMat1,CorMat2, title = "8 hr", textSize = 4, full = T)
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS5C_CorrelationMatrix_MDM2probes_8hr",".pdf"),width = 4.25, height = 4.25)

cormat <- round(cor(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row1 <- cormat[1,]
rowA <- cormat.pval[1,]
cormat <- round(cor(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt0.1_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt0.1_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row2 <- cormat[1,]
rowB <- cormat.pval[1,]
cormat <- round(cor(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt1_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt1_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row3 <- cormat[1,]
rowC <- cormat.pval[1,]
cormat <- round(cor(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt3_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),2)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")], y = df_cpt3_24[,c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")]),3)
row4 <- cormat[1,]
rowD <- cormat.pval[1,]
CorMat1 <- matrix(c(row1,row2,row3,row4), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")))
CorMat2 <- matrix(c(rowA,rowB,rowC,rowD), ncol = 4, byrow = T, dimnames = list(c("0 uM","0.1 uM","1 uM","3.3 uM"),c("TP53","MDM2_23384","MDM2_27496","MDM2_27497")))
plotCor(CorMat1,CorMat2, title = "24 hr", textSize = 4, full = T)
ggsave(file = paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS5C_CorrelationMatrix_MDM2probes_24hr",".pdf"),width = 4.25, height = 4.25)

##################################
########## FIGURE 2 E ############
##################################

# 8 hr
df_cor_cpt0.1_8hr <- left_join(df_med_8 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                               df_cpt0.1_8 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                               by = c("EXP_ID","CELL_ID","SAMPLE_CODE","groups_tp53_all"),suffix = c("_med","_cpt"))
df_cor_cpt1_8hr <- left_join(df_med_8 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                             df_cpt1_8 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                             by = c("EXP_ID","CELL_ID","SAMPLE_CODE","groups_tp53_all"),suffix = c("_med","_cpt"))
df_cor_cpt3_8hr <- left_join(df_med_8 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                             df_cpt3_8 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                             by = c("EXP_ID","CELL_ID","SAMPLE_CODE","groups_tp53_all"),suffix = c("_med","_cpt"))

# 24 hr
df_cor_cpt0.1_24hr <- left_join(df_med_24 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                                df_cpt0.1_24 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                                by = c("EXP_ID","CELL_ID","SAMPLE_CODE","groups_tp53_all"),suffix = c("_med","_cpt"))
df_cor_cpt1_24hr <- left_join(df_med_24 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                              df_cpt1_24 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                              by = c("EXP_ID","CELL_ID","SAMPLE_CODE","groups_tp53_all"),suffix = c("_med","_cpt"))
df_cor_cpt3_24hr <- left_join(df_med_24 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                              df_cpt3_24 %>% dplyr::select(-c(DESIGN,CONCENTRATION, TIMEPOINT,CONC_TIME)),
                              by = c("EXP_ID","CELL_ID","SAMPLE_CODE","groups_tp53_all"),suffix = c("_med","_cpt"))

##################################
########## FIGURE S X ############
##################################

# MDM2
ggplot(df_cor_cpt3_8hr) +
  geom_smooth(aes(x = TP53_med,y = MDM2_cpt),method='lm', formula= y~x, se = F, color = 'grey', lty = "dashed") +
  geom_point(aes(x = TP53_med, y = MDM2_cpt, color = groups_tp53_all)) +
  geom_smooth(aes(x = TP53_med,y = MDM2_cpt, color = groups_tp53_all),method='lm', formula= y~x, se = F) +
  annotate(geom="text", x=6.5, y=8.05, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "High") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "High") %>% pull(MDM2_cpt))$estimate,3))) +
  annotate(geom="text", x=5.75, y=8.15, color = "#35B779FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(MDM2_cpt))$estimate,3))) +
  # annotate(geom="text", x=5, y=8.25, color = "#FDE725FF", parse = T,
  #          label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Low") %>% pull(TP53_med), 
  #                                                        y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Low") %>% pull(MDM2_cpt))$estimate,3))) +
  scale_color_viridis(discrete = T, name = "TP53 expression\ncluster") +
  ylab(bquote(atop("MDM2 expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  xlab("Basal TP53 expression") +
  theme_classic() +
  ggtitle("8 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS5D_MDM2_TP53_correlations_8hr.pdf"),width = 4.5, height = 3)
ggplot(df_cor_cpt3_24hr) +
  geom_smooth(aes(x = TP53_med,y = MDM2_cpt),method='lm', formula= y~x, se = F, color = 'grey', lty = "dashed") +
  geom_point(aes(x = TP53_med, y = MDM2_cpt, color = groups_tp53_all)) +
  geom_smooth(aes(x = TP53_med,y = MDM2_cpt, color = groups_tp53_all),method='lm', formula= y~x, se = F) +
  annotate(geom="text", x=6.5, y=9, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "High") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "High") %>% pull(MDM2_cpt))$estimate,3))) +
  annotate(geom="text", x=6, y=9.5, color = "#35B779FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(MDM2_cpt))$estimate,3))) +
  # annotate(geom="text", x=5.5, y=10, color = "#FDE725FF", parse = T,
  #          label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Low") %>% pull(TP53_med), 
  #                                                        y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Low") %>% pull(MDM2_cpt))$estimate,3))) +
  scale_color_viridis(discrete = T, name = "TP53 expression\ncluster") +
  ylab(bquote(atop("MDM2 expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  xlab("Basal TP53 expression") +
  theme_classic() +
  ggtitle("24 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigS5D_MDM2_TP53_correlations_24hr.pdf"),width = 4.5, height = 3)

# CDKN1A
ggplot(df_cor_cpt3_8hr) +
  geom_smooth(aes(x = TP53_med,y = CDKN1A_cpt),method='lm', formula= y~x, se = F, color = 'grey', lty = "dashed") +
  geom_point(aes(x = TP53_med, y = CDKN1A_cpt, color = groups_tp53_all)) +
  geom_smooth(aes(x = TP53_med,y = CDKN1A_cpt, color = groups_tp53_all),method='lm', formula= y~x, se = F) +
  annotate(geom="text", x=6.25, y=6, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "High") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "High") %>% pull(CDKN1A_cpt))$estimate,3))) +
  annotate(geom="text", x=5.5, y=5.5, color = "#35B779FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(CDKN1A_cpt))$estimate,3))) +
  # annotate(geom="text", x=4.75, y=5, color = "#FDE725FF", parse = T,
  #          label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Low") %>% pull(TP53_med), 
  #                                                        y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Low") %>% pull(CDKN1A_cpt))$estimate,3))) +
  scale_color_viridis(discrete = T, name = "TP53 expression\ncluster") +
  ylab(bquote(atop("CDKN1A expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  xlab("Basal TP53 expression") +
  theme_classic() +
  ggtitle("8 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_CDKN1A_TP53_correlations_8hr.pdf"),width = 4.5, height = 3)
ggplot(df_cor_cpt3_24hr) +
  geom_smooth(aes(x = TP53_med,y = CDKN1A_cpt),method='lm', formula= y~x, se = F, color = 'grey', lty = "dashed") +
  geom_point(aes(x = TP53_med, y = CDKN1A_cpt, color = groups_tp53_all)) +
  geom_smooth(aes(x = TP53_med,y = CDKN1A_cpt, color = groups_tp53_all),method='lm', formula= y~x, se = F) +
  annotate(geom="text", x=6.5, y=6, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "High") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "High") %>% pull(CDKN1A_cpt))$estimate,3))) +
  annotate(geom="text", x=6, y=5.5, color = "#35B779FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(CDKN1A_cpt))$estimate,3))) +
  # annotate(geom="text", x=5.5, y=5, color = "#FDE725FF", parse = T,
  #          label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Low") %>% pull(TP53_med), 
  #                                                        y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Low") %>% pull(CDKN1A_cpt))$estimate,3))) +
  scale_color_viridis(discrete = T, name = "TP53 expression\ncluster") +
  ylab(bquote(atop("CDKN1A expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  xlab("Basal TP53 expression") +
  theme_classic() +
  ggtitle("24 hr") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_CDKN1A_TP53_correlations_24hr.pdf"),width = 4.5, height = 3)

# BTG2
ggplot(df_cor_cpt3_8hr) +
  geom_smooth(aes(x = TP53_med,y = BTG2_cpt),method='lm', formula= y~x, se = F, color = 'grey', lty = "dashed") +
  geom_point(aes(x = TP53_med, y = BTG2_cpt, color = groups_tp53_all)) +
  geom_smooth(aes(x = TP53_med,y = BTG2_cpt, color = groups_tp53_all),method='lm', formula= y~x, se = F) +
  annotate(geom="text", x=6.5, y=8.25, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "High") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "High") %>% pull(BTG2_cpt))$estimate,3))) +
  annotate(geom="text", x=5.75, y=8, color = "#35B779FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(BTG2_cpt))$estimate,3))) +
  # annotate(geom="text", x=5, y=7.75, color = "#FDE725FF", parse = T,
  #          label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Low") %>% pull(TP53_med), 
  #                                                        y = df_cor_cpt3_8hr %>% filter(groups_tp53_all == "Low") %>% pull(BTG2_cpt))$estimate,3))) +
  scale_color_viridis(discrete = T, name = "TP53 expression\ncluster") +
  ylab(bquote(atop("BTG2 expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  xlab("Basal TP53 expression") +
  theme_classic() +
  ggtitle("8 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_BTG2_TP53_correlations_8hr.pdf"),width = 4.5, height = 3)
ggplot(df_cor_cpt3_24hr) +
  geom_smooth(aes(x = TP53_med,y = BTG2_cpt),method='lm', formula= y~x, se = F, color = 'grey', lty = "dashed") +
  geom_point(aes(x = TP53_med, y = BTG2_cpt, color = groups_tp53_all)) +
  geom_smooth(aes(x = TP53_med,y = BTG2_cpt, color = groups_tp53_all),method='lm', formula= y~x, se = F) +
  annotate(geom="text", x=6.75, y=6.75, color = "#440154FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "High") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "High") %>% pull(BTG2_cpt))$estimate,3))) +
  annotate(geom="text", x=6, y=6.5, color = "#35B779FF", parse = T,
           label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(TP53_med), 
                                                         y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Intermediate") %>% pull(BTG2_cpt))$estimate,3))) +
  # annotate(geom="text", x=5.25, y=6.25, color = "#FDE725FF", parse = T,
  #          label= sprintf("rho == %0.3f", round(cor.test(x = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Low") %>% pull(TP53_med), 
  #                                                        y = df_cor_cpt3_24hr %>% filter(groups_tp53_all == "Low") %>% pull(BTG2_cpt))$estimate,3))) +
  scale_color_viridis(discrete = T, name = "TP53 expression\ncluster") +
  ylab(bquote(atop("BTG2 expression at",
                   "3.3 "*mu*"M cisplatin"))) +
  xlab("Basal TP53 expression") +
  theme_classic() +
  ggtitle("24 hr") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputDir_GR,DATE,"_TempOSeq_Figs/FigSX_BTG2_TP53_correlations_24hr.pdf"),width = 4.5, height = 3)

######### ######### ######### ######### ######### 
######### XIII. Save workspace image
######### ######### ######### ######### ######### 

save(list = c("outputDir_GR","DATE",
              "PROTS_OF_INTEREST",
              "count_Filtered", "meta_Filtered",
              "Norm","Normt",
              "hmdata","hmdataMeta","bpdata_phh","bpdata_8hr",
              "df_med_8","df_med_24",
              "df_cpt0.1_8","df_cpt1_8","df_cpt3_8",
              "df_cpt0.1_24","df_cpt1_24","df_cpt3_24"),
     file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/",DATE,"workspace.RData")) 

#save.image(file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/",DATE,"complete_workspace.RData"))
#load(file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/20220321complete_workspace.RData"))
