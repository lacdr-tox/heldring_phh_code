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

rm(list = ls())

load("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/20220321workspace.RData")
DATE <- "20220411_M019"
dir.create(paste0(outputDir_GR,DATE,"_Correlation_Figs"))

## Functions ##

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

# Bootstrap function
correlation <- function(data, indices, probes, decims = 2,  
                        rows_select = "TP53_7287",
                        cols_select = c("TP53_7287","MDM2_sum","CDKN1A_1219","BTG2_13191")) {
  d <- data[indices,probes] # allows boot to select sample
  cormat <- round(cor(d),decims)
  return(cormat[rows_select,cols_select])
}

makeCorGenData <- function(genData, withRNA, varyR, mdm2RNA) {
  # Do correlation analysis per set
  if (withRNA & varyR & (!mdm2RNA)) {
    # Option 1
    print("corGenData Option 1")
    corGenData <- genData %>% group_by(run, rfactor, variation) %>% 
      summarise(
        # 8 hr
        MDM2_med_vs_TP53_med = cor(x = MDM2rna_med, y = TP53_med),
        CDKN1A_med_vs_TP53_med = cor(x = CDKN1A_med, y = TP53_med),
        BTG2_med_vs_TP53_med = cor(x = BTG2rna_med, y = TP53_med),
        MDM2_cpt8_vs_TP53_cpt8 = cor(x = MDM2rna_cpt8, y = TP53_cpt8),
        CDKN1A_cpt8_vs_TP53_cpt8 = cor(x = CDKN1A_cpt8, y = TP53_cpt8),
        BTG2_cpt8_vs_TP53_cpt8 = cor(x = BTG2rna_cpt8, y = TP53_cpt8),
        TP53_cpt8_vs_TP53_med = cor(x = TP53_cpt8, y = TP53_med),
        MDM2_cpt8_vs_TP53_med = cor(x = MDM2rna_cpt8, y = TP53_med),
        CDKN1A_cpt8_vs_TP53_med = cor(x = CDKN1A_cpt8, y = TP53_med),
        BTG2_cpt8_vs_TP53_med = cor(x = BTG2rna_cpt8, y = TP53_med),
        # 24 hr
        MDM2_cpt24_vs_TP53_cpt24 = cor(x = MDM2rna_cpt24, y = TP53_cpt24),
        CDKN1A_cpt24_vs_TP53_cpt24 = cor(x = CDKN1A_cpt24, y = TP53_cpt24),
        BTG2_cpt24_vs_TP53_cpt24 = cor(x = BTG2rna_cpt24, y = TP53_cpt24),
        TP53_cpt24_vs_TP53_med = cor(x = TP53_cpt24, y = TP53_med),
        MDM2_cpt24_vs_TP53_med = cor(x = MDM2rna_cpt24, y = TP53_med),
        CDKN1A_cpt24_vs_TP53_med = cor(x = CDKN1A_cpt24, y = TP53_med),
        BTG2_cpt24_vs_TP53_med = cor(x = BTG2rna_cpt24, y = TP53_med))
  } else if (withRNA & (!varyR) & mdm2RNA) {
    # Option 2
    print("corGenData Option 2")
    corGenData <- genData %>% group_by(run, variation) %>% 
      summarise(
        # 8 hr 
        MDM2_med_vs_TP53_med = cor(x = MDM2rna_med, y = TP53_med),
        MDM2_med_vs_TP53_med = cor(x = MDM2rna_med, y = TP53_med),
        CDKN1A_med_vs_TP53_med = cor(x = CDKN1A_med, y = TP53_med),
        BTG2_med_vs_TP53_med = cor(x = BTG2rna_med, y = TP53_med),
        MDM2_cpt8_vs_TP53_cpt8 = cor(x = MDM2rna_cpt8, y = TP53_cpt8),
        MDM2_cpt8_vs_TP53_cpt8 = cor(x = MDM2rna_cpt8, y = TP53_cpt8),
        CDKN1A_cpt8_vs_TP53_cpt8 = cor(x = CDKN1A_cpt8, y = TP53_cpt8),
        BTG2_cpt8_vs_TP53_cpt8 = cor(x = BTG2rna_cpt8, y = TP53_cpt8),
        TP53_cpt8_vs_TP53_med = cor(x = TP53_cpt8, y = TP53_med),
        MDM2_cpt8_vs_TP53_med = cor(x = MDM2rna_cpt8, y = TP53_med),
        MDM2_cpt8_vs_TP53_med = cor(x = MDM2rna_cpt8, y = TP53_med),
        CDKN1A_cpt8_vs_TP53_med = cor(x = CDKN1A_cpt8, y = TP53_med),
        BTG2_cpt8_vs_TP53_med = cor(x = BTG2rna_cpt8, y = TP53_med),
        # 24 hr
        MDM2_cpt24_vs_TP53_cpt24 = cor(x = MDM2rna_cpt24, y = TP53_cpt24),
        MDM2_cpt24_vs_TP53_cpt24 = cor(x = MDM2rna_cpt24, y = TP53_cpt24),
        CDKN1A_cpt24_vs_TP53_cpt24 = cor(x = CDKN1A_cpt24, y = TP53_cpt24),
        BTG2_cpt24_vs_TP53_cpt24 = cor(x = BTG2rna_cpt24, y = TP53_cpt24),
        TP53_cpt24_vs_TP53_med = cor(x = TP53_cpt24, y = TP53_med),
        MDM2_cpt24_vs_TP53_med = cor(x = MDM2rna_cpt24, y = TP53_med),
        MDM2_cpt24_vs_TP53_med = cor(x = MDM2rna_cpt24, y = TP53_med),
        CDKN1A_cpt24_vs_TP53_med = cor(x = CDKN1A_cpt24, y = TP53_med),
        BTG2_cpt24_vs_TP53_med = cor(x = BTG2rna_cpt24, y = TP53_med))
  } else if (withRNA & (!varyR) & (!mdm2RNA)) {
    # Option 3
    print("corGenData Option 3")
    corGenData <- genData %>% group_by(run, variation) %>% 
      summarise(
        # 8 hr
        MDM2_med_vs_TP53_med = cor(x = MDM2rna_med, y = TP53_med),
        CDKN1A_med_vs_TP53_med = cor(x = CDKN1A_med, y = TP53_med),
        BTG2_med_vs_TP53_med = cor(x = BTG2rna_med, y = TP53_med),
        MDM2_cpt8_vs_TP53_cpt8 = cor(x = MDM2rna_cpt8, y = TP53_cpt8),
        CDKN1A_cpt8_vs_TP53_cpt8 = cor(x = CDKN1A_cpt8, y = TP53_cpt8),
        BTG2_cpt8_vs_TP53_cpt8 = cor(x = BTG2rna_cpt8, y = TP53_cpt8),
        TP53_cpt8_vs_TP53_med = cor(x = TP53_cpt8, y = TP53_med),
        MDM2_cpt8_vs_TP53_med = cor(x = MDM2rna_cpt8, y = TP53_med),
        CDKN1A_cpt8_vs_TP53_med = cor(x = CDKN1A_cpt8, y = TP53_med),
        BTG2_cpt8_vs_TP53_med = cor(x = BTG2rna_cpt8, y = TP53_med),
        # 24 hr
        MDM2_cpt24_vs_TP53_cpt24 = cor(x = MDM2rna_cpt24, y = TP53_cpt24),
        CDKN1A_cpt24_vs_TP53_cpt24 = cor(x = CDKN1A_cpt24, y = TP53_cpt24),
        BTG2_cpt24_vs_TP53_cpt24 = cor(x = BTG2rna_cpt24, y = TP53_cpt24),
        TP53_cpt24_vs_TP53_med = cor(x = TP53_cpt24, y = TP53_med),
        MDM2_cpt24_vs_TP53_med = cor(x = MDM2rna_cpt24, y = TP53_med),
        CDKN1A_cpt24_vs_TP53_med = cor(x = CDKN1A_cpt24, y = TP53_med),
        BTG2_cpt24_vs_TP53_med = cor(x = BTG2rna_cpt24, y = TP53_med))
  } else if ((!withRNA) & (!varyR) & (!mdm2RNA)) {
    # Option 4
    print("corGenData Option 4")
    corGenData <- genData %>% group_by(run, variation) %>% 
      summarise(MDM2_med_vs_p53_med = cor(x = MDM2_med, y = p53_total_med),
                p21_med_vs_p53_med = cor(x = p21_med, y = p53_total_med),
                BTG2_med_vs_p53_med = cor(x = BTG2_med, y = p53_total_med),
                MDM2_cpt_vs_p53_cpt = cor(x = MDM2_cpt, y = p53_total_cpt),
                p21_cpt_vs_p53_cpt = cor(x = p21_cpt, y = p53_total_cpt),
                BTG2_cpt_vs_p53_cpt = cor(x = BTG2_cpt, y = p53_total_cpt),
                p53_cpt_vs_p53_med = cor(x = p53_total_cpt, y = p53_total_med),
                MDM2_cpt_vs_p53_med = cor(x = MDM2_cpt, y = p53_total_med),
                p21_cpt_vs_p53_med = cor(x = p21_cpt, y = p53_total_med),
                BTG2_cpt_vs_p53_med = cor(x = BTG2_cpt, y = p53_total_med))
  } else {
    # Option 5
    warning("Option 5 is not provided. No corGenData created.")
  }
  return(corGenData)
}

changeCols <- function(genData){
  newColNames <- str_replace(
    str_replace(
      str_replace(
        str_replace(colnames(genData), "_init", "_med"), "_8hr", "_cpt8"), "_24hr", "_cpt24"), "P", "p")
  colnames(genData) <- newColNames
  newColNames <- str_replace(
    str_replace(colnames(genData), "p53rna", "TP53"), "p21rna", "CDKN1A")
  colnames(genData) <- newColNames
  return(genData)
}

addNoise <- function(genData, withNoise, groupby, noise = 0.125) {
  if (withNoise){
    genDataNoise <- genData %>% group_by({{groupby}}) %>%
      mutate(TP53_med = exp(rnorm(n(),log(TP53_med),noise)),
             MDM2_med = exp(rnorm(n(),log(MDM2rna_med),noise)),
             CDKN1A_med = exp(rnorm(n(),log(CDKN1A_med),noise)),
             BTG2_med = exp(rnorm(n(),log(BTG2rna_med),noise)),
             
             TP53_cpt8 = exp(rnorm(n(),log(TP53_cpt8),noise)),
             MDM2_cpt8 = exp(rnorm(n(),log(MDM2rna_cpt8),noise)),
             CDKN1A_cpt8 = exp(rnorm(n(),log(CDKN1A_cpt8),noise)),
             BTG2_cpt8 = exp(rnorm(n(),log(BTG2rna_cpt8),noise)),
             
             TP53_cpt24 = exp(rnorm(n(),log(TP53_cpt24), noise)),
             MDM2_cpt24 = exp(rnorm(n(),log(MDM2rna_cpt24), noise)),
             CDKN1A_cpt24 = exp(rnorm(n(),log(CDKN1A_cpt24), noise)),
             BTG2_cpt24 = exp(rnorm(n(),log(BTG2rna_cpt24), noise)))
    genData <- genDataNoise %>% ungroup()
  }
  return(genData)
}

plottingDataFrames <- function(corGenData, mdm2RNA, withRNA){
  # Create dataframes for plotting
  medVSmed <- filter_at(corGenData, "name", any_vars(str_ends(.,"_med_vs_TP53_med")))
  cptVSmed8 <- filter_at(corGenData, "name", any_vars(str_ends(.,"_cpt8_vs_TP53_med")))
  cptVScpt8 <- filter_at(corGenData, "name", any_vars(str_ends(.,"_cpt8_vs_TP53_cpt8")))
  
  cptVSmed24 <- filter_at(corGenData, "name", any_vars(str_ends(.,"_cpt24_vs_TP53_med")))
  cptVScpt24 <- filter_at(corGenData, "name", any_vars(str_ends(.,"_cpt24_vs_TP53_cpt24")))
  
  # Make protein order
  if (mdm2RNA) {
    medVSmed <- mutate_at(medVSmed, vars("axisLab"), factor, c("MDM2rna med vs.\nTP53 med", "MDM2 med vs.\nTP53 med", "CDKN1A med vs.\nTP53 med", "BTG2 med vs.\nTP53 med"))
    cptVSmed8 <- mutate_at(cptVSmed8, vars("axisLab"), factor, c("TP53 cpt8 vs.\nTP53 med", "MDM2rna cpt8 vs.\nTP53 med", "MDM2 cpt8 vs.\nTP53 med", "CDKN1A cpt8 vs.\nTP53 med", "BTG2 cpt8 vs.\nTP53 med"))
    cptVScpt8 <- mutate_at(cptVScpt8, vars("axisLab"), factor, c("MDM2rna cpt8 vs.\nTP53 cpt8", "MDM2 cpt8 vs.\nTP53 cpt8", "CDKN1A cpt8 vs.\nTP53 cpt8", "BTG2 cpt8 vs.\nTP53 cpt8"))
    
    cptVSmed24 <- mutate_at(cptVSmed24, vars("axisLab"), factor, c("TP53 cpt24 vs.\nTP53 med", "MDM2rna cpt24 vs.\nTP53 med", "MDM2 cpt24 vs.\nTP53 med", "CDKN1A cpt24 vs.\nTP53 med", "BTG2 cpt24 vs.\nTP53 med"))
    cptVScpt24 <- mutate_at(cptVScpt24, vars("axisLab"), factor, c("MDM2rna cpt24 vs.\nTP53 cpt24", "MDM2 cpt24 vs.\nTP53 cpt24", "CDKN1A cpt24 vs.\nTP53 cpt24", "BTG2 cpt24 vs.\nTP53 cpt24"))
  } else if (!withRNA) {
    medVSmed <- mutate_at(medVSmed, vars("axisLab"), factor, c("MDM2 med vs.\np53 med", "p21 med vs.\np53 med", "BTG2 med vs.\np53 med"))
    cptVSmed8 <- mutate_at(cptVSmed8, vars("axisLab"), factor, c("p53 cpt8 vs.\np53 med", "MDM2 cpt8 vs.\np53 med", "p21 cpt8 vs.\np53 med", "BTG2 cpt8 vs.\np53 med"))
    cptVScpt8 <- mutate_at(cptVScpt8, vars("axisLab"), factor, c("MDM2 cpt8 vs.\np53 cpt8", "p21 cpt8 vs.\np53 cpt8", "BTG2 cpt8 vs.\np53 cpt8"))
    
    cptVSmed24 <- mutate_at(cptVSmed24, vars("axisLab"), factor, c("p53 cpt24 vs.\np53 med", "MDM2 cpt24 vs.\np53 med", "p21 cpt24 vs.\np53 med", "BTG2 cpt24 vs.\np53 med"))
    cptVScpt24 <- mutate_at(cptVScpt24, vars("axisLab"), factor, c("MDM2 cpt24 vs.\np53 cpt24", "p21 cpt24 vs.\np53 cpt24", "BTG2 cpt24 vs.\np53 cpt24"))
    
    medVSmed <- mutate_at(medVSmed, vars("axisLab2"), factor, c("MDM2", "p21", "BTG2"))
    cptVSmed8 <- mutate_at(cptVSmed8, vars("axisLab2"), factor, c("p53", "MDM2", "p21", "BTG2"))
    cptVScpt8 <- mutate_at(cptVScpt8, vars("axisLab2"), factor, c("MDM2", "p21", "BTG2"))
    
    cptVSmed24 <- mutate_at(cptVSmed24, vars("axisLab2"), factor, c("p53", "MDM2", "p21", "BTG2"))
    cptVScpt24 <- mutate_at(cptVScpt24, vars("axisLab2"), factor, c("MDM2", "p21", "BTG2"))
  } else if (withRNA) {
    medVSmed <- mutate_at(medVSmed, vars("axisLab"), factor, c("MDM2 med vs.\nTP53 med", "CDKN1A med vs.\nTP53 med", "BTG2 med vs.\nTP53 med"))
    cptVSmed8 <- mutate_at(cptVSmed8, vars("axisLab"), factor, c("TP53 cpt8 vs.\nTP53 med", "MDM2 cpt8 vs.\nTP53 med", "CDKN1A cpt8 vs.\nTP53 med", "BTG2 cpt8 vs.\nTP53 med"))
    cptVScpt8 <- mutate_at(cptVScpt8, vars("axisLab"), factor, c("MDM2 cpt8 vs.\nTP53 cpt8", "CDKN1A cpt8 vs.\nTP53 cpt8", "BTG2 cpt8 vs.\nTP53 cpt8"))
    
    cptVSmed24 <- mutate_at(cptVSmed24, vars("axisLab"), factor, c("TP53 cpt24 vs.\nTP53 med", "MDM2 cpt24 vs.\nTP53 med", "CDKN1A cpt24 vs.\nTP53 med", "BTG2 cpt24 vs.\nTP53 med"))
    cptVScpt24 <- mutate_at(cptVScpt24, vars("axisLab"), factor, c("MDM2 cpt24 vs.\nTP53 cpt24", "CDKN1A cpt24 vs.\nTP53 cpt24", "BTG2 cpt24 vs.\nTP53 cpt24"))
    
    medVSmed <- mutate_at(medVSmed, vars("axisLab2"), factor, c("MDM2", "CDKN1A", "BTG2"))
    cptVSmed8 <- mutate_at(cptVSmed8, vars("axisLab2"), factor, c("TP53", "MDM2", "CDKN1A", "BTG2"))
    cptVScpt8 <- mutate_at(cptVScpt8, vars("axisLab2"), factor, c("MDM2", "CDKN1A", "BTG2"))
    
    cptVSmed24 <- mutate_at(cptVSmed24, vars("axisLab2"), factor, c("TP53", "MDM2", "CDKN1A", "BTG2"))
    cptVScpt24 <- mutate_at(cptVScpt24, vars("axisLab2"), factor, c("MDM2", "CDKN1A", "BTG2"))
    
  }
  return(list("medVSmed" = medVSmed, "cptVSmed8" = cptVSmed8, "cptVScpt8" = cptVScpt8, "cptVSmed24" = cptVSmed24, "cptVScpt24" = cptVScpt24))
}

######### ######### ######### ######### ######### ######### 
######### ######### Experimental data # ######### #########
######### ######### ######### ######### ######### ######### 

# Do for 8 hr
## Make correlation plots for counts data (medium vs. medium)
cormat <- round(cor(df_med_8[,PROTS_OF_INTEREST]),3)
cormat.pval <- round(cor.test.p(df_med_8[,PROTS_OF_INTEREST]),3)
medVSmedReal_8 <- cormat["TP53",c("MDM2","CDKN1A","BTG2")]

## Make correlation plots for counts data (cisplatin 3.3 uM vs. cisplatin 3.3 uM)
cormat <- round(cor(df_cpt3_8[,PROTS_OF_INTEREST]),3)
cormat.pval <- round(cor.test.p(df_cpt3_8[,PROTS_OF_INTEREST]),3)
cpt3VScpt3Real_8 <- cormat["TP53",c("MDM2","CDKN1A","BTG2")]

## Make correlation plots for counts data (cisplatin 3 uM vs. medium)
cormat <- round(cor(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt3_8[,PROTS_OF_INTEREST]),3)
cormat.pval <- round(cor.test.p.xy(x = df_med_8[,PROTS_OF_INTEREST], y = df_cpt3_8[,PROTS_OF_INTEREST]),3)
cpt3VSmedReal_8 <- cormat["TP53",c("TP53", "MDM2","CDKN1A","BTG2")]

# Do for 24 hr
## Make correlation plots for counts data (medium vs. medium)
cormat <- round(cor(df_med_24[,PROTS_OF_INTEREST]),3)
cormat.pval <- round(cor.test.p(df_med_24[,PROTS_OF_INTEREST]),3)
medVSmedReal_24 <- cormat["TP53",c("MDM2","CDKN1A","BTG2")]

## Make correlation plots for counts data (cisplatin 3.3 uM vs. cisplatin 3.3 uM)
cormat <- round(cor(df_cpt3_24[,PROTS_OF_INTEREST]),3)
cormat.pval <- round(cor.test.p(df_cpt3_24[,PROTS_OF_INTEREST]),3)
cpt3VScpt3Real_24 <- cormat["TP53",c("MDM2","CDKN1A","BTG2")]

## Make correlation plots for counts data (cisplatin 3 uM vs. medium)
cormat <- round(cor(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt3_24[,PROTS_OF_INTEREST]),3)
cormat.pval <- round(cor.test.p.xy(x = df_med_24[,PROTS_OF_INTEREST], y = df_cpt3_24[,PROTS_OF_INTEREST]),3)
cpt3VSmedReal_24 <- cormat["TP53",c("TP53", "MDM2","CDKN1A","BTG2")]


######### ######### ######### #########
#### Get available data files #########
######### ######### ######### #########

genDataInputDir <- "/data/muriel/Projects/PHH/Modeling/Output/VaryParmsData/"
genDataFiles <- dir(genDataInputDir)
genDataFiles

######### ######### ######### ######### ######### ######### 
#### 0a. Check out parameter variations, without noise ####
######### ######### ######### ######### ######### ######### 

# Get the relevant data
# used for 1st manuscript:"20201116" & "20201201_varyR"
# used for 1st manuscript to co-authors: "20210129_varyR" & "20201114" & "M022"
# used for 1st new run after manuscript came back, negative values due to r=100: "20210331_varyparms_varyR" & "20210329" & "M029"
# used for run with 
variableDate <- "20220322_varyR_VaryParameters" 
modelDate <- "20220314_111222" 
MODEL <- "M019Model_parmset3"

genDataFiles <- genDataFiles[sapply(genDataFiles, grepl, pattern = modelDate)]
genDataFiles

files_tmp <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_w"))] 
files_tmp
parm_pert_file <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_p"))] 
parm_pert_file

withRNA <- T
withNoise <- F
mdm2RNA <- F
varyR <- F

# Printing for information
print(paste("Working on plots for file", paste(files_tmp, sep = " ", collapse = ", ")))

genData_orig <- fread(paste0(genDataInputDir,files_tmp))
pert_orig <- fread(paste0(genDataInputDir,parm_pert_file))

# Check for negative parameter values
# Check for negative parameter values
pert <- pert_orig %>% filter(!ParName %in% c("sf_p53","sf_mdm2","sf_p21","sf_btg2",
                                             "offset_p53","offset_mdm2","offset_p21","offset_btg2"))
cnames <- pert %>% pull(ParName)
pert <- as_tibble(t(pert %>% dplyr::select(-c(V1, ParName))))
colnames(pert) <- cnames
if (min(pert) < 0) {
  if (min(pert %>% dplyr::select(contains("_init"))) < 0) {
    warning("Oops! Negative values detected in steady states!")
    print(min(pert %>% dplyr::select(contains("_init"))))
  } 
  if (min(pert %>% dplyr::select(contains("k"))) < 0) {
    warning("Oops! Negative values detected in parameters!")
    print(min(pert %>% dplyr::select(contains("k"))))
  }
}

# Make tibble data
genData <- as_tibble(genData_orig)

# Change column names
genData <- changeCols(genData)

# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

minGenData <- genData %>% dplyr::select(-c(variation, run, sample)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Add noise
noise = 0.125
genData <- addNoise(genData, withNoise, variation, noise = noise)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(variation, run, sample)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)


corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)

# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

legendlabs <- c("0.001" = "0.001", 
                "0.01" = "0.01",
                "0.05" = "0.05",
                "0.1" = "0.1",
                "0.2" = "0.2",
                "0.3" = "0.3")

if (mdm2RNA) {
  medVSmed <- medVSmed %>% filter(!name == "MDM2_med_vs_p53_med")
  cptVSmed8 <- cptVSmed8 %>% filter(!name == "MDM2_cpt8_vs_p53_med")
}

if (varyR) {
  medVSmed <- medVSmed %>% filter(rfactor == 1)
  cptVSmed8 <- cptVSmed8 %>% filter(rfactor == 1)
  cptVScpt8 <- cptVScpt8 %>% filter(rfactor == 1)
  
  cptVSmed24 <- cptVSmed24 %>% filter(rfactor == 1)
  cptVScpt24 <- cptVScpt24 %>% filter(rfactor == 1)
}

# Do for 8 hours
# bootstrapping with 1000 replicates
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = c("MDM2","CDKN1A","BTG2")) #"TP53", "MDM2","CDKN1A","BTG2"
ci_mdm2 <- boot.ci(boot.out = results_med8, 
                   index = c(1,1),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med8, 
                   index = c(2,2),
                   type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med8, 
                   index = c(3,3),
                   type = c("norm"))

res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med", "CDKN1Amed_vs_TP53_med", "BTG2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med", "CDKN1Amed_vs_TP53_med", "BTG2_med_vs_TP53_med")),
                          mean = c(ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))
  
# Make a boxplot
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8) 
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
ggplot() +
  geom_blank(data = medVSmed , aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_med,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = medVSmed, aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle("Medium, 8 hr") + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10A_genData_boxplot_medium_vs_medium_8_",variableDate,".pdf"),
           width = 5, height = 3)


# bootstrapping with 1000 replicates
# Cisplatin 3.3 uM vs. medium
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("TP53_cpt3","MDM2_cpt3","CDKN1A_cpt3","BTG2_cpt3"))

ci_p53 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(1,1),
                   type = c("norm"))
ci_mdm2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(2,2),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med_cpt3, 
                  index = c(3,3),
                  type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(4,4),
                   type = c("norm"))

res_sum_cpt <- data.frame(name = factor(c("TP53_cpt_vs_TP53_med","MDM2_cpt_vs_TP53_med", "CDKN1Acpt_vs_TP53_med", "BTG2_cpt_vs_TP53_med"),
                                        levels = c("TP53_cpt_vs_TP53_med","MDM2_cpt_vs_TP53_med", "CDKN1Acpt_vs_TP53_med", "BTG2_cpt_vs_TP53_med")),
                          mean = c(ci_p53$t0,ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_p53$normal[2],ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_p53$normal[3],ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))

hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)
ggplot() +
  geom_blank(data = cptVSmed8, aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_cpt,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_cpt, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_cpt, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = cptVSmed8, aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs) + 
  ggtitle("Cisplatin, 8hr")
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10A_genData_boxplot_cisplatin3.3uM_vs_medium_8_",variableDate,".pdf"),
           width = 5, height = 3)


# Do for 24 hours
# bootstrapping with 1000 replicates
results_med24 <- boot(data=df_med_24, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = c("MDM2","CDKN1A","BTG2")) 

ci_mdm2 <- boot.ci(boot.out = results_med24, 
                   index = c(1,1),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med24, 
                  index = c(2,2),
                  type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med24, 
                   index = c(3,3),
                   type = c("norm"))

res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med")),
                          mean = c(ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))

# Make a boxplot
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24) 
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
ggplot() +
  geom_blank(data = medVSmed , aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_med,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = medVSmed, aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle("Medium, 24 hr") + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigSX_genData_boxplot_medium_vs_medium_24_",variableDate,".pdf"),
           width = 5, height = 3)


# bootstrapping with 1000 replicates
# Cisplatin 3.3 uM vs. medium
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("TP53_cpt3","MDM2_cpt3","CDKN1A_cpt3","BTG2_cpt3"))

ci_p53 <- boot.ci(boot.out = results_med_cpt3, 
                  index = c(1,1),
                  type = c("norm"))
ci_mdm2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(2,2),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med_cpt3, 
                  index = c(3,3),
                  type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(4,4),
                   type = c("norm"))

res_sum_cpt <- data.frame(name = factor(c("TP53_cpt_vs_TP53_med","MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med"),
                                        levels = c("TP53_cpt_vs_TP53_med","MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med")),
                          mean = c(ci_p53$t0,ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_p53$normal[2],ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_p53$normal[3],ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))


hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)
ggplot() +
  geom_blank(data = cptVSmed24, aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_cpt,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_cpt, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  geom_boxplot(data = cptVSmed24, aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs) + 
  ggtitle("Cisplatin, 24 hr")
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigSX_genData_boxplot_cisplatin3.3uM_vs_medium_24_",variableDate,".pdf"),
           width = 5, height = 3)

######### ######### ######### ######### ######### ######### 
##### Save data for low, medium and high correlations  #### 
######### ######### ######### ######### ######### ######### 

if (varyR) {
  # Export a low, a medium and a high correlation group
  low <- corGenData %>% group_by(name) %>%
    filter(rfactor == 1,
           variation == 0.2 &
             str_detect(name, "vs_TP53_med")) %>%
    filter(str_starts(name,"TP53", negate = T)) %>%
    filter(value == min(value)) %>% mutate(correlation = "low")
  med <- corGenData %>% group_by(name) %>%
    filter(rfactor == 1,
           variation == 0.2 &
             str_detect(name, "vs_TP53_med")) %>%
    filter(str_starts(name,"TP53", negate = T)) %>%
    filter(abs(value - mean(value)) == min(abs(value - mean(value)))) %>% mutate(correlation = "med")
  high <- corGenData %>% group_by(name) %>%
    filter(rfactor == 1,
           variation == 0.2 &
             str_detect(name, "vs_TP53_med")) %>%
    filter(str_starts(name,"TP53", negate = T)) %>%
    filter(value == max(value)) %>% mutate(correlation = "high")
  all <- bind_rows(list(low,med,high))
  write.table(all,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_CorrelationSetSelection_",variableDate,".csv"))
} else if (!varyR) {
  # Export a low, a medium and a high correlation group
  low <- corGenData %>% group_by(name) %>%
    filter(variation == 0.2 &
             str_detect(name, "vs_TP53_med")) %>%
    filter(str_starts(name,"TP53", negate = T)) %>%
    filter(value == min(value)) %>% mutate(correlation = "low")
  med <- corGenData %>% group_by(name) %>%
    filter(variation == 0.2 &
             str_detect(name, "vs_TP53_med")) %>%
    filter(str_starts(name,"TP53", negate = T)) %>%
    filter(abs(value - mean(value)) == min(abs(value - mean(value)))) %>% mutate(correlation = "med")
  high <- corGenData %>% group_by(name) %>%
    filter(variation == 0.2 &
             str_detect(name, "vs_TP53_med")) %>%
    filter(str_starts(name,"TP53", negate = T)) %>%
    filter(value == max(value)) %>% mutate(correlation = "high")
  all <- bind_rows(list(low,med,high))
  write.table(all,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_CorrelationSetSelection_",variableDate,".csv"))
}

# Save parameter sets low correlations, 8 hr
low_8 <- low %>% filter(str_detect(name,"cpt8"))
runs_low <- unique(paste0("r",low_8$run,"_"))
vars_low <- unique(paste0("_v",low_8$variation))
pert_select_low_8 <- pert_orig %>% dplyr::select(c(ParName,est_value,starts_with(runs_low) & ends_with("rfactor1") & contains(vars_low)))
write.table(pert_select_low_8,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_ParmPertSelection_low_8hr_",variableDate,".csv"))
# Save parameter sets low correlations, 24 hr
low_24 <- low %>% filter(str_detect(name,"cpt24"))
runs_low <- unique(paste0("r",low_24$run,"_"))
vars_low <- unique(paste0("_v",low_24$variation))
pert_select_low_24 <- pert_orig %>% dplyr::select(c(ParName,est_value,starts_with(runs_low) & ends_with("rfactor1") & contains(vars_low)))
write.table(pert_select_low_24,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_ParmPertSelection_low_24hr_",variableDate,".csv"))

# Save parameter sets med correlations, 8 hr
med_8 <- med %>% filter(str_detect(name,"cpt8"))
runs_med <- unique(paste0("r",med_8$run,"_"))
vars_med <- unique(paste0("_v",med_8$variation))
pert_select_med_8 <- pert_orig %>% dplyr::select(c(ParName,est_value,starts_with(runs_med) & ends_with("rfactor1") & contains(vars_med)))
write.table(pert_select_med_8,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_ParmPertSelection_med_8hr_",variableDate,".csv"))
# Save parameter sets med correlations, 24 hr
med_24 <- med %>% filter(str_detect(name,"cpt24"))
runs_med <- unique(paste0("r",med_24$run,"_"))
vars_med <- unique(paste0("_v",med_24$variation))
pert_select_med_24 <- pert_orig %>% dplyr::select(c(ParName,est_value,starts_with(runs_med) & ends_with("rfactor1") & contains(vars_med)))
write.table(pert_select_med_24,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_ParmPertSelection_med_24hr_",variableDate,".csv"))

# Save parameter sets high correlations, 8 hr
high_8 <- high %>% filter(str_detect(name,"cpt8"))
runs_high <- unique(paste0("r",high_8$run,"_"))
vars_high <- unique(paste0("_v",high_8$variation))
pert_select_high_8 <- pert_orig %>% dplyr::select(c(ParName,est_value,starts_with(runs_high) & ends_with("rfactor1") & contains(vars_high)))
write.table(pert_select_high_8,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_ParmPertSelection_high_8hr_",variableDate,".csv"))
# Save parameter sets high correlations, 24 hr
high_24 <- high %>% filter(str_detect(name,"cpt24"))
runs_high <- unique(paste0("r",high_24$run,"_"))
vars_high <- unique(paste0("_v",high_24$variation))
pert_select_high_24 <- pert_orig %>% dplyr::select(c(ParName,est_value,starts_with(runs_high) & ends_with("rfactor1") & contains(vars_high)))
write.table(pert_select_high_24,file = paste0("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/VaryR_ParmSets/",DATE,"_",MODEL,"_ParmPertSelection_high_24hr_",variableDate,".csv"))

######### ######### ######### ######### ######### ######### 
#### 0b. Check out parameter variations, with noise #######
######### ######### ######### ######### ######### ######### 

withNoise <- T

# Make tibble data
genData <- as_tibble(genData_orig)

# Change column names
genData <- changeCols(genData)
colnames(genData)

# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

minGenData <- genData %>% dplyr::select(-c(variation, run, sample)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

noise = 0.125
genData <- addNoise(genData, withNoise, variation, noise = noise)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(variation, run, sample)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)


# Do correlation analysis per set
corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)

# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

# Make data for plotting
l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

legendlabs <- c("0.001" = "0.001", 
                "0.01" = "0.01",
                "0.05" = "0.05",
                "0.1" = "0.1",
                "0.2" = "0.2",
                "0.3" = "0.3")


# Do for 8 hours
# bootstrapping with 1000 replicates
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = c("MDM2","CDKN1A","BTG2")) #"TP53", "MDM2","CDKN1A","BTG2"
#plot(results_med8)
ci_mdm2 <- boot.ci(boot.out = results_med8, 
                   index = c(1,1),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med8, 
                  index = c(2,2),
                  type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med8, 
                   index = c(3,3),
                   type = c("norm"))

res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med")),
                          mean = c(ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))

# Make a boxplot
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8) 
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
ggplot() +
  geom_blank(data = medVSmed %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_med,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = medVSmed %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.7) + 
  theme_classic() + 
  ggtitle("Medium, 8 hr") + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/Fig5B_genData_boxplot_withNoise",noise,"_medium_vs_medium_8_",variableDate,".pdf"),
         width = 3, height = 3)


# bootstrapping with 1000 replicates
# Cisplatin 3.3 uM vs. medium
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("TP53_cpt3","MDM2_cpt3","CDKN1A_cpt3","BTG2_cpt3"))

ci_p53 <- boot.ci(boot.out = results_med_cpt3, 
                  index = c(1,1),
                  type = c("norm"))
ci_mdm2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(2,2),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med_cpt3, 
                  index = c(3,3),
                  type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(4,4),
                   type = c("norm"))

res_sum_cpt <- data.frame(name = factor(c("TP53_cpt_vs_TP53_med","MDM2_cpt_vs_TP53_med", "CDKN1A_cpt_vs_TP53_med", "BTG2_cpt_vs_TP53_med"),
                                        levels = c("TP53_cpt_vs_TP53_med","MDM2_cpt_vs_TP53_med", "CDKN1A_cpt_vs_TP53_med", "BTG2_cpt_vs_TP53_med")),
                          mean = c(ci_p53$t0,ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_p53$normal[2],ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_p53$normal[3],ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))

hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)
ggplot() +
  geom_blank(data = cptVSmed8 %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_cpt,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_cpt, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_cpt, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = cptVSmed8 %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.7) + 
  theme_classic() + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs) + 
  ggtitle("Cisplatin, 8hr")
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/Fig5B_genData_boxplot_withNoise",noise,"_cisplatin3.3uM_vs_medium_8_",variableDate,".pdf"),
         width = 3.5, height = 3)

# Do for 24 hours
# bootstrapping with 1000 replicates
results_med24 <- boot(data=df_med_24, statistic=correlation,
                      R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                      rows_select = "TP53",
                      cols_select = c("MDM2","CDKN1A","BTG2")) #"TP53", "MDM2","CDKN1A","BTG2"

ci_mdm2 <- boot.ci(boot.out = results_med24, 
                   index = c(1,1),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med24, 
                  index = c(2,2),
                  type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med24, 
                   index = c(3,3),
                   type = c("norm"))

res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med")),
                          mean = c(ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))

# Make a boxplot
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24) 
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
ggplot() +
  geom_blank(data = medVSmed  %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_med,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_med, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = medVSmed %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.7) + 
  theme_classic() + 
  ggtitle("Medium, 24 hr") + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigSX_genData_boxplot_withNoise",noise,"_medium_vs_medium_24_",variableDate,".pdf"),
         width = 3, height = 3)


# bootstrapping with 1000 replicates
# Cisplatin 3.3 uM vs. medium
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("TP53_cpt3","MDM2_cpt3","CDKN1A_cpt3","BTG2_cpt3"))

ci_p53 <- boot.ci(boot.out = results_med_cpt3, 
                  index = c(1,1),
                  type = c("norm"))
ci_mdm2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(2,2),
                   type = c("norm"))
ci_p21 <- boot.ci(boot.out = results_med_cpt3, 
                  index = c(3,3),
                  type = c("norm"))
ci_btg2 <- boot.ci(boot.out = results_med_cpt3, 
                   index = c(4,4),
                   type = c("norm"))

res_sum_cpt <- data.frame(name = factor(c("TP53_cpt_vs_TP53_med","MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med"),
                                        levels = c("TP53_cpt_vs_TP53_med","MDM2_med_vs_TP53_med", "CDKN1A_med_vs_TP53_med", "BTG2_med_vs_TP53_med")),
                          mean = c(ci_p53$t0,ci_mdm2$t0,ci_p21$t0,ci_btg2$t0),
                          ci_lower = c(ci_p53$normal[2],ci_mdm2$normal[2],ci_p21$normal[2],ci_btg2$normal[2]),
                          ci_upper = c(ci_p53$normal[3],ci_mdm2$normal[3],ci_p21$normal[3],ci_btg2$normal[3]))


hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)
ggplot() +
  geom_blank(data = cptVSmed24 %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_cpt,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_cpt, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  geom_boxplot(data = cptVSmed24 %>% filter(variation == 0.2), aes(x = axisLab2, y = value, fill = variation), position = "dodge2", width = 0.7) + 
  theme_classic() + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("c"))),labels = legendlabs) + 
  ggtitle("Cisplatin, 24 hr")
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigSX_genData_boxplot_withNoise",noise,"_cisplatin3.3uM_vs_medium_24_",variableDate,".pdf"),
         width = 3.5, height = 3)

######### ######### ######### ######### ######### ######### 
##### Ia. Do analysis on kd_p53(p)_mdm2, without noise #### 
######### ######### ######### ######### ######### ######### 

# Vary R parm plot

variableDate <- "20220321_varyR_rFactor_MDM2Feedback" 

files_tmp <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_w"))] 
files_tmp
parm_pert_file <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_p"))] 
parm_pert_file

varyR <- T
withNoise <- F

# Printing for information
print(paste("Working on plots for file", paste(files_tmp, sep = " ", collapse = ", ")))

genData_orig <- fread(paste0(genDataInputDir,files_tmp))
pert_orig <- fread(paste0(genDataInputDir,parm_pert_file))

# Check for negative parameter values
pert <- pert_orig %>% filter(!ParName %in% c("sf_p53","sf_mdm2","sf_p21","sf_btg2",
                                             "offset_p53","offset_mdm2","offset_p21","offset_btg2"))
cnames <- pert %>% pull(ParName)
pert <- as_tibble(t(pert %>% dplyr::select(-c(V1, ParName))))
colnames(pert) <- cnames
if (min(pert) < 0) {
  if (min(pert %>% dplyr::select(contains("_init"))) < 0) {
    warning("Oops! Negative values detected in steady states!")
    print(min(pert %>% dplyr::select(contains("_init"))))
  } 
  if (min(pert %>% dplyr::select(contains("k"))) < 0) {
    warning("Oops! Negative values detected in parameters!")
    print(min(pert %>% dplyr::select(contains("k"))))
  }
}

# Make tibble data
genData <- as_tibble(genData_orig)

# Change column names
genData <- changeCols(genData)

# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Remove outliers
print("Error!")
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

genData <- addNoise(genData, withNoise, rfactor, noise)
  
# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))


# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)

# Do correlation analysis per set
corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)


# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

# Make med and cpt bootstrap data
# Create dataframes for plotting
# Make data for plotting
l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

# Do for 8 hours
# bootstrapping with 1000 replicates
# Medium
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med8,
                   type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed8 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("kd"["p53 mdm2"]), " & ",italic("kd"["p53p mdm2"]), ", 8 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10B_kd_p53p_mdm2_genData_boxplot_varyR_8hr_",variableDate,".pdf"),
           width = 5, height = 3)



# Do for 24 hours
# bootstrapping with 1000 replicates
# Medium
results_med24 <- boot(data=df_med_24, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med24,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed24 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("kd"["p53 mdm2"]), " & ",italic("kd"["p53p mdm2"]), ", 24 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10B_kd_p53p_mdm2_genData_boxplot_varyR_24hr_",variableDate,".pdf"),
           width = 5, height = 3)

######### ######### ######### ######### ######### ######### 
##### Ib. Do analysis on kd_p53(p)_mdm2, with noise #### 
######### ######### ######### ######### ######### ######### 

withNoise <- T

# Make tibble data
genData <- as_tibble(genData_orig)
#varyR <- T

# Change column names
genData <- changeCols(genData)


# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Remove outliers
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add noise
genData <- addNoise(genData, withNoise, rfactor, noise)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)

# Do correlation analysis per set
corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)


# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

# Make med and cpt bootstrap data
# Make data for plotting
l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

# Do for 8 hours
# bootstrapping with 1000 replicates
# Medium
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med8,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed8 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("kd"["p53 mdm2"]), " & ",italic("kd"["p53p mdm2"]), ", 8 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/Fig5C_kd_p53p_mdm2_genData_boxplot_withNoise",noise,"_varyR_8hr_",variableDate,".pdf"),
         width = 5, height = 3)


# Do for 24 hours
# bootstrapping with 1000 replicates
# Medium
results_med24 <- boot(data=df_med_24, statistic=correlation,
                      R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                      rows_select = "TP53",
                      cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med24,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed24 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("kd"["p53 mdm2"]), " & ",italic("kd"["p53p mdm2"]), ", 24 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigSX_kd_p53p_mdm2_genData_boxplot_withNoise",noise,"_varyR_24hr_",variableDate,".pdf"),
         width = 5, height = 3)

######### ######### ######### ######### ######### ######### 
######## IIa. Do analysis on k_dp, without noise  ######### 
######### ######### ######### ######### ######### ######### 

# Vary R parm plot

variableDate <- "20220322_varyR_dephos" 

files_tmp <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_w"))] 
files_tmp
parm_pert_file <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_p"))] 
parm_pert_file

varyR <- T
withNoise <- F
withRNA <- T

# Printing for information
print(paste("Working on plots for file", paste(files_tmp, sep = " ", collapse = ", ")))

genData_orig <- fread(paste0(genDataInputDir,files_tmp))
pert_orig <- fread(paste0(genDataInputDir,parm_pert_file))

# Check for negative parameter values
pert <- pert_orig %>% filter(!ParName %in% c("sf_p53","sf_mdm2","sf_p21","sf_btg2",
                                             "offset_p53","offset_mdm2","offset_p21","offset_btg2"))
cnames <- pert %>% pull(ParName)
pert <- as_tibble(t(pert %>% dplyr::select(-c(V1, ParName))))
colnames(pert) <- cnames
if (min(pert) < 0) {
  if (min(pert %>% dplyr::select(contains("_init"))) < 0) {
    warning("Oops! Negative values detected in steady states!")
    print(min(pert %>% dplyr::select(contains("_init"))))
  } 
  if (min(pert %>% dplyr::select(contains("k"))) < 0) {
    warning("Oops! Negative values detected in parameters!")
    print(min(pert %>% dplyr::select(contains("k"))))
  }
}

# Make tibble data
genData <- as_tibble(genData_orig)

# Change column names
genData <- changeCols(genData)



# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Remove outliers
print("Error!")
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add noise
genData <- addNoise(genData, withNoise, rfactor, noise)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))


# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)

# Do correlation analysis per set
corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)


# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

# Make med and cpt bootstrap data
# Make data for plotting
l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

# Do for 8 hours
# bootstrapping with 1000 replicates
# Medium
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med8,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed8 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("k"["dp"]), ", 8 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10C_dephos_genData_boxplot_varyR_8hr_",variableDate,".pdf"),
       width = 5, height = 3)



# Do for 24 hours
# bootstrapping with 1000 replicates
# Medium
results_med24 <- boot(data=df_med_24, statistic=correlation,
                      R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                      rows_select = "TP53",
                      cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med24,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed24 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("k"["dp"]), ", 24 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10C_dephos_genData_boxplot_varyR_24hr_",variableDate,".pdf"),
       width = 5, height = 3)


######### ######### ######### ######### ######### ######### 
######### IIb. Do analysis on dephos, with noise  ######### 
######### ######### ######### ######### ######### ######### 

withNoise <- T

# Make tibble data
genData <- as_tibble(genData_orig)
#varyR <- T

# Change column names
genData <- changeCols(genData)


# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Remove outliers
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add noise
genData <- addNoise(genData, withNoise, rfactor, noise)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))


# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)

# Do correlation analysis per set
corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)


# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

# Make med and cpt bootstrap data
# Make data for plotting
l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

# Do for 8 hours
# bootstrapping with 1000 replicates
# Medium
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med8,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed8 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("k"["dp"]), ", 8 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/Fig5C_dephos_genData_boxplot_withNoise",noise,"_varyR_8hr_",variableDate,".pdf"),
       width = 5, height = 3)


# Do for 24 hours
# bootstrapping with 1000 replicates
# Medium
results_med24 <- boot(data=df_med_24, statistic=correlation,
                      R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                      rows_select = "TP53",
                      cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med24,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed24 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("k"["dp"]), ", 24 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigSX_dephos_genData_boxplot_withNoise",noise,"_varyR_24hr_",variableDate,".pdf"),
       width = 5, height = 3)


######### ######### ######### ######### ######### ######### 
#### IIIa. Do analysis on ks_mdm2_p53p, without noise ######## 
######### ######### ######### ######### ######### ######### 

# Vary R parm plot

variableDate <- "20220322_varyR_MDM2sythesis" 

files_tmp <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_w"))] 
files_tmp
parm_pert_file <- genDataFiles[sapply(genDataFiles, grepl, pattern = paste0(variableDate,"_",modelDate,"_",MODEL,"_1000times50SampleSimulations_p"))] 
parm_pert_file

varyR <- T
withNoise <- F
withRNA <- T

# Printing for information
print(paste("Working on plots for file", paste(files_tmp, sep = " ", collapse = ", ")))

genData_orig <- fread(paste0(genDataInputDir,files_tmp))
pert_orig <- fread(paste0(genDataInputDir,parm_pert_file))

# Check for negative parameter values
pert <- pert_orig %>% filter(!ParName %in% c("sf_p53","sf_mdm2","sf_p21","sf_btg2",
                                             "offset_p53","offset_mdm2","offset_p21","offset_btg2"))
cnames <- pert %>% pull(ParName)
pert <- as_tibble(t(pert %>% dplyr::select(-c(V1, ParName))))
colnames(pert) <- cnames
if (min(pert) < 0) {
  if (min(pert %>% dplyr::select(contains("_init"))) < 0) {
    warning("Oops! Negative values detected in steady states!")
    print(min(pert %>% dplyr::select(contains("_init"))))
  } 
  if (min(pert %>% dplyr::select(contains("k"))) < 0) {
    warning("Oops! Negative values detected in parameters!")
    print(min(pert %>% dplyr::select(contains("k"))))
  }
}

# Make tibble data
genData <- as_tibble(genData_orig)

# Change column names
genData <- changeCols(genData)

# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Remove outliers
print("Error!")
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add noise
genData <- addNoise(genData, withNoise, rfactor, noise)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))


# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)

# Do correlation analysis per set
corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)


# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

# Make med and cpt bootstrap data
# Make data for plotting
l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

# Do for 8 hours
# bootstrapping with 1000 replicates
# Medium
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med8,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed8 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("ks"["mdm2 p53p"]), ", 8 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10D_ks_mdm2_p53p_genData_boxplot_varyR_8hr_",variableDate,".pdf"),
       width = 5, height = 3)

# Do for 24 hours
# bootstrapping with 1000 replicates
# Medium
results_med24 <- boot(data=df_med_24, statistic=correlation,
                      R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                      rows_select = "TP53",
                      cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med24,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed24 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("ks"["mdm2 p53p"]), ", 24 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigS10D_ks_mdm2_p53p_genData_boxplot_varyR_24hr_",variableDate,".pdf"),
       width = 5, height = 3)



######### ######### ######### ######### ######### ######### 
#### IIIb. Do analysis on ks_mdm2_p53p, with noise ######## 
######### ######### ######### ######### ######### ######### 

withNoise <- T

# Make tibble data
genData <- as_tibble(genData_orig)

# Change column names
genData <- changeCols(genData)


# Split information
genData <- genData %>% mutate(run = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",1),"r")),
                              sample = sapply(strsplit(V1,"_"),"[",2),
                              variation = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",3),"v")))
if (varyR) {
  genData <- genData %>% mutate(rfactor = as.factor(str_remove(sapply(strsplit(V1,"_"),"[",4),"rfactor")))
}
genData <- dplyr::select(genData, -V1)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Remove outliers
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))

# Add noise
genData <- addNoise(genData, withNoise, rfactor, noise)

# Check for negative values
minGenData <- genData %>% dplyr::select(-c(sample,variation, run, rfactor)) %>% summarise_all(min)
if (min(minGenData) < 0) {
  warning("Oops! Negative values detected!")
}

# Filter out outliers
# Remove NAs
genData <- genData %>% filter_at(vars(matches(c("med","cpt"))), any_vars(!is.na(.)))
# Remove outliers
genData <- genData %>% dplyr::select(-contains("noise")) %>% filter_at(vars(matches(c("med","cpt"))), all_vars(. > 0))


# Add p53 total
genData <- genData %>% mutate("p53_total_med" = p53_med + p53P_med,
                              "p53_total_cpt8" = p53_cpt8 + p53P_cpt8,
                              "p53_total_cpt24" = p53_cpt24 + p53P_cpt24)

# Do correlation analysis per set
corGenData <- makeCorGenData(genData, withRNA, varyR, mdm2RNA)


# Pivot data frame
corGenData <- pivot_longer(corGenData, cols = contains("vs"))
corGenData <- corGenData %>% mutate(axisLab = str_replace_all(name, "_"," "),
                                    axisLab = str_replace_all(axisLab, "vs ","vs.\n"),
                                    axisLab2 = sapply(name, function(n){str_split(n,"_")[[1]][1]}))

# Make med and cpt bootstrap data
# Make data for plotting
l <- plottingDataFrames(corGenData, mdm2RNA, withRNA)
medVSmed <- l[["medVSmed"]]
cptVSmed8 <- l[["cptVSmed8"]]
cptVScpt8 <- l[["cptVScpt8"]]
cptVSmed24 <- l[["cptVSmed24"]]
cptVScpt24 <- l[["cptVScpt24"]]

# Do for 8 hours
# bootstrapping with 1000 replicates
# Medium
results_med8 <- boot(data=df_med_8, statistic=correlation,
                     R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                     rows_select = "TP53",
                     cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med8,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_8)), rho=medVSmedReal_8)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_8,
                         df_cpt3_8,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_8)), rho=cpt3VSmedReal_8) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed8 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("ks"["mdm2 p53p"]), ", 8 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/Fig5C_ks_mdm2_p53p_genData_boxplot_withNoise",noise,"_varyR_8hr_",variableDate,".pdf"),
       width = 5, height = 3)


# Do for 24 hours
# bootstrapping with 1000 replicates
# Medium
results_med24 <- boot(data=df_med_24, statistic=correlation,
                      R=1000, probes = PROTS_OF_INTEREST, decims = 2,
                      rows_select = "TP53",
                      cols_select = "MDM2")
ci_mdm2_med <- boot.ci(boot.out = results_med24,
                       type = c("norm"))
res_sum_med <- data.frame(name = factor(c("MDM2_med_vs_TP53_med"),
                                        levels = c("MDM2_med_vs_TP53_med")),
                          mean = c(ci_mdm2_med$t0),
                          ci_lower = c(ci_mdm2_med$normal[2]),
                          ci_upper = c(ci_mdm2_med$normal[3]))
hline_med <- data.frame(axisLab2=seq(1,length(medVSmedReal_24)), rho=medVSmedReal_24)
res_sum_med$axisLab2 <- as.numeric(res_sum_med$name)
# Cisplatin
merged_med_cpt3 <- merge(df_med_24,
                         df_cpt3_24,
                         by = "CELL_ID",
                         suffixes = c("_med","_cpt3"))
results_med_cpt3 <- boot(data=merged_med_cpt3, statistic=correlation,
                         R=1000, probes = c(paste0(PROTS_OF_INTEREST,"_med"),paste0(PROTS_OF_INTEREST,"_cpt3")), 
                         decims = 3,
                         rows_select = "TP53_med",
                         cols_select = c("MDM2_cpt3"))
ci_mdm2_cpt <- boot.ci(boot.out = results_med_cpt3,
                       type = c("norm"))
res_sum_cpt <- data.frame(name = factor(c("MDM2_cpt_vs_TP53_med"),
                                        levels = c("MDM2_cpt_vs_TP53_med")),
                          mean = c(ci_mdm2_cpt$t0),
                          ci_lower = c(ci_mdm2_cpt$normal[2]),
                          ci_upper = c(ci_mdm2_cpt$normal[3]))
hline_cpt <- data.frame(axisLab2=seq(1,length(cpt3VSmedReal_24)), rho=cpt3VSmedReal_24) 
res_sum_cpt$axisLab2 <- as.numeric(res_sum_cpt$name)

res_sum_rdata <- bind_rows(res_sum_med %>% filter(name == "MDM2_med_vs_TP53_med") %>% mutate(axisLab2 = 1),
                           res_sum_cpt %>% filter(name == "MDM2_cpt_vs_TP53_med") %>% mutate(axisLab2 = 2))

hline_rdata <- bind_rows(hline_med[1,],
                         hline_cpt[2,])
# Make r factor data
rdata <- bind_rows(medVSmed %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "medium"),
                   cptVSmed24 %>% filter(axisLab2 == "MDM2", variation == 0.2) %>% mutate(comparison = "cisplatin")) %>%
  mutate(comparison = factor(comparison, levels = c("medium","cisplatin")))

ggplot() +
  geom_blank(data = rdata , aes(x = comparison, y = value, fill = rfactor)) +
  geom_hline(aes(yintercept = 0)) +
  geom_rect(data = res_sum_rdata,
            aes(xmin = axisLab2-0.45,
                xmax = axisLab2+0.45,
                ymin = ci_lower,
                ymax = ci_upper,
                group = axisLab2),
            fill = "grey", alpha = 0.5) +
  geom_segment(data = res_sum_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=mean, yend=mean), color = "grey70", linetype = "solid") +
  #  geom_segment(data = hline_rdata, aes(x=axisLab2-0.45, xend=axisLab2+0.45, y=rho, yend=rho), color = "grey40", linetype = "dashed") +
  geom_boxplot(data = rdata , aes(x = comparison, y = value, fill = rfactor), position = "dodge2", width = 0.8) + 
  theme_classic() + 
  ggtitle(expression(paste(italic("ks"["mdm2 p53p"]), ", 24 hr"))) + 
  ylab("Pearson correlation") + ylim(c(-1,1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(discrete = T, name = expression(paste("factor ",italic("r"))),labels = legendlabs)
ggsave(file = paste0(outputDir_GR,DATE,"_Correlation_Figs/FigSX_ks_mdm2_p53p_genData_boxplot_withNoise",noise,"_varyR_24hr_",variableDate,".pdf"),
       width = 5, height = 3)

