# Description:
# R script to make data for modeling purpose. 
# Filtering of data.
# Calculate mean of technical replicates.
# 1) p53, Mdm2, Btg2 and p21 reporters
# 2) up to 4 replicates, and 
# 3) GFP channels and AnV/PI channels

# Last modified: 2022-04-13
# Written by Muriel Heldring

# R version 4.1.3 (2022-03-10) -- "One Push-Up"

# Clear environment
rm(list = ls())
gc()

# Load the packages
library(tidyverse)
library(rlang)
library(lemon)
library(data.table)
library(splines)

# Choose project
PROJECT <- "DDP" 

# Load the global variables
load(file = paste0("/data/muriel/Projects/",PROJECT,"/DataAnalysis/Dynamics/Output/RData/00_constants.RData"))
# Load the functions
source(paste0("/data/muriel/Projects/","DDP","/DataAnalysis/Dynamics/","RScripts_ImagingDataAnalysis/01_Functions.R"))

# Choose project
PROJECT_PATH <- DDP_FOLDER_PATH 

DATE <- "20220411" 
outputDir_GR <- "/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/Figures/"

REPLICATE_ORDER <- c("1","2","3","4")
PROTEIN_ORDER <- c("p53","MDM2","p21","BTG2")

dir.create(paste0(outputDir_GR,DATE,"_GFPdata_Figs"))

# -----------------------------------------------------------------------------
# ----------------------------MAKE SUMMARY DATA--------------------------------
# -----------------------------------------------------------------------------

# Save the summary data as it is
# load(file = paste0(PROJECT_PATH, OUTPUT_PATH,"RData/04_", "summaryDataPerImage.RData"))  # Load the data frame
# write_csv(summaryDataPerImage, file = "/data/muriel/Projects/PHH/DataAnalysis/GFPdata/SummaryData.csv")
summaryDataPerImage <- read_csv(file = "/data/muriel/Projects/PHH/DataAnalysis/GFPdata/SummaryData.csv")

# -------------------MAKE PLOTS OF THE TECHNICAL REPLICATES PER WELL-------------------
dirPath <- paste0(PROJECT_PATH,OUTPUT_PATH,"Figures/",DATE,"_TechnicalReplicates/")
dir.create(file.path(dirPath))

for (p in PROTS_OF_INTEREST) {
  for (t in COMPS_OF_INTEREST) {
    dftmp <- summaryDataPerImage %>% filter(protein == p & 
                                              treatment == t) 
    for (r in c(1,2,3,4)) {
      if (r %in% levels(droplevels(dftmp$replID))) {
        dftmp2 <- dftmp %>% filter(replID == r) %>% select(timeID,well,imageNumber, GM10Nuclei_Integrated_Intensity, GM10Cytoplasm_Integrated_Intensity,imageCountParentObj)
        dftmp_piv <- dftmp2 %>% pivot_wider(names_from = imageNumber, values_from = c(GM10Nuclei_Integrated_Intensity,GM10Cytoplasm_Integrated_Intensity,imageCountParentObj))
        
        if (p %in% c("p53","p21","MDM2")) {
          gplottmp <- ggplot(dftmp_piv) + geom_point(aes_string(x = paste("GM10Nuclei_Integrated_Intensity","1",sep = "_"),
                                                                y = paste("GM10Nuclei_Integrated_Intensity","2",sep = "_"), color = "well")) +
            xlab("GFP in technical repliclate 1") + ylab("GFP in technical repliclate 2") 
          ggsave(paste0(dirPath, "NucleiIntensity_",p,"_",t,"_rep",r,".pdf"), plot = gplottmp,width = 7, height = 5)}
        else {
          gplottmp <- ggplot(dftmp_piv) + geom_point(aes_string(x = paste("GM10Cytoplasm_Integrated_Intensity","1",sep = "_"),
                                                                y = paste("GM10Cytoplasm_Integrated_Intensity","2",sep = "_"), color = "well")) +
            xlab("GFP in technical repliclate 1") + ylab("GFP in technical repliclate 2") 
          ggsave(paste0(dirPath, "CytoplamsIntensity_",p,"_",t,"_rep",r,".pdf"), plot = gplottmp,width = 7, height = 5)}
        
        gplottmp <- ggplot(dftmp_piv) + geom_point(aes_string(x = paste("imageCountParentObj","1",sep = "_"),
                                                              y = paste("imageCountParentObj","2",sep = "_"), color = "well")) +
          xlab("Cell count in technical repliclate 1") + ylab("Cell count in technical repliclate 2") 
        ggsave(paste0(dirPath, "CellCount",p,"_",t,"_rep",r,".pdf"), plot = gplottmp,width = 7, height = 5)
      }}}}

# -------------------REMOVE TECHNICAL REPLICATES-------------------
summaryDataPerImageFiltered <- summaryDataPerImage %>% 
  filter(!((protein == "p21" & treatment == "DMEM" & replID == 4 & locationID == "F23_1") | 
             (protein == "p53" & treatment == "ETO" & replID == 3 & locationID == "H19_2") |
             (protein %in% c("MDM2","p21","BTG2") & replID == 1 & treatment == "MYT")))

# -------------------CALCULATE AVERAGES PER WELL-------------------
summaryData <- summaryDataPerImageFiltered %>% group_by(treatment,dose_uM,cell_line,plateID,protein,
                                                        replID,timeID,timeAfterExposure,well) %>% 
  summarise(imageCountParentObj = mean(imageCountParentObj, na.rm = T),
            Cytoplasm_Mean_Intensity = ifelse(all(is.na(GM10Cytoplasm_Mean_Intensity)),NA,mean(GM10Cytoplasm_Mean_Intensity, na.rm = T)),
            Cytoplasm_Integrated_Intensity = ifelse(all(is.na(GM10Cytoplasm_Integrated_Intensity)),NA,mean(GM10Cytoplasm_Integrated_Intensity, na.rm = T)),
            Nuclei_Mean_Intensity = ifelse(all(is.na(GM10Nuclei_Mean_Intensity)),NA,mean(GM10Nuclei_Mean_Intensity, na.rm = T)),
            Nuclei_Integrated_Intensity = ifelse(all(is.na(GM10Nuclei_Integrated_Intensity)),NA,mean(GM10Nuclei_Integrated_Intensity, na.rm = T)),
            Nuclei_AreaShape_Area = mean(Nuclei_AreaShape_Area, na.rm = T),
            propNonRespCytoplasm_Mean_Intensity = ifelse(all(is.na(propNonRespCytoplasm_Mean_Intensity)),NA,mean(propNonRespCytoplasm_Mean_Intensity, na.rm = T)),
            propNonRespCytoplasm_Integrated_Intensity = ifelse(all(is.na(propNonRespCytoplasm_Integrated_Intensity)),NA,mean(propNonRespCytoplasm_Integrated_Intensity, na.rm = T)),
            propNonRespNuclei_Mean_Intensity = ifelse(all(is.na(propNonRespNuclei_Mean_Intensity)),NA,mean(propNonRespNuclei_Mean_Intensity, na.rm = T)),
            propNonRespNuclei_Integrated_Intensity = ifelse(all(is.na(propNonRespNuclei_Integrated_Intensity)),NA,mean(propNonRespNuclei_Integrated_Intensity, na.rm = T)),
            cellCount_AnVPI = mean(cellCount_AnVPI, na.rm = T),
            propCount_PI = mean(propCount_PI, na.rm = T),
            propCount_AnV = mean(propCount_AnV, na.rm = T))

# -------------------ADD RELATIVE CELL COUNTS TO SUMMARYDATA-------------------
summaryData <- summaryData %>% group_by(replID,treatment,protein,dose_uM) %>% 
  dplyr::mutate(relativeCellCount = imageCountParentObj/imageCountParentObj[1])

# -------------------RELEVEL THE REPLICATES AND PROTEINS-----------------------
summaryData$replID <- factor(summaryData$replID, levels = REPLICATE_ORDER)
summaryData$protein <- factor(summaryData$protein, levels = PROTEIN_ORDER)

# In project DDP; set the AnV and PI zero's in p21 replicate 3 at timeID >= 36 to NA
summaryData[summaryData$protein == "p21" & 
              summaryData$replID == 3 &
              summaryData$timeID > 35,c("propCount_AnV","propCount_PI")] <-  NA

# Refactorize the dose_um and dose_uMadj columns
summaryData$dose_uM <- factor(as.numeric(as.character(summaryData$dose_uM)))

# Summarise controls

dfTmp <- summaryData %>% filter(treatment %in% c("CDDP","DMEM")) %>% 
  mutate(dose_uMadj = factor(ifelse(as.character(dose_uM) == "100", "0", as.character(dose_uM)),levels = c(0,1,2.5,5,10,15,20,25,50)),
         dataForModel = ifelse(protein %in% c("p53", "p21","MDM2"),
                               Nuclei_Integrated_Intensity,
                               Cytoplasm_Integrated_Intensity)) %>%
  group_by(treatment, dose_uMadj, cell_line, plateID, protein, replID, timeID, timeAfterExposure) %>%
  summarise(dataForModel = ifelse(treatment == "DMEM",
                                  mean(dataForModel, na.rm = T),
                                  dataForModel)) %>%
  ungroup()

repnames <- c("1" = "Repl 1",
              "2" = "Repl 2",
              "3" = "Repl 3",
              "4" = "Repl 4",
              "p53" = "p53",
              "MDM2" = "MDM2",
              "p21" = "p21",
              "BTG2" = "BTG2")

# Make a plot
ggplot() + geom_point(data = dfTmp, aes(x = timeAfterExposure, y = dataForModel, color = dose_uMadj), size = 0.7) + 
  theme_classic() +
  scale_color_viridis_d(name = expression(paste("Cisplatin ("*mu*"M)"))) +
  ylab("GFP intensity (a.u.)") +
  xlab("Time (h)") +
  facet_grid(protein~replID, labeller = as_labeller(repnames)) +
  theme(axis.title = element_text(size = 18),
        strip.text = element_text(size = 14))
ggsave(filename = paste0(outputDir_GR,DATE,"_GFPdata_Figs/","FigS6A_CDDP_data_unnormalised.pdf"), width = 6, height = 6)

# -----------------------------------------------------------------------------
# -----------------------SUBTRACT CONTROL CONDITION PER PLATE------------------
# -----------------------------------------------------------------------------

df <- summaryData
refTreatment <- c("DMEM","DMSO")
refConcentration <- c(100,0.2)

# Get control summary
dfControl <- df[df$treatment %in% refTreatment & df$dose_uM %in% refConcentration,] %>% 
  select(-c(well)) %>%
  group_by(treatment, dose_uM, cell_line, plateID, protein, replID, timeID, timeAfterExposure) %>% 
  summarise_all(list(mean), na.rm = T) %>% ungroup()

# Make summary per control condition
dfRef <- df[df$treatment %in% refTreatment & df$dose_uM %in% refConcentration,] %>% 
  select(-c(well,imageCountParentObj,
            Nuclei_AreaShape_Area,
            propNonRespCytoplasm_Mean_Intensity,
            propNonRespCytoplasm_Integrated_Intensity,
            propNonRespNuclei_Mean_Intensity,
            propNonRespNuclei_Integrated_Intensity,
            cellCount_AnVPI,
            propCount_PI,
            propCount_AnV,
            relativeCellCount)) %>%
  group_by(treatment, dose_uM, cell_line, plateID, protein, replID, timeID, timeAfterExposure) %>% 
  summarise_all(list(control = mean), na.rm = T) %>% ungroup()

dfRefDMEM <- dfRef %>% filter(treatment == "DMEM") %>% select(-c(treatment,dose_uM))
dfRefDMSO <- dfRef %>% filter(treatment == "DMSO") %>% select(-c(treatment,dose_uM))

# Check on NAs
if (!any(is.na(dfRef))) {print("No NAs found!")}

# Select treated condition
dfCDDP <- df %>% filter(treatment == "CDDP")
dfDMEM <- dfControl %>% filter(treatment == "DMEM")

# Merge dfTreatments with dfControl
dfMergeCDDP <- left_join(dfCDDP,dfRefDMEM, by = c("cell_line", "plateID", "protein", "replID", "timeID", "timeAfterExposure"))
dfMergeDMEM <- left_join(dfDMEM,dfRefDMEM, by = c("cell_line", "plateID", "protein", "replID", "timeID", "timeAfterExposure"))

# Bind the data frames
dfMerge <- bind_rows(dfMergeCDDP,dfMergeDMEM)

# Convert NaN to NA
dfMergeCDDP <- dfMergeCDDP %>% mutate_if(is.numeric, list(~na_if(., NaN)))

# Do background subtraction per protein
summaryDataBackgroundSub <- dfMerge %>% mutate(
  Cytoplasm_Mean_Intensity = Cytoplasm_Mean_Intensity - Cytoplasm_Mean_Intensity_control,
  Cytoplasm_Integrated_Intensity = Cytoplasm_Integrated_Intensity - Cytoplasm_Integrated_Intensity_control,
  Nuclei_Mean_Intensity = Nuclei_Mean_Intensity - Nuclei_Mean_Intensity_control,
  Nuclei_Integrated_Intensity = Nuclei_Integrated_Intensity - Nuclei_Integrated_Intensity_control)

# Make a plot
dfTmp <- summaryDataBackgroundSub %>% filter(treatment %in% c("CDDP","DMEM"), dose_uM %in% c(100,1,2.5,5))
ggplot() + geom_point(data = dfTmp, aes(x = timeID, y = Nuclei_Integrated_Intensity, color = dose_uM)) + 
  facet_grid(protein~replID)

# -----------------------------------------------------------------------------
# ---------------------------DO MIN-MAX NORMALISATION--------------------------
# -----------------------------------------------------------------------------

# Do min-max normalisation for all conditions
summaryDataNorm <- doMinMaxNormAllConditions(
  df = summaryDataBackgroundSub,
  normPerCompCtrlSet = F,
  columnsForNormalisation = list(c("Cytoplasm_Mean_Intensity","Nuclei_Mean_Intensity"),
                                 c("Cytoplasm_Integrated_Intensity", "Nuclei_Integrated_Intensity")),
  by = c("protein","replID"), combineCytAndNuc = T)

# Make a plot
colnames(summaryDataNorm)
dfTmp <- summaryDataNorm %>% filter(treatment %in% c("CDDP","DMEM"), dose_uM %in% c(100,1,2.5,5))
ggplot() + geom_point(data = dfTmp, aes(x = timeID, y = Nuclei_Integrated_Intensity_Norm, color = dose_uM)) + 
  facet_grid(protein~replID)

# Replace the dose_uM value for the control concentration with 0
summaryDataNorm$dose_uM <- as.numeric(as.character(droplevels(summaryDataNorm$dose_uM)))
summaryDataNorm$dose_uMadj <- sapply(1:nrow(summaryDataNorm), function(rownr){
  ifelse(summaryDataNorm$dose_uM[rownr] == 100 & 
           summaryDataNorm$treatment[rownr] == "DMEM", 0, summaryDataNorm$dose_uM[rownr])})

# Refactorize the dose_uM and dose_uMadj columns
summaryDataNorm$dose_uMadj <- factor(as.numeric(as.character(summaryDataNorm$dose_uMadj)))
summaryDataNorm$dose_uM <- factor(as.numeric(as.character(summaryDataNorm$dose_uM)))

# Interpolate the separate replicates
summaryDataNorm <- as.data.frame(summaryDataNorm %>% group_by(treatment, dose_uM, replID, protein) %>% mutate(
  timepoints = seq(1,TIME_BETWEEN_FRAMES*length(timeID),TIME_BETWEEN_FRAMES),
  interpol_Nuclei_Integrated_Intensity_Norm = predict(lm(Nuclei_Integrated_Intensity_Norm ~ bs(timeAfterExposure,df = 6, degree = 3,intercept = F)), 
                                                      data.frame(timeAfterExposure = timepoints)),
  interpol_Nuclei_Mean_Intensity_Norm = predict(lm(Nuclei_Mean_Intensity_Norm ~ bs(timeAfterExposure,df = 6, degree = 3,intercept = F)), 
                                                data.frame(timeAfterExposure = timepoints)),
  interpol_Cytoplasm_Integrated_Intensity_Norm = predict(lm(Cytoplasm_Integrated_Intensity_Norm ~ bs(timeAfterExposure,df = 6, degree = 3,intercept = F)), 
                                                         data.frame(timeAfterExposure = timepoints)),
  interpol_Cytoplasm_Mean_Intensity_Norm = predict(lm(Cytoplasm_Mean_Intensity_Norm ~ bs(timeAfterExposure,df = 6, degree = 3,intercept = F)), 
                                                   data.frame(timeAfterExposure = timepoints)),
  interpol_relativeCellCount = predict(lm(relativeCellCount ~ bs(timeAfterExposure,df = 6, degree = 3,intercept = F)), 
                                       data.frame(timeAfterExposure = timepoints)),
  interpol_absoluteCellCount = predict(lm(imageCountParentObj ~ bs(timeAfterExposure,df = 6, degree = 3,intercept = F)), 
                                       data.frame(timeAfterExposure = timepoints))))

# Make a plot
colnames(summaryDataNorm)
dfTmp <- summaryDataNorm %>% filter(treatment %in% c("CDDP","DMEM"), dose_uM %in% c(100,1,2.5,5))
ggplot() + geom_point(data = dfTmp, aes(x = timepoints, y = Nuclei_Integrated_Intensity_Norm, color = dose_uM)) + 
  facet_grid(protein~replID)

dfTmp <- summaryDataNorm %>% filter(treatment %in% c("CDDP","DMEM"), dose_uM %in% c(100,1,2.5,5))
ggplot() + geom_point(data = dfTmp, aes(x = timepoints, y = interpol_Nuclei_Integrated_Intensity_Norm, color = dose_uM)) + 
  facet_grid(protein~replID)

# Make data for parameter estimation.
# Make one data column to be used for parameter estimation, i.e. the GFP expression 
# of data in nucleus and cytoplasm according to their location of expression 
summaryDataNorm$data4modelInterpol <- sapply(1:nrow(summaryDataNorm), function(rowNr){
  ifelse(summaryDataNorm[rowNr,"protein"] %in% c("p53", "p21","MDM2"),
         summaryDataNorm[rowNr,"interpol_Nuclei_Integrated_Intensity_Norm"], 
         summaryDataNorm[rowNr,"interpol_Cytoplasm_Integrated_Intensity_Norm"])})
summaryDataNorm$data4modelReal <- sapply(1:nrow(summaryDataNorm), function(rowNr){
  ifelse(summaryDataNorm[rowNr,"protein"] %in% c("p53", "p21","MDM2"),
         summaryDataNorm[rowNr,"Nuclei_Integrated_Intensity_Norm"], 
         summaryDataNorm[rowNr,"Cytoplasm_Integrated_Intensity_Norm"])})

ggplot(summaryDataNorm %>% filter(protein %in% c("p53","p21","MDM2","BTG2"))) + 
  geom_point(aes(x = timeAfterExposure, y = data4modelReal, color = dose_uMadj)) + 
  geom_line(aes(x = timeAfterExposure, y = data4modelInterpol, color = dose_uMadj)) + 
  facet_grid(~protein)


# Split the summary data frame in seperate data frames for different compounds.
CDDPdata <- summaryDataNorm[summaryDataNorm$treatment %in% c("CDDP","DMEM"),]
DMEMdata <- summaryDataNorm[summaryDataNorm$treatment == "DMEM",]

# Make data frames including cell counts
CDDP4Model <- prepareD4M(CDDPdata, absolute = T)
DMEM4Model <- prepareD4M(DMEMdata, absolute = T)

# Save the data frames to a csv file
save(CDDP4Model, file = paste0(PROJECT_PATH, OUTPUT_PATH,"RData/08a_CDDPdata_",PROJECT,".RData"))
save(DMEM4Model, file = paste0(PROJECT_PATH, OUTPUT_PATH,"RData/08d_DMEMdata_",PROJECT,".RData"))

# # Write the data frames to a csv file
# write_csv(CDDP4Model, path = paste0(PROJECT_PATH, OUTPUT_PATH,"DataTables/",DATE,"_MH_CDDPdata_","All_",PROJECT,".csv"))
# write_csv(DMEM4Model, path = paste0(PROJECT_PATH, OUTPUT_PATH,"DataTables/",DATE,"_MH_DMEMdata_","All_",PROJECT,".csv"))


# Create Summary Figures (delay and peak)

# Change to tibble
CDDP4Model <- as_tibble(CDDP4Model)

# Change column name protein
CDDP4Model <- dplyr::rename(CDDP4Model, c(StateVar = protein))
# Add a column
CDDP4Model <- mutate(CDDP4Model, timeID_orig = timeID)

# Make a plot
CDDP4Model <- CDDP4Model %>% filter(StateVar %in% c("p53","MDM2","BTG2","p21") & dose_uMadj %in% c(0,1,2.5,5))
CDDP4Model <- CDDP4Model %>% mutate(StateVar = factor(StateVar,levels = c("p53","MDM2","p21","BTG2")))

ggplot() + 
  geom_line(data = CDDP4Model, aes(x = timepoints, y = data4modelInterpol, color = dose_uMadj)) + 
  geom_point(data = CDDP4Model, aes(x = timepoints, y = data4modelReal, color = dose_uMadj), size = 0.5) + 
  facet_grid(StateVar~replID)

# Get the maximum expression timepoint
nTP <- CDDP4Model %>% 
  group_by(timepoints, StateVar,dose_uMadj) %>% 
  summarise(nTP = n())
mergeCDDP <- left_join(CDDP4Model,nTP, by = c("dose_uMadj","StateVar","timepoints"))

# Do for everything
meanCDDP <- mergeCDDP %>% filter(!StateVar == "N") %>%
  filter(nTP >= 3) %>%
  group_by(timepoints, StateVar,dose_uMadj) %>%
  summarise(data4modelInterpolSD = sd(data4modelInterpol),
            data4modelInterpol = mean(data4modelInterpol),
            data4modelReal = mean(data4modelReal)) %>% ungroup()
sdCONTROL <- meanCDDP %>% group_by(StateVar) %>% filter(dose_uMadj == 0) %>%
  mutate(med3sd = data4modelInterpol * 1.5) %>% 
  dplyr::select(-c(data4modelInterpol,data4modelInterpolSD,data4modelReal,dose_uMadj))
meanCDDP <- left_join(meanCDDP,sdCONTROL, by = c("timepoints","StateVar"))
ggplot(meanCDDP %>% filter(StateVar %in% c("p53","p21","MDM2","BTG2"))) + 
  geom_point(aes(x = timepoints, y = data4modelReal, color = dose_uMadj)) + 
  geom_line(aes(x = timepoints, y = data4modelInterpol, color = dose_uMadj)) + 
  geom_line(aes(x = timepoints, y = med3sd), color = 'black') + 
  facet_grid(~StateVar)

dynchar <- meanCDDP %>% filter(!dose_uMadj == 0) %>% group_by(dose_uMadj, StateVar) %>% 
  summarise(maxI = max(data4modelInterpol),
            maxR = max(data4modelReal),
            tpod = timepoints[which(data4modelInterpol > med3sd)[1]],
            peakI = timepoints[which(data4modelInterpol == max(data4modelInterpol))],
            peakR = timepoints[which(data4modelReal == max(data4modelReal))])
dynchar_diff <- left_join(dynchar %>% filter(!StateVar == "p53"),
                          dynchar %>% filter(StateVar == "p53") %>% dplyr::select(-c(StateVar)),
                          by = c("dose_uMadj"), suffix = c("_dst","_p53"))
dynchar_diff <- dynchar_diff %>% mutate(diff_tpod = tpod_dst - tpod_p53,
                                        diff_peak = peakI_dst - peakI_p53)
dynchar_mean <- dynchar %>% group_by(StateVar) %>%
  summarise(mean_tpod = mean(tpod),
            mean_peak = mean(peakI))
dynchar_diff_mean <- dynchar_diff %>% group_by(StateVar) %>%
  summarise(mean_tpod = mean(diff_tpod),
            mean_peak = mean(diff_peak))

# Do for replicates separately
sdCONTROL_perRep <- mergeCDDP %>% group_by(StateVar,replID) %>% filter(dose_uMadj == 0) %>%
  mutate(med3sd = data4modelInterpol * 1.5) %>%
  dplyr::select(StateVar,replID,med3sd, timepoints)
downstrCDDP <- left_join(mergeCDDP,sdCONTROL_perRep, by = c("timepoints","StateVar","replID"))
dynchar_downstr <- downstrCDDP %>% filter(nTP >= 3,
                                          !dose_uMadj == 0) %>%
  group_by(dose_uMadj, StateVar, replID) %>% 
  summarise(maxI = max(data4modelInterpol),
            maxR = max(data4modelReal),
            tpod = timepoints[which(data4modelInterpol > med3sd)[1]],
            peakI = timepoints[which(data4modelInterpol == max(data4modelInterpol))],
            peakR = timepoints[which(data4modelReal == max(data4modelReal))])
dynchar_p53 <- dynchar %>% filter(StateVar == "p53")
dynchar_diff_perRep <- left_join(dynchar_downstr,
                                 dynchar_p53 %>% dplyr::select(-c(StateVar)),
                                 by = c("dose_uMadj"), suffix = c("_dst","_p53"))
dynchar_diff_perRep <- dynchar_diff_perRep %>% mutate(diff_tpod = tpod_dst - tpod_p53,
                                                      diff_peak = peakI_dst - peakI_p53)
# Option 1: do for replicates separately
ggplot()+
  geom_beeswarm(data = dynchar_diff_perRep %>% filter(!StateVar == "p53"), aes(x = StateVar, y = diff_peak, color = dose_uMadj), size = 2.5, cex = 5) +
  theme_classic() + 
  scale_color_viridis(discrete = T, name = expression("Conc. ("*mu*"M)")) +
  ylab("Peak delay with respect to\nmean p53 peak (h)") + xlab("")
ggsave(paste0(outputDir_GR,DATE,"_GFPdata_Figs/Fig3D_PeakDelay.pdf"),width = 4, height = 2.5)

ggplot() +
  geom_beeswarm(data = dynchar_downstr, aes(x = dose_uMadj, y = tpod, color = StateVar), size = 2.5, cex = 5) +
  theme_classic() + scale_color_viridis(discrete = T, name = "Protein") +
  ylab("Initial response latency (h)") + xlab(expression("Concentration ("*mu*"M)")) +
  ylim(0,45)
ggsave(paste0(outputDir_GR,DATE,"_GFPdata_Figs/Fig3E_InitialResponseLatency.pdf"),width = 4, height = 2.5)