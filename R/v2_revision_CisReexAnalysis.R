# Description:
# R script to analyse data of cisplatin reexposure experiment. 

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
library(readr)
library(viridis)

# library(sqldf)
# library(gridExtra)
# library(grid)
# library(prodlim)
# library(mixtools)
# library(ggpubr)

# -----------------------------------------------------------------------------
# --------------------------LOAD GLOBAL VARIABLES------------------------------
# -----------------------------------------------------------------------------

# Load the functions
source(paste0("/data/muriel/Projects/","DDP","/DataAnalysis/Dynamics/RScripts_ImagingDataAnalysis/01_Functions.R"))

INPUT_PATH <- "/data/muriel/Projects/DDP/DataAnalysis/Dynamics/Exp011_HepG2_CisplatinReexposure/H5CP_output/"
outputDir_GR <- "/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/Figures/"
DATE <- "20220411"

FILES <- dir(INPUT_PATH)
FILES

getWell <- function(locid) {
  wells <- sapply(locid,function(id){
    l <- strsplit(id,"_")
    l[[1]][1]})
}
OUTPUT_FIGPATH <- paste0(outputDir_GR,DATE,"_CisReex_Figs/")

dir.create(OUTPUT_FIGPATH)

# -----------------------------------------------------------------------------
# ----------------------------DO DATA ANALYSIS---------------------------------
# -----------------------------------------------------------------------------

# # Data loading
# data_p53 <- read_tsv(paste0(INPUT_PATH,FILES[2]))
# 
# # Make factors of time, dose_uM and replID
# data_p53 <- data_p53 %>% mutate(dose_uM = as.factor(dose_uM),
#                                 replID = as.factor(replID),
#                                 timeID = as.factor(timeID))
# 
# # Filter out unrelevant columns
# data_p53 <- data_p53 %>% select(-c(imageNumber,groupNumber,groupInd,imageID,Nuclei_AreaShape_Area,
#                                    groupNumber,plateWellID,control,Nuclei_Number_Object_Number))
# 
# # Calculate summaries per technical condition, i.e. make population data
# data_p53_ps <- data_p53 %>% group_by(treatment, dose_uM, replID, plateID, cell_line, 
#                                      timeID, timeAfterExposure, locationID) %>% 
#   summarise(IntegratedIntensity_GM = geomMean10(Nuclei_Intensity_IntegratedIntensity_image_gfp),
#             MeanIntensity_GM = geomMean10(Nuclei_Intensity_MeanIntensity_image_gfp),
#             Cell_count = imageCountParentObj[1]) %>% ungroup()
# 
# # Save the data files
# write_csv(data_p53_ps, file = "/data/muriel/Projects/PHH/DataAnalysis/GFPdata/CisReexpData.csv")

# Load summary data
data_p53_ps <- read_csv(file = "/data/muriel/Projects/PHH/DataAnalysis/GFPdata/CisReexpData.csv",
                         col_types = "cffccfdcdd")
data_p53_ps <- data_p53_ps %>%
  mutate(dose_uM = factor(dose_uM, levels = c(0,1,2.5,5,10,15,20,25,50)))

# Make a figure
ggplot(data = data_p53_ps) + 
  geom_point(aes(x = timeAfterExposure, y = IntegratedIntensity_GM, color = locationID)) +
  geom_line(aes(x = timeAfterExposure, y = IntegratedIntensity_GM, group = locationID, color = locationID)) +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  scale_color_viridis(discrete = T) + theme_classic() +
  facet_wrap(~dose_uM)

ggplot(data = data_p53_ps) + 
  geom_point(aes(x = timeAfterExposure, y = Cell_count, color = locationID)) +
  geom_line(aes(x = timeAfterExposure, y = Cell_count, group = locationID, color = locationID)) +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  scale_color_viridis(discrete = T) + theme_classic() +
  facet_wrap(~dose_uM)

# Calculate summaries per biological replicate, i.e. calculate average of 2 technical replicates
data_p53_ts <- data_p53_ps %>% group_by(treatment, dose_uM, replID, plateID, cell_line, 
                                        timeID, timeAfterExposure) %>% 
  summarise(IntegratedIntensity = mean(IntegratedIntensity_GM),
            MeanIntensity = mean(MeanIntensity_GM),
            Cell_count = Cell_count[1]) %>% ungroup

# Do background subtraction
control_p53 <- data_p53_ts %>% filter(treatment == "DMEM") %>% 
  dplyr::rename(DMEM_IntInt = IntegratedIntensity,
                DMEM_MeanInt = MeanIntensity) %>% select(-c(Cell_count,dose_uM, treatment)) 
data_p53_ts <- left_join(data_p53_ts,control_p53,
                         by = c("replID", "plateID", "cell_line", 
                                "timeID", "timeAfterExposure"), suffix = c("",""))
data_p53_ts <- data_p53_ts %>% mutate(IntIntCorr = IntegratedIntensity - DMEM_IntInt,
                                      MeanIntCorr = MeanIntensity - DMEM_MeanInt)
data_p53_ts <- data_p53_ts %>% group_by(treatment, dose_uM, replID, plateID, cell_line) %>% 
  mutate(IntIntCorr_tp1 = IntIntCorr - IntIntCorr[1],
         Cell_count_tp1 = Cell_count / Cell_count[1])

data_p53_ts$doseLabels <- factor(data_p53_ts$dose_uM, labels = c(expression(paste("0 ", mu,"M")),
                                                                 expression(paste("1 ", mu,"M")),
                                                                 expression(paste("2.5 ", mu,"M")),
                                                                 expression(paste("5 ", mu,"M")),
                                                                 expression(paste("10 ", mu,"M")),
                                                                 expression(paste("15 ", mu,"M")),
                                                                 expression(paste("20 ", mu,"M")),
                                                                 expression(paste("25 ", mu,"M")),
                                                                 expression(paste("50 ", mu,"M"))))

# Make a figure
ggplot(data = data_p53_ts) + 
  geom_point(aes(x = timeAfterExposure, y = IntegratedIntensity, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = IntegratedIntensity, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_tech_sum_raw_p53_IntInt.pdf"), width = 10, height = 1.75)

# Make a figure
ggplot(data_p53_ts %>% filter(dose_uM %in% c(0,1,2.5,5))) + 
  geom_point(aes(x = timeAfterExposure, y = IntIntCorr, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = IntIntCorr, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_tech_sum_norm_p53_IntInt_selection.pdf"), width = 8, height = 1.75)

# Make a figure
ggplot(data_p53_ts%>% filter(!dose_uM == 0)) + 
  geom_point(aes(x = timeAfterExposure, y = IntIntCorr, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = IntIntCorr, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"Fig3G_tech_sum_norm_p53_IntInt.pdf"), width = 8, height = 1.75)

# Make a figure
ggplot(data_p53_ts %>% filter(dose_uM %in% c(0,1,2.5,5))) + 
  geom_point(aes(x = timeAfterExposure, y = IntIntCorr_tp1, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = IntIntCorr_tp1, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_tech_sum_norm_to_t1_p53_IntInt_selection.pdf"), width = 10, height = 1.75)

ggplot(data_p53_ts %>% filter(!dose_uM == 0)) + 
  geom_point(aes(x = timeAfterExposure, y = IntIntCorr_tp1, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = IntIntCorr_tp1, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_tech_sum_norm_to_t1_p53_IntInt.pdf"), width = 10, height = 1.75)

# Make a figure
ggplot(data_p53_ts %>% filter(!dose_uM == 0)) + 
  geom_point(aes(x = timeAfterExposure, y = Cell_count_tp1, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = Cell_count_tp1, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  ylab("Cell count") + xlab("Time (h)") +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_tech_sum_norm_to_t1_p53_CellCount.pdf"), width = 10, height = 1.75)

# Calculate summaries per experiment replicate, i.e. calculate average of 3 biological replicates
data_p53_bs <- data_p53_ts %>% group_by(treatment, dose_uM, plateID, cell_line, 
                                        timeID, timeAfterExposure) %>% 
  summarise(IntInt = mean(IntegratedIntensity),
            MeanInt = mean(MeanIntensity),
            IntIntCorr = mean(IntIntCorr),
            MeanIntCorr = mean(MeanIntCorr),
            Cell_count = mean(Cell_count))
data_p53_bs <- data_p53_bs %>% group_by(treatment, dose_uM, plateID, cell_line) %>% 
  mutate(IntInt_tp1 = IntIntCorr - IntIntCorr[1],
            Cell_count_tp1 = Cell_count / Cell_count[1])

# Make a figure
ggplot() + 
  geom_point(data = data_p53_bs, aes(x = timeAfterExposure, y = IntInt, color = dose_uM)) +
  geom_line(data = data_p53_bs, aes(x = timeAfterExposure, y = IntInt, group = dose_uM, color = dose_uM)) +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  scale_color_viridis(discrete = T, name = expression("Concentration ("*mu*"M)")) + theme_classic()
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_biol_sum_raw_p53_IntInt.pdf"), width = 4, height = 3)

ggplot() + 
  geom_point(data = data_p53_bs, aes(x = timeAfterExposure, y = IntIntCorr, color = dose_uM)) +
  geom_line(data = data_p53_bs, aes(x = timeAfterExposure, y = IntIntCorr, group = dose_uM, color = dose_uM)) +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  scale_color_viridis(discrete = T, name = expression("Concentration ("*mu*"M)")) + theme_classic()
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_biol_sum_norm_p53_IntInt.pdf"), width = 4, height = 3)

ggplot() + 
  geom_point(data = data_p53_bs, aes(x = timeAfterExposure, y = Cell_count, color = dose_uM)) +
  geom_line(data = data_p53_bs, aes(x = timeAfterExposure, y = Cell_count, group = dose_uM, color = dose_uM)) +
  ylab("Cell count") + xlab("Time (h)") +
  scale_color_viridis(discrete = T, name = expression("Concentration ("*mu*"M)")) + theme_classic()
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_biol_sum_raw_p53_CellCount.pdf"), width = 4, height = 3)

ggplot() + 
  geom_point(data = data_p53_bs, aes(x = timeAfterExposure, y = Cell_count_tp1, color = dose_uM)) +
  geom_line(data = data_p53_bs, aes(x = timeAfterExposure, y = Cell_count_tp1, group = dose_uM, color = dose_uM)) +
  ylab("Cell count") + xlab("Time (h)") +
  scale_color_viridis(discrete = T, name = expression("Concentration ("*mu*"M)")) + theme_classic()
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_biol_sum_norm_to_t1_p53_CellCount.pdf"), width = 4, height = 3)

ggplot() + 
  geom_point(data = data_p53_bs, aes(x = timeAfterExposure, y = IntInt_tp1, color = dose_uM)) +
  geom_line(data = data_p53_bs, aes(x = timeAfterExposure, y = IntInt_tp1, group = dose_uM, color = dose_uM)) +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  scale_color_viridis(discrete = T, name = expression("Concentration ("*mu*"M)")) + theme_classic()
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_biol_sum_norm_to_t1_p53_IntInt.pdf"), width = 4, height = 3)

ggplot(data_p53_bs %>% filter(dose_uM %in% c(0,1,2.5,5,10,15))) + 
  geom_point(aes(x = timeAfterExposure, y = IntInt_tp1, color = dose_uM)) +
  geom_line(aes(x = timeAfterExposure, y = IntInt_tp1, group = dose_uM, color = dose_uM)) +
  #geom_smooth(aes(x = timeAfterExposure, y = IntInt_tp1, group = dose_uM, color = dose_uM), se = F) +
  ylab("p53-GFP intensity (a.u.)") + xlab("Time (h)") +
  scale_color_viridis(discrete = T, name = expression("Concentration ("*mu*"M)")) + theme_classic()
ggsave(paste0(OUTPUT_FIGPATH,"FigSX_biol_sum_norm_to_t1_p53_IntInt_selection.pdf"), width = 4, height = 3)


#### Do for pi ####
data_pi <- data_pi %>% select(-c(imageNumber,groupNumber,groupInd,imageID,
                                 groupNumber,plateWellID,control))
colnames(data_pi)
data_pi <- pivot_wider(data_pi, id_cols = c("treatment","dose_uM","plateID","timeID","timeAfterExposure","cell_line","replID","locationID"), 
                       names_from = "variable", values_from = "value")
data_pi <- data_pi %>% mutate(well = getWell(locationID))
data_pi_p53 <- data_pi %>% filter(cell_line == "HepG2 P53-GFP") %>% mutate(dose_uM = as.factor(dose_uM),
                                                                           replID = as.factor(replID))

ggplot(data_pi_p53) + 
  geom_point(aes(x = timeAfterExposure, y = count_PI_object_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0_, color = locationID)) +
  geom_line(aes(x = timeAfterExposure, y = count_PI_object_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0_, color = locationID)) +
  scale_color_viridis(discrete = T) + theme_classic() +
  ylim(0,1) +
  facet_wrap(~dose_uM)

# Calculate tech rep mean
data_pi_p53_ts <- data_pi_p53 %>% group_by(treatment, dose_uM, replID, plateID, cell_line, 
                                        timeID, timeAfterExposure, well) %>% 
  summarise(count_PI_positive = mean(count_PI_object_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0_),
            Cell_count = mean(imageCountParentObj, na.rm = T)) %>% ungroup()

ggplot(data_pi_p53_ts) + 
  geom_point(aes(x = timeAfterExposure, y = count_PI_positive, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = count_PI_positive, color = replID)) +
  scale_color_viridis(discrete = T) + theme_classic() +
  ylim(0,1) +
  facet_wrap(~dose_uM)

ggplot(data_pi_p53_ts) + 
  geom_point(aes(x = timeAfterExposure, y = Cell_count, color = replID)) +
  geom_line(aes(x = timeAfterExposure, y = Cell_count, color = replID)) +
  scale_color_viridis(discrete = T) + theme_classic() +
  facet_wrap(~dose_uM)

# Calculate biol rep mean
data_pi_p53_bs <- data_pi_p53_ts %>% group_by(treatment, dose_uM, plateID, cell_line, 
                                              timeID, timeAfterExposure) %>% 
  summarise(count_PI_positive = mean(count_PI_positive),
            Cell_count = mean(Cell_count)) %>% ungroup()

ggplot(data_pi_p53_bs) + # %>% filter(dose_uM %in% c(0,1,2.5,5))) + 
  geom_point(aes(x = timeAfterExposure, y = count_PI_positive, color = dose_uM)) +
  geom_line(aes(x = timeAfterExposure, y = count_PI_positive, color = dose_uM)) +
  theme_classic() +
  scale_color_viridis(discrete = T)

ggplot(data_pi_p53_bs) + 
  geom_point(aes(x = timeAfterExposure, y = Cell_count, color = dose_uM)) +
  geom_line(aes(x = timeAfterExposure, y = Cell_count, color = dose_uM)) +
  theme_classic() +
  scale_color_viridis(discrete = T)

  