# Description:
# R script to do the Analysis on the data with Nutlin exposure
# 1) p53, Mdm2 reporters
# 2) 3 replicates, and 
# 3) GFP and PI channels

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

# -----------------------------------------------------------------------------
# --------------------------LOAD GLOBAL VARIABLES------------------------------
# -----------------------------------------------------------------------------

# Load the functions
source(paste0("/data/muriel/Projects/","DDP","/DataAnalysis/Dynamics/","RScripts_ImagingDataAnalysis/01_Functions.R"))

INPUT_PATH <- "/data/muriel/Projects/PHH/DataAnalysis/NutlinExperiment/Input/"
OUTPUT_FIGPATH <- "/data/muriel/Projects/PHH/DataAnalysis/NutlinExperiment/Output/"
FILES <- dir(INPUT_PATH)
START_TIME <- 0 

# -----------------------------------------------------------------------------
# ----------------------------DO DATA ANALYSIS---------------------------------
# -----------------------------------------------------------------------------

# # Data loading
# all_data <- read_tsv(paste0(INPUT_PATH,"H5CP_Single_Cell.txt"))
# 
# # Add a low value (0.00001) to 0 measurements, which is necessary to calculate the geometric mean
# all_data <- all_data %>% 
#   mutate(Nuclei_Intensity_IntegratedIntensity_image_gfp = ifelse(Nuclei_Intensity_IntegratedIntensity_image_gfp == 0, 0.00001,Nuclei_Intensity_IntegratedIntensity_image_gfp),
#          Nuclei_Intensity_MeanIntensity_image_gfp = ifelse(Nuclei_Intensity_MeanIntensity_image_gfp == 0, 0.00001,Nuclei_Intensity_MeanIntensity_image_gfp))
# 
# # Summarise data into population summary
# all_data_ps <- all_data %>% select(-c(imageNumber,groupNumber,groupInd,imageID,Nuclei_AreaShape_Area,
#                                       groupNumber,plateWellID,Nuclei_Number_Object_Number)) %>%
#   mutate(time_h = timeAfterExposure + START_TIME) %>% 
#   group_by(treatment, dose_uM, replID, plateID, cell_line, 
#            timeID, time_h,timeAfterExposure, locationID, control) %>% 
#   summarise(IntegratedIntensity_GM = geomMean10(Nuclei_Intensity_IntegratedIntensity_image_gfp),
#             MeanIntensity_GM = geomMean10(Nuclei_Intensity_MeanIntensity_image_gfp),
#             Cell_count_mean = mean(imageCountParentObj),
#             Cell_count_sd = sd(imageCountParentObj),
#             PI_count = sum(count_Masked_PI_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0_)) %>% ungroup()
# # Check if imageCountParentObj sd = 0
# if (!all(all_data_ps %>% pull(Cell_count_sd) == 0)) {
#   warnings("Warning! Not all imageCountParentObj sd's are 0, which means that the data was not grouped correctly")
# } else {
#   all_data_ps <- all_data_ps %>% dplyr::rename("Cell_count" = "Cell_count_mean") %>% select(-c(Cell_count_sd))
# }
# 
# 
# all_data_ps <- all_data_ps %>% mutate(Nutlin_dose = sapply(control, function(x){str_split(x, "_")[[1]][2]}),
#                                       well = sapply(locationID, function(x){str_split(x, "_")[[1]][1]}),
#                                       dose_uM = ifelse(well %in% c("B08","C08","D08","E08","F08","G08",
#                                                                    "H08","I08","J08","K08","L08","M08",
#                                                                    "B17","C17","D17","E17","F17","G17"), 2.5, dose_uM),
#                                       plateRow = substr(well,1,1),
#                                       nutlin = ifelse(plateRow %in% c("H","I","J","K","L"), "Control", "Nutlin"))
# 
# # Save the summary file
# write_csv(all_data_ps, file = paste0(OUTPUT_FIGPATH,"NutlinSummary.csv"))
# 
# # Remove raw data
# rm(list = c("all_data"))

# Read in summary data
all_data_ps <- read_csv(file = paste0(OUTPUT_FIGPATH,"NutlinSummary.csv"))


# Make factors of time, dose_uM and replID
all_data_ps <- all_data_ps %>% mutate(dose_uM = as.factor(dose_uM),
                                      replID = as.factor(replID),
                                      timeID = as.factor(timeID),
                                      nutlin = as.factor(nutlin))

all_data_ps <- all_data_ps %>% filter(!treatment == "Cis1-test well-notuse")

# Calculate summaries per biological replicate, i.e. calculate average of 2 technical replicates
all_data_ts <- all_data_ps %>% group_by(treatment, dose_uM, replID, plateID, cell_line, 
                                        timeID, time_h, nutlin, well, Nutlin_dose) %>% 
  summarise(IntegratedIntensity = mean(IntegratedIntensity_GM),
            MeanIntensity = mean(MeanIntensity_GM),
            Cell_count = mean(Cell_count),
            PI_count = mean(PI_count)) %>% ungroup

# Do background subtraction
control <- all_data_ts %>% filter(nutlin == "Control") %>% 
  dplyr::rename(noNutlin_IntInt = IntegratedIntensity,
                noNutlin_MeanInt = MeanIntensity) %>% select(-c(Cell_count,nutlin,well, Nutlin_dose)) 
all_data_ts <- left_join(all_data_ts,control,
                         by = c("replID", "plateID", "cell_line",
                                "timeID", "time_h","dose_uM", "treatment"), suffix = c("",""))
all_data_ts <- all_data_ts %>% mutate(IntIntCorr = IntegratedIntensity - noNutlin_IntInt,
                                      MeanIntCorr = MeanIntensity - noNutlin_MeanInt)

all_data_ts <- all_data_ts %>% mutate(doseLabels = factor(dose_uM, labels = c(expression(paste("0 ", mu,"M")),
                                                                              expression(paste("1 ", mu,"M")),
                                                                              expression(paste("2.5 ", mu,"M")),
                                                                              expression(paste("5 ", mu,"M")),
                                                                              expression(paste("10 ", mu,"M")),
                                                                              expression(paste("15 ", mu,"M")),
                                                                              expression(paste("20 ", mu,"M")),
                                                                              expression(paste("25 ", mu,"M")),
                                                                              expression(paste("50 ", mu,"M")))),
                                      Nutlin_dose = factor(Nutlin_dose, levels = c(0,5,10)))

# Make a overview figure for p53
ggplot() + 
  geom_point(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_TP53"), aes(x = time_h, y = IntegratedIntensity, color = replID)) +
  geom_line(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_TP53"), aes(x = time_h, y = IntegratedIntensity, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  xlab("Time (h)") + 
  ylab("p53-GFP intensity (a.u.)") +
  facet_grid(Nutlin_dose~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"tech_sum_raw_p53_intint.pdf"), width = 10, height = 6)

# Make a overview figure for MDM2
ggplot() + 
  geom_point(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_MDM2"), aes(x = time_h, y = IntegratedIntensity, color = replID)) +
  geom_line(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_MDM2"), aes(x = time_h, y = IntegratedIntensity, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  xlab("Time (h)") + 
  ylab("p53-GFP intensity (a.u.)") +
  facet_grid(Nutlin_dose~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"tech_sum_raw_MDM2_intint.pdf"), width = 10, height = 6)

# Make a overview figure for p53
ggplot() + 
  geom_point(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_TP53"), aes(x = time_h, y = MeanIntensity, color = replID)) +
  geom_line(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_TP53"), aes(x = time_h, y = MeanIntensity, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  xlab("Time (h)") + 
  ylab("p53-GFP intensity (a.u.)") +
  facet_grid(Nutlin_dose~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"tech_sum_raw_p53_meanint.pdf"), width = 10, height = 6)

# Make a overview figure for MDM2
ggplot() + 
  geom_point(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_MDM2"), aes(x = time_h, y = MeanIntensity, color = replID)) +
  geom_line(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_MDM2"), aes(x = time_h, y = MeanIntensity, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  xlab("Time (h)") + 
  ylab("p53-GFP intensity (a.u.)") +
  facet_grid(Nutlin_dose~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"tech_sum_raw_MDM2_meanint.pdf"), width = 10, height = 6)

# Make a overview figure of the normalised data for MDM2
ggplot() + 
  geom_point(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_MDM2"), aes(x = time_h, y = MeanIntCorr, color = replID)) +
  geom_line(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_MDM2"), aes(x = time_h, y = MeanIntCorr, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  xlab("Time (h)") + 
  ylab("p53-GFP intensity (a.u.)") +
  facet_grid(Nutlin_dose~doseLabels, labeller = label_parsed)

# Make a overview figure of the normalised data for p53
ggplot() + 
  geom_point(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_TP53"), aes(x = time_h, y = MeanIntCorr, color = replID)) +
  geom_line(data = all_data_ts %>% filter(cell_line == "HepG2_GFP_TP53"), aes(x = time_h, y = MeanIntCorr, group = replID, color = replID)) +
  scale_color_viridis(discrete = T, name = "Exp.") + theme_classic() +
  xlab("Time (h)") + 
  ylab("p53-GFP intensity (a.u.)") +
  facet_grid(Nutlin_dose~doseLabels, labeller = label_parsed)

# Calculate summaries per experiment replicate, i.e. calculate average of 3 biological replicates
all_data_bs <- all_data_ts %>% group_by(treatment, dose_uM, doseLabels, plateID, cell_line, 
                                        timeID, time_h,nutlin, Nutlin_dose) %>% 
  summarise(IntInt = mean(IntegratedIntensity),
            MeanInt = mean(MeanIntensity),
            IntIntCorr = mean(IntIntCorr, na.rm = T),
            MeanIntCorr = mean(MeanIntCorr, na.rm = T),
            Cell_count = mean(Cell_count),
            sdIntInt = sd(IntegratedIntensity),
            sdMeanInt = sd(MeanIntensity),
            sdIntIntCorr = sd(IntIntCorr, na.rm = T),
            sdMeanIntCorr = sd(MeanIntCorr, na.rm = T),
            sdCell_count = sd(Cell_count)) %>% ungroup()

# Make a figure for P53
ggplot(data = all_data_bs %>% filter(cell_line == "HepG2_GFP_TP53")) + 
  geom_point(aes(x = time_h, y = IntInt, color = Nutlin_dose), size = 2) +
  geom_errorbar(aes(x = time_h, ymin = IntInt - sdIntInt, ymax = IntInt + sdIntInt, color = Nutlin_dose),
                width = 2) +
  geom_line(aes(x = time_h, y = IntInt, group = Nutlin_dose, color = Nutlin_dose)) +
  xlab(" Time (h)") + ylab("p53-GFP intensity (a.u.)") +
  scale_color_viridis(discrete = T, name = expression(paste("Nutlin ("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(0,3.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"biol_sum_raw_p53_intint.pdf"), width = 10, height = 2)

# Make a figure for P53
ggplot(data = all_data_bs %>% filter(cell_line == "HepG2_GFP_TP53", 
                                     dose_uM %in% c(0,1,2.5,5))) + 
  geom_point(aes(x = time_h, y = IntInt, color = Nutlin_dose), size = 2) +
  geom_errorbar(aes(x = time_h, ymin = IntInt - sdIntInt, ymax = IntInt + sdIntInt, color = Nutlin_dose),
                width = 2) +
  geom_line(aes(x = time_h, y = IntInt, group = Nutlin_dose, color = Nutlin_dose)) +
  xlab(" Time (h)") + ylab("p53-GFP intensity (a.u.)") +
  scale_color_viridis(discrete = T, name = expression(paste("Nutlin ("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(0,3.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"biol_sum_raw_p53_intint_selection.pdf"), width = 5, height = 2)

# Make a figure for MDM2
ggplot(data = all_data_bs %>% filter(cell_line == "HepG2_GFP_MDM2")) + 
  geom_point(aes(x = time_h, y = IntInt, color = Nutlin_dose), size = 2) +
  geom_errorbar(aes(x = time_h, ymin = IntInt - sdIntInt, ymax = IntInt + sdIntInt, color = Nutlin_dose),
                width = 2) +
  geom_line(aes(x = time_h, y = IntInt, group = Nutlin_dose, color = Nutlin_dose)) +
  xlab(" Time (h)") + ylab("MDM2-GFP intensity (a.u.)") +
  scale_color_viridis(discrete = T, name = expression(paste("Nutlin ("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(0,3.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"biol_sum_raw_MDM2_intint.pdf"), width = 10, height = 2)

# Make a figure for MDM2
ggplot(data = all_data_bs %>% filter(cell_line == "HepG2_GFP_MDM2", 
                                     dose_uM %in% c(0,1,2.5,5))) + 
  geom_point(aes(x = time_h, y = IntInt, color = Nutlin_dose), size = 2) +
  geom_errorbar(aes(x = time_h, ymin = IntInt - sdIntInt, ymax = IntInt + sdIntInt, color = Nutlin_dose),
                width = 2) +
  geom_line(aes(x = time_h, y = IntInt, group = Nutlin_dose, color = Nutlin_dose)) +
  xlab(" Time (h)") + ylab("MDM2-GFP intensity (a.u.)") +
  scale_color_viridis(discrete = T, name = expression(paste("Nutlin ("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(0,3.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"biol_sum_raw_MDM2_intint_selection.pdf"), width = 5, height = 2)


# Make a figure for P53
ggplot(data = all_data_bs %>% filter(cell_line == "HepG2_GFP_TP53", 
                                     dose_uM %in% c(0,1,2.5,5))) + 
  geom_point(aes(x = time_h, y = IntIntCorr, color = Nutlin_dose), size = 2) +
  geom_errorbar(aes(x = time_h, ymin = IntIntCorr - sdIntIntCorr, ymax = IntIntCorr + sdIntIntCorr, color = Nutlin_dose),
                width = 2) +
  geom_line(aes(x = time_h, y = IntIntCorr, group = Nutlin_dose, color = Nutlin_dose)) +
  xlab(" Time (h)") + ylab("p53-GFP intensity (a.u.)") +
  scale_color_viridis(discrete = T, name = expression(paste("Nutlin ("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(0,2.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"biol_sum_norm_p53_intint_selection.pdf"), width = 5, height = 2)


# Make a figure for MDM2
ggplot(data = all_data_bs %>% filter(cell_line == "HepG2_GFP_MDM2", 
                                     dose_uM %in% c(0,1,2.5,5))) + 
  geom_point(aes(x = time_h, y = IntIntCorr, color = Nutlin_dose), size = 2) +
  geom_errorbar(aes(x = time_h, ymin = IntIntCorr - sdIntIntCorr, ymax = IntIntCorr + sdIntIntCorr, color = Nutlin_dose),
                width = 2) +
  geom_line(aes(x = time_h, y = IntIntCorr, group = Nutlin_dose, color = Nutlin_dose)) +
  xlab(" Time (h)") + ylab("MDM2-GFP intensity (a.u.)") +
  scale_color_viridis(discrete = T, name = expression(paste("Nutlin ("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(0,2.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~doseLabels, labeller = label_parsed)
ggsave(paste0(OUTPUT_FIGPATH,"biol_sum_norm_MDM2_intint_selection.pdf"), width = 5, height = 2)

ggplot() + 
  geom_point(data = data_p53_bs, aes(x = time_h, y = IntInt, color = nutlin)) +
  geom_line(data = data_p53_bs, aes(x = time_h, y = IntInt, group = nutlin, color = nutlin)) +
  scale_color_viridis(discrete = T) + theme_classic() +
  facet_grid(~doseLabels, labeller = label_parsed)

ggplot() + 
  geom_point(data = data_p53_bs, aes(x = time_h, y = IntIntCorr, color = nutlin)) +
  geom_line(data = data_p53_bs, aes(x = time_h, y = IntIntCorr, group = nutlin, color = nutlin)) +
  scale_color_viridis(discrete = T) + theme_classic() +
  facet_grid(~doseLabels, labeller = label_parsed)

ggplot() + 
  geom_point(data = data_p53_bs, aes(x = time_h, y = MeanIntCorr, color = nutlin)) +
  geom_line(data = data_p53_bs, aes(x = time_h, y = MeanIntCorr, group = nutlin, color = nutlin)) +
  scale_color_viridis(discrete = T) + theme_classic() +
  facet_grid(~doseLabels, labeller = label_parsed)

