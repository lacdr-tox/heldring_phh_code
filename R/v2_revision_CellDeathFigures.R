# Description:
# R script to make PI and AnV cell death plots
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

outputDir_GR <- "/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/Figures/"

DATE <- "20220411"

# Load the data
summaryDataNorm <- read_csv(file = "/data/muriel/Projects/PHH/DataAnalysis/GFPdata/SummaryDataNorm.csv",
                            col_types = paste0("ffccffddf",paste(rep("d",41-9),collapse = "")))

# Make factors
summaryDataNorm <- summaryDataNorm %>%
  mutate(dose_uMadj = factor(dose_uMadj, levels = c(0,1,2.5,5,10,15,20,25,50)),
         protein = factor(protein, levels = c("p53","MDM2","p21","BTG2")))

# Make exploring figure
ggplot(summaryDataNorm %>% filter(treatment %in% c("CDDP","DMEM"))) + 
  geom_point(aes(x = timeAfterExposure, y = propCount_PI, color = dose_uMadj)) +
  scale_color_viridis_d() +
  facet_grid(replID ~ protein)

ggplot(summaryDataNorm %>% filter(treatment %in% c("CDDP","DMEM"))) + 
  geom_point(aes(x = timeAfterExposure, y = propCount_AnV, color = dose_uMadj)) +
  scale_color_viridis_d() +
  facet_grid(replID ~ protein)

# Select first replicates
subset <- summaryDataNorm %>% filter(treatment %in% c("CDDP","DMEM")) %>% #, (replID == 1 | (replID == 2 & protein == "p53"))) %>% 
  group_by(dose_uMadj, protein, replID, treatment) %>% 
  summarise(max_prop_pi = max(propCount_PI),
            max_prop_anv = max(propCount_AnV)) %>%
  filter(max_prop_pi > 0 & !is.nan(max_prop_pi),
         max_prop_anv > 0 & !is.nan(max_prop_pi))

subset <- subset %>% group_by(dose_uMadj) %>%
  summarise(mean_prop_pi = mean(max_prop_pi, na.rm = T),
            sd_prop_pi = sd(max_prop_pi, na.rm = T),
            mean_prop_anv = mean(max_prop_anv, na.rm = T),
            sd_prop_anv = sd(max_prop_anv, na.rm = T))

ggplot(subset) +
  geom_point(aes(x = dose_uMadj, y = mean_prop_pi)) +
  geom_errorbar(aes(x = dose_uMadj, ymin = mean_prop_pi - sd_prop_pi, ymax = mean_prop_pi + sd_prop_pi), width = 0.2) +
  labs(x = expression("Cisplatin concentration ("*mu*"M)"), 
       y = "Proportion PI-positive") +
  theme_classic()
ggsave(paste0(outputDir_GR,DATE,"_GFPdata_Figs/PI_cell_death.pdf"),width = 3, height = 2.5)

ggplot(subset) +
  geom_point(aes(x = dose_uMadj, y = mean_prop_anv)) +
  geom_errorbar(aes(x = dose_uMadj, ymin = mean_prop_anv - sd_prop_anv, ymax = mean_prop_anv + sd_prop_anv), width = 0.2) +
  labs(x = expression("Cisplatin concentration ("*mu*"M)"), 
       y = "Proportion AnV-positive") +
  theme_classic()
ggsave(paste0(outputDir_GR,DATE,"_GFPdata_Figs/AnV_cell_death.pdf"),width = 3, height = 2.5)




