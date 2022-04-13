# Description:
# R script to analyse yH2AX foci data

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

INPUT_PATH <- "/data/muriel/Projects/PHH/DataAnalysis/yH2AX/Input/"
OUTPUT_FIGPATH <- "/data/muriel/Projects/PHH/DataAnalysis/yH2AX/Output/"
FILES <- dir(INPUT_PATH)

file1 <- read_tsv(file = paste0(INPUT_PATH,FILES[1]))
file2 <- read_tsv(file = paste0(INPUT_PATH,FILES[2]))

data <- bind_rows(file1,file2)
data <- pivot_wider(data = data,names_from = "variable", values_from = "value", 
                    id_cols = c("treatment","dose_uM","plateID","timeID","cell_line","replID","locationID"))
cis_data <- data %>% filter(treatment == "CIS") %>%
  mutate(dose_uM = factor(as.character(dose_uM), levels = c(0,2.5,5,10,25)))

# Normalise data
control_data <- cis_data %>% filter(dose_uM == 0) %>% select(-c("locationID", "dose_uM"))
cis_data <- left_join(cis_data, control_data, by=c("treatment","plateID","timeID","cell_line","replID"), suffix = c("","_control"))
cis_data <- cis_data %>% mutate(normFociCount = obj_nc_Children_obj_foci_Count-obj_nc_Children_obj_foci_Count_control)


# Get technical summary
cis_data_ts <- cis_data %>% group_by(treatment,dose_uM,plateID,timeID,cell_line,replID) %>%
  summarise(meanFociCountNorm = mean(normFociCount),
            sdFociCountNorm = sd(normFociCount),
            meanFociCount = mean(obj_nc_Children_obj_foci_Count),
            sdFociCount = sd(obj_nc_Children_obj_foci_Count))


# Make figure raw data
ggplot(data = cis_data_ts) + 
  geom_point(aes(x = timeID, y = meanFociCount, color = dose_uM), size = 1) +
  # geom_errorbar(aes(x = timeID, ymin = meanFociCount - sdFociCount, ymax = meanFociCount + sdFociCount, color = dose_uM),
  #               width = 2) +
  geom_line(aes(x = timeID, y = meanFociCount, group = dose_uM, color = dose_uM)) +
  xlab(" Time (h)") + ylab("Mean number of foci/cell") +
  scale_color_viridis_d(name = expression(paste("Cisplatin ("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(-1,10)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(OUTPUT_FIGPATH,"Foci_count_raw.pdf"), width = 3, height = 2)

# Make figure norm data
ggplot(data = cis_data_ts) + 
  geom_point(aes(x = timeID, y = meanFociCountNorm, color = dose_uM), size = 1) +
  # geom_errorbar(aes(x = timeID, ymin = meanFociCountNorm - sdFociCountNorm, ymax = meanFociCountNorm + sdFociCountNorm, color = dose_uM),
  #               width = 2) +
  geom_line(aes(x = timeID, y = meanFociCountNorm, group = dose_uM, color = dose_uM)) +
  xlab(" Time (h)") + ylab("Mean number of foci/cell") +
  scale_color_viridis_d(name = expression(paste("Cisplatin\n("*mu*"M)"))) + 
  theme_classic() +
  ggtitle("") +
  ylim(c(-1,10)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(OUTPUT_FIGPATH,"Foci_count_norm.pdf"), width = 3, height = 2)

colnames(data)
