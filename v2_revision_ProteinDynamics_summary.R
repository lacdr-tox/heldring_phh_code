## Script for BioSpyder TempO-seq analysis and figure preparation for PHH manuscript
## Author: Muriel Heldring
## Date: 18 March, 2022

## Requirements:
## - Protein expression data (in 08a_CDDPdata_DDP.RData)
## - 


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

rm(list=ls())
gc()

# Working directories #
setwd("/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/")

DATE <- "20220411"

#Give path to output dir
outputDir_GR <- '/data/muriel/Projects/DDP/DataAnalysis/BioSpyder/Exp009_BioSpyder_Marije/DEGAnalysis/Output/Figures/'

# Input
PATH_TO_PROTDATA <- "/data/muriel/Projects/DDP/DataAnalysis/Dynamics/Output/RData/"

dir.create(paste0(outputDir_GR,DATE,"_GFPdata_Figs"))

###### INPUT USER ######


## Assembly of Data ##

# Load protein data
load(file = paste0(PATH_TO_PROTDATA,"08a_CDDPdata_DDP.RData"))

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

