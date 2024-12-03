# scripts to evaulate cell-type proportion differences in GeoMx following RCTD cell-type deconvolution for iNM manuscript

# libraries
library(edgeR)
library(dplyr)
library(plyr)
library(EnhancedVolcano)
library(ggplot2)
library(ggridges)
library(MASS)
library(ggplot2)
library(viridis)
library(ggpubr)
library(EnvStats)
library(rstatix)
library(ggridges)
library(lme4)
library(multcomp)
library(randomForest)
library(e1071)
library(caret)
library(nnet)
library(forcats)
library(ggcorrplot)
library(plotly)

# inputs
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/")
rctd_kamath = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/meta_rctd_Kamath_n63.Rdata"
rctd_webber = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124/meta_rctd_Webber_n10.Rdata"


ROI_col = c("SNV" = "purple",
            "SNM" = viridis(6)[2],
            "SND" = viridis(6)[3],
            "SNL" = viridis(6)[4],
            "VTA" = viridis(6)[5])

CA_col = c("CALB1_CALCR" = "purple",
                  "CALB1_CRYM_CCDC68" = "purple",
                  "CALB1_GEM" = "purple",   
                  "CALB1_PPP1R17" = "purple",
                  "CALB1_RBP4" = "purple",
                  "CALB1_TRHR" = "purple",      
                  "SOX6_AGTR1" = "dodgerblue",
                  "SOX6_DDT" = "dodgerblue",
                  "SOX6_GFRA2" = "dodgerblue",    
                  "SOX6_PART1" = "dodgerblue",
                  "NE" = "forestgreen")


# combine data from Kamath decon SN/VTA/RN and Webber decon LC
load(rctd_webber)
m_lc <- meta_rctd
load(rctd_kamath)
m_k <- meta_rctd

meta_dat <- merge(m_lc,m_k,by="row.names",keep = TRUE)
row.names(meta_dat) <- meta_dat$Row.names
meta_dat <- meta_dat[,-grep("*\\.y$", colnames(meta_dat))]
meta_dat <- gsub("*\\.x$","",colnames(meta_dat))

