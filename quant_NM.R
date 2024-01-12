# Scripts to quantify iNM and eNM from NM brightfield image analyses using 'TrueAI' as part of VS200 software

library(readxl)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(MASS)
theme_set(theme_minimal())

###############
### Inputs  ###
###############
analysis_dir <- "/Users/zacc/USyd/NM_analysis/"
NM_data_dir = "/Users/zacc/USyd/NM_analysis/NM_data_new"
setwd(analysis_dir)

############################################
### Part 1: Extract NM quantification's  ###
############################################
# list of NM quant files
file.names <- list.files(path=NM_data_dir,full.names=T)
ag_list <- list()
for (i in 1:length(file.names)){
  print(length(file.names)-i)
  print(file.names[i])
  df <- read_xlsx(file.names[i], 1)
  # remove rows without ROI eg. summary rows.
  df <- df[!is.na(df$ROI),]
  # add Brainbank_ID
  df$Brainbank_ID <- rep(gsub("_.*","",basename(file.names[i])),nrow(df))
  # select cols
  #df <- df[,c("Area [µm²]","ROI","Perimeter [µm]","Mean (Radius) [µm]","Mean (Gray Intensity Value)","Brainbank_ID")]
  df <- df[,c("Area [µm²]","ROI","Mean (Gray Intensity Value)","Brainbank_ID")]
  
  # process totals 
  df_roi <- read_xlsx(file.names[i], 3)
  df_roi <-  df_roi[!is.na(df_roi$ROI),]
  # select cols
  df_roi <- df_roi[,c("ROI","Sum (Area) [µm²]","Object Area Fraction ROI [%]","Mean (Mean (Gray Intensity Value))","Mean (Shape Factor)","ROI Area [µm²]")]
  # aggregate over ROI
  tmp <- aggregate(cbind(`Sum (Area) [µm²]`,`Mean (Mean (Gray Intensity Value))`,`ROI Area [µm²]`) ~ ROI, data = df_roi, FUN = sum, na.rm = TRUE)
  # mapvalues
  df$`Sum (Area) [µm²]` <- mapvalues(df$ROI,tmp$ROI,tmp$`Sum (Area) [µm²]`)
  df$`Mean (Mean (Gray Intensity Value))` <- mapvalues(df$ROI,tmp$ROI,tmp$`Mean (Mean (Gray Intensity Value))`)
  df$`ROI Area [µm²]` <- mapvalues(df$ROI,tmp$ROI,tmp$`ROI Area [µm²]`)
  
  # add to list
  ag_list[[i]] <- df
}

# collapse to dataframe
df_agg <- as.data.frame(do.call("rbind", ag_list))

# format
df_agg$log10_Area <- log10(df_agg$`Area [µm²]`)

######################################################
### Part 2: Assigning intra and extra-cellular NM ###
#####################################################
x <- log10(df_agg$`Area [µm²]`)

# k-means clustering on log10 reads
set.seed(20)
clusters <- kmeans(x, 2)
df_agg$cluster <- clusters$cluster

# calculate 95% confidence intervals for normal distribution of highest K-cluster by log10(unique reads)
# fit normal distribution
extra_NM_cluster = 2

fit <- fitdistr(x[df_agg$cluster == extra_NM_cluster], "normal")
para <- fit$estimate
a <- para[[1]]
s <- para[[2]]

# select 95% CI threshold
left_extra <- a-(1.96 * s)
right_extra <- a+(1.96 * s)

# fit normal distribution
intra_NM_cluster = 1

fit <- fitdistr(x[df_agg$cluster == intra_NM_cluster], "normal")
para <- fit$estimate
a <- para[[1]]
s <- para[[2]]

# select 95% CI threshold
left_intra <- a-(1.96 * s)
right_intra <- a+(1.96 * s)

# assign clusters
df_agg$intra.extra[df_agg$cluster == 2 & x < left_intra] <- "eNM"
df_agg$intra.extra[df_agg$cluster == 1 & x > right_extra] <- "iNM"

# plot data
p3 <- df_agg %>%
  ggplot(aes(log10(`Area [µm²]`), ..scaled.., fill=factor(intra.extra))) +
  geom_histogram(bins=20,aes(y=..count../sum(..count..))) +
  scale_fill_manual(values=c("grey","red","black")) +
  theme_minimal() + 
  xlim(0.5,3.5) + 
  labs(fill="") + geom_vline(xintercept = left_intra, color = "red", lty=2) +
  geom_vline(xintercept = right_extra, color = "grey", lty=2) +
  geom_density(aes(group=1)) 

ggsave("dist_iNMeNM.png", p3)

# save to file for downstream analysis
saveRDS(df_agg,file="df_iNMeNM.Rdata")
