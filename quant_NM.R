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
NM_data_dir = "/Users/zacc/USyd/NM_analysis/NM_data_230124"
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

  # process totals 
  df_roi <- read_xlsx(file.names[i], 3)
  df_roi <-  df_roi[!is.na(df_roi$ROI),]

   # aggregate over ROI
  tmp <- aggregate(cbind(`Object Area Fraction ROI [%]`,	`Mean (Mean (Saturation))`,	`Mean (Variance (Saturation))`,	`Mean (Mean (Hue))`,	`Mean (Variance (Hue))`,	`Mean (Mean (Red))`,	`Mean (Variance (Red))`,	`Mean (Mean (Green))`,	`Mean (Variance (Green))`,	`Mean (Mean (Blue))`,	`Mean (Variance (Blue))`,	`ROI Area [µm²]`,	`Object Count [#]`) ~ ROI, data = df_roi, FUN = sum, na.rm = TRUE)
  
  # mapvalues
  df$`ROI Area [µm²]` <- mapvalues(df$ROI,tmp$ROI,tmp$`ROI Area [µm²]`)
  
  # add to list
  ag_list[[i]] <- df
}

# collapse to dataframe
df_agg <- as.data.frame(do.call("rbind", ag_list))
# format
df_agg$log10_Area <- log10(df_agg$`Area [µm²]`)

# ## OR ##
# # from one single file with each tab as data
# file.names <- "/Users/zacc/USyd/NM_analysis/NM_data_170124/IntracellularNMMidbrain 17:01:24.xlsx"
# sheets <- excel_sheets(path = file.names)
# ag_list <- list()
# for (i in 1:length(sheets)){
#   print(i)
#   df <- read_xlsx(file.names, sheets[i])
#   # remove rows without ROI eg. summary rows.
#   df <- df[!is.na(df$ROI),]
#   # add Brainbank_ID
#   df$Brainbank_ID <- sheets[i]
#   
#   # add to list
#   ag_list[[i]] <- df
# }
# ag_list1 <- ag_list
# 
# file.names <- "/Users/zacc/USyd/NM_analysis/NM_data_170124/IntracellularNMLC 17:01:24.xlsx"
# sheets <- excel_sheets(path = file.names)
# ag_list <- list()
# for (i in 1:length(sheets)){
#   print(i)
#   df <- read_xlsx(file.names, sheets[i])
#   # remove rows without ROI eg. summary rows.
#   df <- df[!is.na(df$ROI),]
#   # add Brainbank_ID
#   df$Brainbank_ID <- sheets[i]
#   
#   # add to list
#   ag_list[[i]] <- df
# }
# 
# # collapse to dataframe
# df_agg1 <- as.data.frame(do.call("rbind", ag_list1))
# df_agg <- as.data.frame(do.call("rbind", ag_list))
# df_agg <- rbind(df_agg,df_agg1)
# 
# # format
# df_agg$log10_Area <- log10(df_agg$`Area [µm²]`)

#########################################################################
### Part 2: Evaluating intra and extra-cellular NM thresholds per ROI ###
#########################################################################
df_agg_run <- df_agg

cutoff_eNM <- list()
cutoff_iNM <- list()
max_iNM <- list()
min_eNM <- list()

roi_uniq <- unique(df_agg_run$ROI)
for (i in 1:length(roi_uniq)){
  # k-means clustering on log10 reads
  set.seed(20)
  # get roi
  df_agg <- df_agg_run[df_agg_run$ROI == roi_uniq[i],]
  
  # get values
  x <- log10(df_agg$`Area [µm²]`)
  clusters <- kmeans(x, 2)
  df_agg$cluster <- clusters$cluster
  
  # get clusters mean
  ag_mean <- aggregate(x=x, by=list(df_agg$cluster ),
            FUN=mean)
  
  # calculate 95% confidence intervals for normal distribution of highest K-cluster by log10(unique reads)
  # fit normal distribution
  extra_NM_cluster = ag_mean$Group.1[ag_mean$x == min(ag_mean$x)]
  
  fit <- fitdistr(x[df_agg$cluster == extra_NM_cluster], "normal")
  para <- fit$estimate
  a <- para[[1]]
  s <- para[[2]]
  
  # select 95% CI threshold
  left_extra <- a-(1.96 * s)
  right_extra <- a+(1.96 * s)
  
  # fit normal distribution
  intra_NM_cluster = ag_mean$Group.1[ag_mean$x != min(ag_mean$x)]
  
  fit <- fitdistr(x[df_agg$cluster == intra_NM_cluster], "normal")
  para <- fit$estimate
  a <- para[[1]]
  s <- para[[2]]
  
  # select 95% CI threshold
  left_intra <- a-(1.96 * s)
  right_intra <- a+(1.96 * s)
 
  # assign clusters
  df_agg$intra.extra[df_agg$cluster == extra_NM_cluster & x < left_intra] <- "eNM"
  df_agg$intra.extra[df_agg$cluster == intra_NM_cluster & x > right_extra] <- "iNM"
# 
#   print(table(df_agg$cluster == intra_NM_cluster & x > log10(67)))
#   print(table(df_agg$cluster == intra_NM_cluster & x > log10(73)))
#   print(table(df_agg$cluster == intra_NM_cluster & x > log10(78)))
  #print(table(df_agg$cluster == extra_NM_cluster & x < log10(48)))
  #print(table(df_agg$cluster == extra_NM_cluster & x < log10(57)))
  print(table(df_agg$cluster == extra_NM_cluster & x < log10(75.4)))
  
  # assign list variables
  cutoff_eNM[[i]] <- left_intra
  cutoff_iNM[[i]] <- right_extra
  max_iNM[[i]] <- max(x[df_agg$intra.extra == "iNM"],na.rm = T)
  min_eNM[[i]] <- min(x[df_agg$intra.extra == "eNM"],na.rm = T)
  
  # plot data
  p3 <- df_agg %>%
    ggplot(aes(log10(`Area [µm²]`), ..scaled.., fill=factor(intra.extra))) +
    geom_histogram(bins=20,aes(y=..count../sum(..count..))) +
    scale_fill_manual(values=c("grey","red","black")) +
    theme_minimal() + 
    xlim(0.5,3.5) + 
    labs(fill="") + geom_vline(xintercept = left_intra, color = "red", lty=2) +
    geom_vline(xintercept = right_extra, color = "grey", lty=2) +
    geom_density(aes(group=1)) + ggtitle(roi_uniq[i])
  
  ggsave(paste0(roi_uniq[i],"dist_iNMeNM.png"), p3)

}

dat <- as.data.frame(cbind(ROI = roi_uniq, cutoff_eNM = 10^unlist(cutoff_eNM),
             cutoff_iNM = 10^unlist(cutoff_iNM), max_iNM = 10^unlist(max_iNM),
             min_eNM = 10^unlist(min_eNM)))

write.table(dat, file="dist_metrics.txt",sep="\t", quote = F, row.names = F)

######################################################
### Part 3: Assigning intra and extra-cellular NM ###
#####################################################
# We use the max threshold identified for iNM and min threshold to eNM
# k-means clustering on log10 reads
set.seed(20)
# get roi
df_agg <- df_agg_run

# get values
x <- log10(df_agg$`Area [µm²]`)
clusters <- kmeans(x, 2)
df_agg$cluster <- clusters$cluster

# get clusters mean
ag_mean <- aggregate(x=x, by=list(df_agg$cluster ),
                     FUN=mean)

# calculate 95% confidence intervals for normal distribution of highest K-cluster by log10(unique reads)
# fit normal distribution
extra_NM_cluster = ag_mean$Group.1[ag_mean$x == min(ag_mean$x)]

# fit normal distribution
intra_NM_cluster = ag_mean$Group.1[ag_mean$x != min(ag_mean$x)]

# assign clusters
df_agg$intra.extra <- NA
#df_agg$intra.extra[df_agg$cluster == extra_NM_cluster & x < log10(min(as.numeric(dat$cutoff_eNM)))] <- "eNM"
df_agg$intra.extra[df_agg$cluster == extra_NM_cluster & x < log10(sort(as.numeric(dat$cutoff_eNM))[2])] <- "eNM" # second highest eNM cutoff
df_agg$intra.extra[df_agg$cluster == intra_NM_cluster & x > log10(max(as.numeric(dat$cutoff_iNM)))] <- "iNM"

# plot data
p3 <- df_agg %>%
  ggplot(aes(log10(`Area [µm²]`), ..scaled.., fill=factor(intra.extra))) +
  geom_histogram(bins=200,aes(y=..count../sum(..count..) * 100)) +
  scale_fill_manual(values=c("grey","red","black")) +
  theme_minimal() +
  xlim(0.5,3.5) +
  labs(fill="") + geom_vline(xintercept = log10(max(as.numeric(dat$cutoff_iNM))), color = "red", lty=2) +
  geom_vline(xintercept = log10(sort(as.numeric(dat$cutoff_eNM))[2]), color = "grey", lty=2) +
  geom_density(aes(group=1))

ggsave("threshold_dist_iNMeNM.png", p3)

# save to file for downstream analysis
save(df_agg,file="df_iNMeNM_230124.Rdata")

