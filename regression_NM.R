# linear regression of GeoMx gene and cell-types x Neuromelanin quantification.

library(ggpp)
library(readxl)
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
library(tidyr)

theme_set(theme_minimal())

## inputs
analysis_dir <- "/Users/zacc/USyd/NM_analysis"
setwd(analysis_dir)
quant_data = "/Users/zacc/USyd/NM_analysis/df_iNMeNM_230124.Rdata"
meta_data = "/Users/zacc/USyd/NM_analysis/NM_data_170124/ASAP1-IHC cohort.xlsx"

##############################################################
### Part 1: Extract NM quantification's  ###
##############################################################
# load NM quant data
load(quant_data)

# load metadata
meta <- read_xlsx(meta_data,1)
meta$Diagnosis[meta$Diagnosis == "Ct"] <- "CTR"

# # merge and format
df_agg_merged <- merge(df_agg,meta,by="Brainbank_ID",all.x=TRUE)
df_agg_merged$Diagnosis <- factor(df_agg_merged$Diagnosis,levels=c("CTR","ILBD", "ePD","lPD"))
df_agg_merged <- df_agg_merged[df_agg_merged$ROI != "ROI 1",]
df_agg_merged$log10_Area <- log10(df_agg_merged$`Area [µm²]`)
df_agg_merged$Sex[df_agg_merged$Sex == "1"] <- "M"
df_agg_merged$Sex[df_agg_merged$Sex == "2"] <- "F"
df_agg_merged$Diagnosis_stage <- as.character(df_agg_merged$Diagnosis)
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "CTR"] <- 0
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "ILBD"] <- 1
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "ePD"] <- 2
df_agg_merged$Diagnosis_stage[df_agg_merged$Diagnosis_stage == "lPD"] <- 3
df_agg_merged$Brainregion <- df_agg_merged$ROI
df_agg_merged$Brainregion[df_agg_merged$Brainregion %in% c("LC")] <- "A6"
df_agg_merged$Brainregion[df_agg_merged$Brainregion %in% c("SND","SNM","SNL","SNV")] <- "A9"
df_agg_merged$Brainregion[df_agg_merged$Brainregion %in% c("VTA")] <- "A10"
df_agg_merged$Direction.Red <- df_agg_merged$`Mean (Red)`/ (df_agg_merged$`Mean (Red)` + df_agg_merged$`Mean (Green)` + df_agg_merged$`Mean (Blue)`)*100
df_agg_merged$Direction.Green <- df_agg_merged$`Mean (Green)`/ (df_agg_merged$`Mean (Red)` + df_agg_merged$`Mean (Green)` + df_agg_merged$`Mean (Blue)`)*100
df_agg_merged$Direction.Blue <- df_agg_merged$`Mean (Blue)`/ (df_agg_merged$`Mean (Red)` + df_agg_merged$`Mean (Green)` + df_agg_merged$`Mean (Blue)`)*100
df_agg_merged <- df_agg_merged[complete.cases(df_agg_merged$log10_Area),]
df_agg_merged$Estimated.gray.value <- (0.299*df_agg_merged$`Mean (Red)`) + (0.587*df_agg_merged$`Mean (Green)`) + (0.114 * df_agg_merged$`Mean (Blue)`)
colnames(df_agg_merged) <- make.names(colnames(df_agg_merged))


# # format age group
# df_agg_merged$Age <- as.numeric(df_agg_merged$Age)
# dplot <- df_agg_merged
# dplot <- dplot %>%
#   mutate(
#     # Create categories
#     age_group = dplyr::case_when(
#       Age <= 70            ~ "<70",
#       Age > 70 & Age <= 80 ~ "70-80",
#       Age > 80 & Age <= 90 ~ "80-90",
#       Age > 90             ~ "> 90"
#     ),
#     # Convert to factor
#     age_group = factor(
#       age_group,
#       level = c("<70", "70-80","80-90", "> 90")
#     )
#   )
# df_agg_merged <- dplot

# format age group
df_agg_merged$Age <- as.numeric(df_agg_merged$Age)
dplot <- df_agg_merged
dplot <- dplot %>%
  mutate(
    # Create categories
    age_group = dplyr::case_when(
      Age <= 80 ~ "<80",
      Age > 80 & Age <= 90 ~ "80-90",
      Age > 90             ~ "> 90"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("<80","80-90", "> 90")
    )
  )
df_agg_merged <- dplot


# colour palettes
Diagnosis_col = c("CTR"= "grey","ILBD" = "#00AFBB", "ePD" = "#E7B800","lPD" = "red")
age_group_col = magma(4)
ROI_col = c("SNV" = "purple",
            "SNM" = viridis(6)[2],
            "SND" = viridis(6)[3],
            "SNL" = viridis(6)[4],
            "VTA" = viridis(6)[5],
            "LC" = viridis(6)[6],
            "SNV_young" = "dodgerblue")

RGB_col = c("Red" = "red",
            "Blue" = "blue",
            "Green" = "green")

# histograms of cohort
#dplot <- df_agg_merged[df_agg_merged$Age > 57,]
dplot <- df_agg_merged

pa <- ggplot(dplot, aes(x = log10(`Area [µm²]`))) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
  geom_density(alpha=.2, fill="grey") + xlim(0,4)

#dplot <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57),]
dplot <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" ),]

p1 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=fct_reorder(age_group, log10(`Area [µm²]`), .fun = mean))) +
  geom_density() + theme_minimal() + theme(legend.position = "none") + ylab("Age Group") + 
  geom_density_ridges(aes(fill = age_group)) + scale_fill_manual( values = age_group_col) + xlab("Log10(Area [µm²])")

p2 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=fct_reorder(Diagnosis, log10(`Area [µm²]`), .fun = mean))) +
  geom_density() + theme_minimal() + theme(legend.position = "none") + ylab("Diagnosis") + 
  geom_density_ridges(aes(fill = Diagnosis)) + scale_fill_manual( values = Diagnosis_col) + xlab("Log10(Area [µm²])")

p3 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=fct_reorder(ROI, log10(`Area [µm²]`), .fun = mean))) +
  geom_density() + theme_minimal() + theme(legend.position = "none") + ylab("Brain Region") + 
  geom_density_ridges(aes(fill = ROI)) + scale_fill_manual( values = ROI_col) + xlab("Log10(Area [µm²])")


arrange <- ggarrange(plotlist=list(p1,p2,p3), nrow=2, ncol=3, widths = c(2,2))
ggsave("distributions_iNMarea.pdf", arrange, width = 8, height = 6)


# group-wise mean and sd
data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57,]
colnames(data_table) <- make.names(colnames(data_table))

sd1 <- data_table %>%
  group_by(Diagnosis,ROI) %>%
  summarise_at(vars(Mean..Red.,Mean..Green.,Mean..Blue.,Direction.Red,Direction.Green, Direction.Blue, Mean..Hue.,Mean..Saturation.), list(sd=sd))


m1 <- data_table %>%
  group_by(Diagnosis,ROI) %>%
  summarise_at(vars(Mean..Red.,Mean..Green.,Mean..Blue.,Direction.Red,Direction.Green, Direction.Blue, Mean..Hue.,Mean..Saturation.), list(mean=mean))

tmp <- merge(m1,sd1,by.x=c("Diagnosis","ROI"),by.y=c("Diagnosis","ROI"))
write.table(tmp,file="iNM_mean.sd.txt",sep='\t',quote=F,row.names = F)

#####################################################
### Part 2.1: Evaluate the iNM in Controls by ROI ###
#####################################################
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57 ,]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

# run shapiro - not normal
shapiro.test(data_table$log10_Area[1:5000])

# linear mixed-effects model of iNM size v region
y_list <- c("log10_Area")
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = "Log10(Area [µm²])"
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black") +
    ylim(1.2,4.2)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  
  # write summary stats to file
  summary_stats <- data_table %>%
    group_by(ROI) %>%
    summarise(
      mean = mean(y, na.rm = TRUE),
      sd = sd(y, na.rm = TRUE)
    )
  
  write.table(as.data.frame(summary_stats),
              file=paste0("summary_stats_",make.names(y_lab),".txt"),
              quote=FALSE, row.names =FALSE)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
}

g1 <- bxp

# evaluate counts/ um2
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57 ,]
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, ROI,Age,Sex,PMD,`ROI Area [µm²]`, sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.mm2 <- (dplot$n / as.numeric(dplot$`ROI Area [µm²]`)) * 1000000
data_table <- dplot[complete.cases(dplot$ROI),]
colnames(data_table) <- make.names(colnames(data_table))


y_list <- c("n_per.mm2")
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = "iNM (No/mm²)"
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                       geom="pointrange", color="black") +
    ylim(0,130)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # write summary stats to file
  summary_stats <- data_table %>%
    group_by(ROI) %>%
    summarise(
      mean = mean(y, na.rm = TRUE),
      sd = sd(y, na.rm = TRUE)
    )
  
  write.table(as.data.frame(summary_stats),
              file=paste0("summary_stats_",make.names(y_lab),".txt"),
              quote=FALSE, row.names =FALSE)
  
  # # get sig values
  # if(length(stat.test$p.adj < 0.05) > 1){
  #   # add p-values to plot and save the result
  #   stat.test <- stat.test[stat.test$p.adj < 0.05, ]
  #   bxp <- bxp + stat_pvalue_manual(stat.test,
  #                                   y.position = max(data_table$y) * 1.6, step.increase = 0.1,
  #                                   label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  # }
  
  # save plot
  #ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
}

arrange <- ggarrange(plotlist=list(g1,bxp), nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_ctr_iNM.ROI.pdf", arrange,width = 8, height = 6)



#############################################################
### Part 2.2: Evaluate the iNM by PD group in Brainregion ###
#############################################################
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57),]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

data_table$total <- 1

# linear mixed-effects model of iNM size between diagnosis in each Brain region; area
y_list <- c("log10_Area")
res <- list()
roi_uniq <- c("A10","A9","A6")
main_lab <- c("VTA","SNpc","LC")
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$Brainregion == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = "Log10(Area [µm²])"
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(main_lab[z]) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(1.5,4)
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      print("test")
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = 3.8, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM.Brainregion.pdf", arrange, width = 8, height = 6)
#ggsave("Vln_iNM.ROI.pdf", arrange)


# linear mixed-effects model of iNM size between diagnosis in each Brain region; n/um2
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57),]
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, ROI, Brainregion,Age,Sex,PMD,ROI.Area..µm.., sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.mm2 <- (dplot$n / as.numeric(dplot$ROI.Area..µm..)) * 1000000
data_table <- dplot[complete.cases(dplot$Brainregion),]
colnames(data_table) <- make.names(colnames(data_table))

y_list <- c("n_per.mm2")
res <- list()
roi_uniq <- c("A10","A9","A6")
main_lab <- c("VTA","SNpc","LC")

roi_uniq <- unique(data_table$ROI)
main_lab <- roi_uniq
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$Brainregion == roi_cont,]
 
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = "iNM (No./mm²)"
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(main_lab[z]) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(0,200)
  
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = 250, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM_n_per.mm2.Brainregion.pdf", arrange,width = 8, height = 6)
#ggsave("Vln_iNM_n_per.mm2.ROI.pdf", arrange)



### linear mixed-effects model of iNM size between diagnosis in each ROI for PD total; n/um2 as % of CTR
# select and format data
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57),]
data_table$Diagnosis2 <- as.character(data_table$Diagnosis)
data_table$Diagnosis2[data_table$Diagnosis2 %in% c("ePD","lPD")] <- "PD"
data_table <- data_table[data_table$Diagnosis2 != "ILBD",]

# format measurements
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis2,Brainbank_ID_recode, ROI, Brainregion,Age,Sex,PMD,ROI.Area..µm.., sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.mm2 <- (dplot$n / as.numeric(dplot$ROI.Area..µm..)) * 1000000
data_table <- dplot[complete.cases(dplot$Brainregion),]
colnames(data_table) <- make.names(colnames(data_table))

# select columns
columns_to_calculate <- c("n_per.mm2")

# group by ROI and Diagnosis, then calculate the mean for each group
df <- data_table
group_means <- df %>%
  group_by(ROI, Diagnosis2) %>%
  summarise(across(all_of(columns_to_calculate), mean, na.rm = TRUE), .groups = 'drop')

# filter for CTR means
ctr_means <- group_means %>%
  filter(Diagnosis2 == "CTR")

# join the CTR means back with the original dataset to have control means for comparison
df_with_ctr_means <- df %>%
  left_join(ctr_means, by = c("ROI"), suffix = c("", "_ctr_mean"))

# calculate the percentage of the control mean for each Diagnosis within each segment and ROI
df_with_pct_of_ctr <- df_with_ctr_means %>%
  mutate(across(all_of(columns_to_calculate),
                ~(.x / get(paste0(cur_column(), "_ctr_mean")) * 100),
                .names = "pct_of_ctr_{.col}"))

# test correct prc have been calculated 
mean(df_with_pct_of_ctr$pct_of_ctr_n_per.mm2[df_with_pct_of_ctr$ROI == "SNV" & 
                                               df_with_pct_of_ctr$Diagnosis2 == "CTR"])

# assign to data_table
data_table <- df_with_pct_of_ctr

# regression between CTR and PD for each ROI
y_list <- c("pct_of_ctr_n_per.mm2")
res_contrasts <- list()
roi_uniq <- unique(data_table$ROI)
main_lab <- roi_uniq
for(z in 1:length(roi_uniq)){
  print(roi_uniq[z])
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis2"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = "iNM (No./mm²)"
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis2 = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(main_lab[z]) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(0,200)
    
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create dummy stat.test df
    stat.test <- data_table2 %>%
      t_test(y ~ x) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
    stat.test <- stat.test %>%
      add_xy_position(fun = "mean_sd", x = "x", dodge = 0.8) 
    
    # replace p-values of stat.test with linear mixed-effects model
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    tmp <- as.data.frame(cbind(group1,group2,p.adj))
    
    # get sig values
    if(all(stat.test$group1 == tmp$group1 & stat.test$group2 == tmp$group2)){
      stat.test <- stat.test[,-which(colnames(stat.test) %in% c("p"))]
      stat.test$p.adj <- as.numeric(tmp$p.adj)
      stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                       cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                                       symbols = c("***", "**", "*", ".", " "))
      
      summary_df <- data_table2 %>%
        group_by(x) %>%
        summarise(
          mean_y = mean(y, na.rm = TRUE),  # Calculate mean, remove NA values if any
          sd_y = sd(y, na.rm = TRUE)       # Calculate standard deviation, remove NA values if any
        )
      stat.test$mean_group1 <- as.numeric(mapvalues(stat.test$group1,summary_df$x,summary_df$mean_y))
      stat.test$mean_group2 <- as.numeric(mapvalues(stat.test$group2,summary_df$x,summary_df$mean_y))
      stat.test$sd_group1 <- as.numeric(mapvalues(stat.test$group1,summary_df$x,summary_df$sd_y))
      stat.test$sd_group2 <- as.numeric(mapvalues(stat.test$group2,summary_df$x,summary_df$sd_y))
      stat.test$contrast <- roi_cont
    }
    
    # assign to list
    res_contrasts[[z]] <- stat.test
  }
}

stat.summary <- do.call(rbind,res_contrasts)

# write summary to file
df <- apply(stat.summary,2,as.character)
write.table(df, file="Density.iNM_ctr.pd_stat.summary.txt", sep="\t",row.names = F, quote = F)

## plot in single violin plot
tmp <- aggregate(pct_of_ctr_n_per.mm2~ ROI , data_table, mean)
fact_lvls <- tmp[order(-tmp[,"pct_of_ctr_n_per.mm2"]),][,"ROI"]
data_table[,"ROI"] <- factor(data_table[,"ROI"], levels = fact_lvls)

g1 <- ggviolin(
  data_table, x = "ROI", y = "pct_of_ctr_n_per.mm2",
  fill = "Diagnosis2",
  palette = c("grey","red"),
  add = "mean_sd", size = 0.2) +
  theme(legend.position = "none") +
  geom_hline(yintercept=100, linetype="dashed", color = "grey") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10),
        plot.title = element_text(size = 10)) +
  ylab("iNM (No./mm²) % Controls") + xlab("") +
  ggtitle("") + ylim(0,450) 

arrange <- ggarrange(plotlist=list(g1), nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_iNM_n_per.mm2.ROI_CTRPD.pdf", arrange,width = 8, height = 6)
#ggsave("Vln_iNM_n_per.mm2.ROI.pdf", arrange)

#############################################################
### Part 2.3: Evaluate the eNM by PD group in Brainregion ###
#############################################################
# linear mixed-effects model of iNM size between diagnosis in each Brain region; n/um2
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "eNM"),]
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, Brainregion,Age,Sex,PMD,`ROI Area [µm²]`, sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.µm2 <- dplot$n / as.numeric(dplot$`ROI Area [µm²]`) * 1000000
data_table <- dplot[complete.cases(dplot$Brainregion),]
colnames(data_table) <- make.names(colnames(data_table))

y_list <- c("n_per.µm2")
res <- list()
roi_uniq <- unique(data_table$Brainregion)
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$Brainregion == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = "eNM (No./mm²)"
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(roi_cont) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(0,170)
    
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = max(data_table2$y) * 1.6, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    ggsave(paste("violin_eNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=3, ncol=3, widths = c(2,2))
ggsave("Vln_eNM_n_per.mm2.Brainregion.pdf", arrange,width = 8, height = 6)


####################################
### Part 3.1: Evaluate grey OD  ###
###################################
# histograms of cohort
#dplot <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57),]
dplot <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" ),]


pa <- ggplot(dplot, aes(x = Estimated.gray.value)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
  geom_density(alpha=.2, fill="grey") 


p1 <- ggplot(dplot, aes(x = Estimated.gray.value, y=fct_reorder(age_group, Estimated.gray.value, .fun = mean))) +
  geom_density() + theme_minimal() + theme(legend.position = "none") + ylab("Age Group") + xlab("Grey OD") +
  geom_density_ridges(aes(fill = age_group)) + scale_fill_manual( values = age_group_col)
  

p2 <- ggplot(dplot, aes(x = Estimated.gray.value, y=fct_reorder(Diagnosis, Estimated.gray.value, .fun = mean))) +
  geom_density() + theme_minimal() + theme(legend.position = "none") + ylab("Diagnosis") + xlab("Grey OD") +
  geom_density_ridges(aes(fill = Diagnosis)) + scale_fill_manual( values = Diagnosis_col)

p3 <- ggplot(dplot, aes(x = Estimated.gray.value, y=fct_reorder(ROI, Estimated.gray.value, .fun = mean))) +
  geom_density() + theme_minimal() + theme(legend.position = "none") + ylab("Brain Region") + xlab("Grey OD") +
  geom_density_ridges(aes(fill = ROI)) + scale_fill_manual( values = ROI_col)
  

arrange <- ggarrange(plotlist=list(pa,p1,p2,p3), nrow=2, ncol=2, widths = c(2,2))
ggsave("distributions_greyOD.pdf", arrange)


# linear mixed-effects model of grey OD
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57),]

y_list <- c("Estimated.gray.value")
res <- list()
#roi_uniq <- unique(data_table$ROI)
roi_uniq <- c("A10","A9","A6")
main_lab <- c("VTA","SNpc","LC")
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$Brainregion== roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = "Grey OD"
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(main_lab[z]) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(30,280)
    
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = 250, step.increase = 0.18,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
#ggsave("Vln_iNM.greyOD.Brainregion.pdf", arrange,width = 8, height = 6)
#ggsave("Vln_iNM.greyOD.ROI.pdf", arrange)


# linear mixed-effects model of grey OD in controls
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$Diagnosis == "CTR" & df_agg_merged$Age > 57),]
y_list <- c("Estimated.gray.value")
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = "Grey OD"
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black") + ylim(20,250)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # # get sig values
  # if(length(stat.test$p.adj < 0.05) > 1){
  #   # add p-values to plot and save the result
  #   stat.test <- stat.test[stat.test$p.adj < 0.05, ]
  #   bxp <- bxp + stat_pvalue_manual(stat.test,
  #                                   y.position = max(data_table$y) * 1.6, step.increase = 0.1,
  #                                   label = "p.adj.signif") 
  # }
  # 
  # save plot
  #ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
}

g1 <- bxp
arrange <- ggarrange(plotlist=list(g1), nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM.greyOD.ROI.CTR.pdf", arrange)

arrange <- ggarrange(plotlist=list(g1,res[[1]],NULL,res[[2]],res[[3]]), nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM.greyOD.ROI.CTR_Dx.Brainregion.pdf", arrange,width = 8, height = 6)

#############################################################
### Part 4.1: Evaluate color variables in Controls by ROI ###
#############################################################
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" & df_agg_merged$Age > 57 ,]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

# run shapiro - not normal
shapiro.test(data_table$log10_Area[1:5000])

# linear mixed-effects model of RGB v region in controls
y_list <- c("Mean..Red.","Mean..Green.","Mean..Blue.")
y_labs <- c("Mean Red (a.u.)","Mean Green (a.u.)","Mean Blue (a.u.)")
res <- list()
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = y_labs[i]
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black") + ylim(0,700)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = 250, step.increase = 0.16,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
  # assign plot
  res[[i]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM_RGB.CTR.ROI.pdf", arrange,width = 8, height = 6)


var.test(data_table$Mean..Green.,data_table$Mean..Red.)
var.test(data_table$Mean..Blue.,data_table$Mean..Green.)
var.test(data_table$Mean..Red.,data_table$Mean..Blue.)

var(data_table$Mean..Red.)
var(data_table$Mean..Green.)
var(data_table$Mean..Blue.)

# linear mixed-effects model of Direction RGB v region
y_list <- c("Direction.Red","Direction.Green","Direction.Blue")
res <- list()
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = paste0(gsub("Direction.","",y_variable), " in iNM (%)")
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black") +
    ylim(min(data_table[, y_variable])* 0.8,max(data_table[, y_variable] * 1.2))
    #ylim(28.5,37)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  #ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
  # assign plot
  res[[i]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM_directionRGB.CTR.ROI_1.pdf", arrange,width = 8, height = 6)

var.test(data_table$Direction.Red,data_table$Direction.Green)
var.test(data_table$Direction.Red,data_table$Direction.Blue)
var.test(data_table$Direction.Blue,data_table$Direction.Green)

# linear mixed-effects model of Hue/ Saturation v region
y_list <- c("Mean..Hue.","Mean..Saturation.")
res <- list()
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = paste0(gsub("Direction.","",y_variable), " in iNM (%)")
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black") +
    ylim(0,)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
  # assign plot
  res[[i]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM_HueSat.CTR.ROI.pdf", arrange,width = 8, height = 6)


#############################################################
### Part 4.1b: Evaluate color variables in Controls and young CTRs in SNV and LC ###
#############################################################
# select cases
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% c("SNV","LC") ,]

# format data table
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]
data_table$ROI_dx <- data_table$ROI
data_table$ROI_dx[data_table$Age < 57] <- "SNV_young"
colnames(data_table) <- gsub("Mean..","",colnames(data_table))
colnames(data_table)  <- gsub("\\.","",colnames(data_table))

# Fit a linear mixed-effects model and collect stats for plotting and writing to file
y_list <- c("Red","Green","Blue")
res_contrasts <- list()
for(i in 1:length(y_list)){
  print(y_list[i])
  # define variables
  x_variable = "ROI_dx"
  y_variable = y_list[i]
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI_dx = "Tukey"))
  summary(posthoc)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  
  # create dummy stat.test df
  stat.test <- data_table %>%
    t_test(y ~ x) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = "x", dodge = 0.8) 
  
  
  # replace p-values of stat.test with linear mixed-effects model
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  tmp <- as.data.frame(cbind(group1,group2,p.adj))
  
  # get sig values
  if(all(stat.test$group1 == tmp$group1 & stat.test$group2 == tmp$group2)){
      stat.test <- stat.test[,-which(colnames(stat.test) %in% c("p"))]
      stat.test$p.adj <- as.numeric(tmp$p.adj)
      stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                       cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                                       symbols = c("***", "**", "*", ".", " "))
      
      summary_df <- data_table %>%
        group_by(x) %>%
        summarise(
          mean_y = mean(y, na.rm = TRUE),  # Calculate mean, remove NA values if any
          sd_y = sd(y, na.rm = TRUE)       # Calculate standard deviation, remove NA values if any
        )
      stat.test$mean_group1 <- as.numeric(mapvalues(stat.test$group1,summary_df$x,summary_df$mean_y))
      stat.test$mean_group2 <- as.numeric(mapvalues(stat.test$group2,summary_df$x,summary_df$mean_y))
      stat.test$sd_group1 <- as.numeric(mapvalues(stat.test$group1,summary_df$x,summary_df$sd_y))
      stat.test$sd_group2 <- as.numeric(mapvalues(stat.test$group2,summary_df$x,summary_df$sd_y))
      stat.test$contrast <- y_list[i]
  }
  
  # assign to list
  res_contrasts[[i]] <- stat.test
}

stat.summary <- do.call(rbind,res_contrasts)

# write summary to file
df <- apply(stat.summary,2,as.character)
write.table(df, file="Mean.color_ctr_snv.lc.stat.summary.txt", sep="\t",row.names = F, quote = F)

# prepare data for plotting
data_table_long <- data_table %>%
  pivot_longer(
    cols = c("Red","Green","Blue"),  names_to = "Condition", values_to = "Value")

data_table_long$Condition <- as.factor(data_table_long$Condition)

# Create the violin plot
g1 <- ggboxplot(
  data_table_long, x = "Condition", y = "Value",
  fill = "ROI_dx",
  palette = ROI_col,
  add = "mean_sd", size = 0.5) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10,colour = c("Blue","Green","Red")),
        plot.title = element_text(size = 10)) +
  ylab("Mean Color (a.u)") + xlab("") +
  ggtitle("") + ylim(0,250) 


## Pie chart 
res_pie <- list()
y_list <- unique(data_table_long$ROI_dx)
for(i in 1:length(y_list)){
  group1 <- y_list[i]
  df <- data_table_long[data_table_long$ROI_dx == group1, ]
  df <- aggregate(df,list(df$Condition),mean)
  df$Color <- df$Group.1
  df$Value <- round(df$Value,0)
  
  res_pie[[i]] <- ggplot(df, aes(x = "", y = Value, fill = Color)) +
    geom_col(color = "black") +
    geom_text(aes(label = Value),
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    scale_fill_brewer() + 
    theme_void() + scale_fill_manual(values=c("blue","green","red")) +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) + ggtitle(group1)
  }

arrange <- ggarrange(plotlist=res_pie, nrow=2, ncol=3, widths = c(2,2))
ggsave("Pie_iNM_mean.colors.pdf", arrange,width = 8, height = 6)



## Direction
# select cases
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% c("SNV","LC") ,]

# format data table
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]
data_table$ROI_dx <- data_table$ROI
data_table$ROI_dx[data_table$Age < 57] <- "SNV_young"
colnames(data_table) <- gsub("Direction.","",colnames(data_table))

# Fit a linear mixed-effects model and collect stats for plotting and writing to file
y_list <- c("Red","Green","Blue")
res_contrasts <- list()
for(i in 1:length(y_list)){
  SCNAprint(y_list[i])
  # define variables
  x_variable = "ROI_dx"
  y_variable = y_list[i]
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI_dx = "Tukey"))
  summary(posthoc)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  
  # create dummy stat.test df
  stat.test <- data_table %>%
    t_test(y ~ x) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = "x", dodge = 0.8) 
  
  
  # replace p-values of stat.test with linear mixed-effects model
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  tmp <- as.data.frame(cbind(group1,group2,p.adj))
  
  # get sig values
  if(all(stat.test$group1 == tmp$group1 & stat.test$group2 == tmp$group2)){
    stat.test <- stat.test[,-which(colnames(stat.test) %in% c("p"))]
    stat.test$p.adj <- as.numeric(tmp$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    
    summary_df <- data_table %>%
      group_by(x) %>%
      summarise(
        mean_y = mean(y, na.rm = TRUE),  # Calculate mean, remove NA values if any
        sd_y = sd(y, na.rm = TRUE)       # Calculate standard deviation, remove NA values if any
      )
    stat.test$mean_group1 <- as.numeric(mapvalues(stat.test$group1,summary_df$x,summary_df$mean_y))
    stat.test$mean_group2 <- as.numeric(mapvalues(stat.test$group2,summary_df$x,summary_df$mean_y))
    stat.test$sd_group1 <- as.numeric(mapvalues(stat.test$group1,summary_df$x,summary_df$sd_y))
    stat.test$sd_group2 <- as.numeric(mapvalues(stat.test$group2,summary_df$x,summary_df$sd_y))
    stat.test$contrast <- y_list[i]
  }
  
  # assign to list
  res_contrasts[[i]] <- stat.test
}

stat.summary <- do.call(rbind,res_contrasts)

# write summary to file
df <- apply(stat.summary,2,as.character)
write.table(df, file="Direction.color_ctr_snv.lc.stat.summary.txt", sep="\t",row.names = F, quote = F)

# prepare data for plotting
data_table_long <- data_table %>%
  pivot_longer(
    cols = c("Red","Green","Blue"),  names_to = "Condition", values_to = "Value")

data_table_long$Condition <- as.factor(data_table_long$Condition)


# Create the violin plot
g2 <- ggboxplot(
  data_table_long, x = "Condition", y = "Value",
  fill = "ROI_dx",
  palette = ROI_col,
  add = "mean_sd", size = 0.5) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10,colour = c("Blue","Green","Red")),
        plot.title = element_text(size = 10)) +
  ylab("Color Direction (%)") + xlab("") +
  ggtitle("") + ylim(15,60)  


arrange <- ggarrange(plotlist=list(g1,NULL,g2), nrow=2, ncol=2, widths = c(2,2))
ggsave("Boxplots_iNM_colors.pdf", arrange,width = 8, height = 6)


## Pie chart 
res_pie <- list()
y_list <- unique(data_table_long$ROI_dx)
for(i in 1:length(y_list)){
  group1 <- y_list[i]
  df <- data_table_long[data_table_long$ROI_dx == group1, ]
  df <- aggregate(df,list(df$Condition),mean)
  df$Color <- df$Group.1
  df$Value <- round(df$Value,0)
  
  res_pie[[i]] <- ggplot(df, aes(x = "", y = Value, fill = Color)) +
    geom_col(color = "black") +
    geom_text(aes(label = Value),
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    scale_fill_brewer() + 
    theme_void() + scale_fill_manual(values=c("blue","green","red")) +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) + ggtitle(group1)
}

arrange <- ggarrange(plotlist=res_pie, nrow=2, ncol=3, widths = c(2,2))
ggsave("Pie_iNM_direction.colors.pdf", arrange,width = 8, height = 6)


## Plot correlation between colors and Area
g3 <- ggplot(data_table_long, aes(x = log10_Area, y = Value, color = Condition)) +
  geom_point(alpha=.1,shape = 16) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +  # Add linear regression lines
  scale_color_manual(values=c("blue", "green", "red")) +
  xlab("Log10(Area [µm²])") + ylab("Color Direction (%)") +
  theme_minimal() +   theme(legend.position = "none")

arrange <- ggarrange(plotlist=list(g3), nrow=2, ncol=2, widths = c(2,2))
ggsave("Plot_direction_size_ctrs.snv.lc.pdf", arrange,width = 8, height = 6)


cor.test(data_table$Red,data_table$log10_Area)[[3]]
cor.test(data_table$Red,data_table$log10_Area)[[4]]^2

cor.test(data_table$Green,data_table$log10_Area)[[3]]
cor.test(data_table$Green,data_table$log10_Area)[[4]]^2

cor.test(data_table$Blue,data_table$log10_Area)[[3]]
cor.test(data_table$Blue,data_table$log10_Area)[[4]]^2



# linear mixed-effects model of RGB v region in controls
y_list <- c("Mean..Red.","Mean..Green.","Mean..Blue.")

y_labs <- c("Mean Red (a.u.)","Mean Green (a.u.)","Mean Blue (a.u.)")
res <- list()
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = y_labs[i]
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black") + ylim(0,700)
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = 250, step.increase = 0.16,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
  # assign plot
  res[[i]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM_RGB.CTR.ROI.pdf", arrange,width = 8, height = 6)



############################################################################
### Part 4.2: Evaluate color variables in Controls vs PD by Brain region ###
############################################################################
# linear mixed-effects model of color variables and Diagnosis within each ROI; RGB
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM"  & df_agg_merged$Age > 57),]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

# linear mixed-effects model of RGB between Ctr and PD groups within each ROI; RGB
#y_list <- c("Mean..Red.","Mean..Green.","Mean..Blue.","Mean..Hue.","Mean..Saturation.","Direction.Red","Direction.Green","Direction.Blue")
y_list <- c("Direction.Red")
#y_list <- c("Direction.Blue")
#y_list <- c("Direction.Green")
res <- list()
res_stat <- list()

roi_uniq <- c("A10","A9","A6")
main_lab <- c("VTA","SNpc","LC")
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$Brainregion == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = paste0(gsub("Direction.","",y_variable), " in iNM (%)")
    colour_palette = Diagnosis_col
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(main_lab[z]) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(30,56)
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      print("test")
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = 36, step.increase = 0.15,
                                      label = "p.adj.signif") 
    }
    
    # save stat.test as dataframe
    stat.test$ROI <- roi_cont
    stat.test$y_variable <- y_variable
    res_stat[[z]] <- stat.test
    
    # save plot
    #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

# tmp <- do.call("rbind", res_stat)
# write.table(tmp,file = paste0("stat.test_",y_variable,".txt"), quote = F, sep="\t", row.names = F)

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_iNM_Direction.Red.Diagnosis_Brainregion.pdf", arrange, width = 8, height = 6)


# linear mixed-effects model of color variables and Diagnosis within each ROI; area
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM"  & df_agg_merged$Age > 57),]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]
data_table$CTR <- as.factor(data_table$Diagnosis == "CTR")

# linear mixed-effects model of RGB between Ctr and PD groups within each ROI; RGB
y_list <- c("Direction.Red")
#y_list <- c("Direction.Blue")
res <- list()
roi_uniq <- unique(data_table$ROI)
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "CTR"
    y_variable = y_list[i]
    x_lab = "CTR"
    y_lab = paste0("Weight of ",gsub("Direction.","",y_variable), " in ",roi_cont, " iNM (%)")
    colour_palette = c("grey","red")
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(CTR = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(roi_cont) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(25,60)
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      print("test")
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = 55, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
   #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
ggsave("Vln_Direction.Blue.CTRvNonCTR_ROI.pdf", arrange)


# # linear mixed-effects model of RGB by disease stage within each ROI; area
# data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM"),]
# colnames(data_table) <- make.names(colnames(data_table))
# data_table <- data_table[complete.cases(data_table$ROI),]
# data_table$Diagnosis_stage <- as.numeric(data_table$Diagnosis_stage)
# 
# # linear mixed-effects model of RGB between Ctr and ILBD within each ROI; area
# y_list <- c("Mean..Red.","Mean..Green.","Mean..Blue.")
# res <- list()
# roi_uniq <- unique(data_table$ROI)
# for(z in 1:length(roi_uniq)){
#   print(z)
#   roi_cont <- roi_uniq[z]
#   data_table2 <- data_table[data_table$ROI == roi_cont,]
#   for(i in 1:length(y_list)){
#     print(i)
#     # define variables
#     x_variable = "Diagnosis_stage"
#     y_variable = y_list[i]
#     x_lab = "Diagnosis stage"
#     y_lab = y_variable
#     colour_palette = Diagnosis_col
#     
#     # Fit a linear mixed-effects model
#     model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | PMD) + (1 | log10_Area)")), data = data_table2)
#     
#     # perform post-hoc tests to compare different regions using Tukey method
#     posthoc <- glht(model, linfct = mcp(Diagnosis_stage = "Tukey"))
#     summary(posthoc)
#     
#     # format for plotting
#     tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
#     fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#     data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
#     
#     # make violin plot 
#     bxp <- ggviolin(
#       data_table2, x = x_variable, y = y_variable, 
#       fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(roi_cont) + 
#       ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black")
#     
#     # set x and y in data_table2 for plotting
#     data_table2$x <- data_table2[, x_variable]
#     data_table2$y <- data_table2[, y_variable]
#     
#     # create stat.test df
#     summary_posthoc <- summary(posthoc)
#     group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#     group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#     p.adj = as.numeric(summary_posthoc$test$pvalues)
#     stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#     stat.test$p.adj <- as.numeric(stat.test$p.adj)
#     stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                      symbols = c("***", "**", "*", ".", " "))
#     print(stat.test)
#     # get sig values
#     if(any(stat.test$p.adj < 0.05)){
#       print("test")
#       # add p-values to plot and save the result
#       stat.test <- stat.test[stat.test$p.adj < 0.05, ]
#       bxp <- bxp + stat_pvalue_manual(stat.test,
#                                       y.position = max(data_table2$y) * 1.6, step.increase = 0.1,
#                                       label = "p.adj.signif") 
#     }
#     
#     
#     # save plot
#     ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
#   }
#   res[[z]] <- bxp
# }
# 
# arrange <- ggarrange(plotlist=res, nrow=2, ncol=2, widths = c(2,2))
# ggsave("Vln_RGB.Diagnosis_ROI.png", arrange)
# 


############################################################################
### Part 5: Principle component analysis of iNM size and color variables ###
############################################################################
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)


# select data and variables
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% c("SNV","LC") ) ,]
colnames(data_table) <- make.names(colnames(data_table))
#data_table$CTR <- as.factor(data_table$Diagnosis == "CTR")
# d2 <- data_table[,c("Mean (Green)","Mean (Blue)","Mean (Red)","Mean (Hue)",
#                     "Mean (Saturation)","log10_Area","Estimated.gray.value","Area [µm²]","Perimeter [µm]","Nearest Neighbor Distance [µm]")]
# d2 <- data_table[,c("Mean (Blue)","Mean (Hue)",
#                     "Mean (Saturation)","log10_Area","Nearest Neighbor Distance [µm]")]
d2 <- data_table[,c("log10_Area","Direction.Blue","Mean..Saturation.","Direction.Green","Direction.Red","Mean..Blue.","Mean..Red." )]


# # normalize data
# data_normalized <- scale(d2)
# 
# # compute correlation matrix
# corr_matrix <- cor(data_normalized)
# g1 <- ggcorrplot(corr_matrix)
# ggsave("Corr.matrix_PCA-min.pdf", g1)

# apply PCA
data.pca <- prcomp(d2, scale = TRUE)
summary(data.pca)

# plot explained variance
fviz_eig(data.pca, addlabels = TRUE)

groups <- data_table$Diagnosis
g2 <- fviz_pca_ind(data.pca,
             geom = c("point"),
             col.ind = groups, # color by groups
             palette = Diagnosis_col,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "CTR",
             repel = TRUE,
             alpha.ind = 0.5
)

ggsave("PCA_CTR_ROI_snvlc.pdf", g2)


#################################################
### Part 6: Correlations with SNCA expression ###
#################################################
# ST data
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata")
# extract expression values
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta <- as.data.frame(gxdat_s@meta.data)
meta$SNCA <- exp_dat["SNCA",]
meta <- meta[,c("Brainbank_ID","ROI","SNCA")]

# summarise SNCA per ROI per subject
df1 <- meta %>%
  group_by(Brainbank_ID, ROI) %>%
  summarize_at(vars(-group_cols()), mean, na.rm = TRUE)

df1$id <- paste(df1$Brainbank_ID,sep=".",df1$ROI)

# summarise iNM per ROI per subject
data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI == "SNV" & df_agg_merged$Age > 57 ,c("Brainbank_ID","ROI","Mean (Green)","Mean (Blue)","Mean (Red)","Mean (Hue)",
                                                                  "Mean (Saturation)","log10_Area","Estimated.gray.value","Area [µm²]","Perimeter [µm]",
                                                                  "Nearest Neighbor Distance [µm]")]
df <- data_table %>%
  group_by(Brainbank_ID, ROI) %>%
  summarize_at(vars(-group_cols()), mean, na.rm = TRUE)

df$id <- paste(df$Brainbank_ID,sep=".",df$ROI)

# merge
mdat <- merge(df,df1,by=c("id"))

# normalize data
data_normalized <- scale(mdat[,c(4:13,16)])

# compute correlation matrix
corr_matrix <- cor(data_normalized)
g1 <- ggcorrplot(corr_matrix)
ggsave("Corr.matrix_iNM.SNCA_snv.pdf", g1)


##########################################################
## Part 7: Evaluate color variables in Controls by ROI ###
##########################################################
data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" & df_agg_merged$Diagnosis == "CTR" & df_agg_merged$ROI %in% c("SNV","VTA"),]
colnames(data_table) <- make.names(colnames(data_table))
df <- data_table[complete.cases(data_table$ROI),]

# Convert categorical variables to factors
df$Sex <- as.factor(df$Sex)
df$ROI <- as.factor(df$ROI)
df$Brainbank_ID_recode <- as.factor(df$Brainbank_ID_recode)

df$Age <- as.numeric(df$Age)
df$PMD <- as.numeric(df$PMD)

# Split data into training and testing sets
set.seed(123)
training_indices <- sample(1:nrow(df), 0.7*nrow(df))
train_data <- df[training_indices, ]
test_data <- df[-training_indices, ]

# Multinomial logistic regression
multinom_model1 <- multinom(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + log10_Area + Mean..Hue. + Mean..Saturation. + 
                              Estimated.gray.value + Direction.Red + Direction.Green + Direction.Blue + Variance..Green. + 
                              Variance..Blue. + Variance..Red.  + Variance..Hue. + Variance..Saturation. , data = train_data)
multinom_model2 <- multinom(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + Mean..Hue. + Mean..Saturation. + log10_Area, data = train_data)
multinom_model3 <- multinom(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + Mean..Hue. + Mean..Saturation., data = train_data)
multinom_model4 <- multinom(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue., data = train_data)
multinom_model5 <- multinom(ROI ~ log10_Area, data = train_data)

model_summary <- summary(multinom_model1)
coefficients <- model_summary$coefficients / model_summary$standard.errors # Extract coefficients and their standard errors
p_values <- (1 - pnorm(abs(coefficients), 0, 1)) * 2 # Calculating p-values
as.data.frame(p_values)

# Random Forest
rf_model1 <- randomForest(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + log10_Area + Mean..Hue. + Mean..Saturation. + 
                            Estimated.gray.value + Direction.Red + Direction.Green + Direction.Blue + Variance..Green. + 
                            Variance..Blue. + Variance..Red.  + Variance..Hue. + Variance..Saturation. , data = train_data)
rf_model2 <- randomForest(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + Mean..Hue. + Mean..Saturation. + log10_Area, data = train_data)
rf_model3 <- randomForest(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + Mean..Hue. + Mean..Saturation., data = train_data)
rf_model4 <- randomForest(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue., data = train_data)
rf_model5 <- randomForest(ROI ~ log10_Area, data = train_data)

#importance(rf_model1)

# get df of accuracies
accuracy_data_list <- list()
multinom_mode_list <- list(multinom_model1,multinom_model2,multinom_model3,multinom_model4,multinom_model5)
rf_mode_list <- list(rf_model1 ,rf_model2,rf_model3 ,rf_model4 ,rf_model5)
for (i in 1:length(multinom_mode_list)){
  
  multinom_model <- multinom_mode_list[[i]]
  rf_model <- rf_mode_list[[i]]
  
  # Predictions
  multinom_predictions <- predict(multinom_model, test_data)
  rf_predictions <- predict(rf_model, test_data)
  
  # Evaluation Metrics
  log_conf <- confusionMatrix(multinom_predictions, test_data$ROI)
  rf_conf <- confusionMatrix(rf_predictions, test_data$ROI)
  
  log_accuracy <- sum(diag(log_conf$table)) / sum(log_conf$table)
  rf_accuracy <- sum(diag(rf_conf$table)) / sum(rf_conf$table)
  
  ## dataframe for plotting
  accuracy_data_list[[i]] <- data.frame(
    Model = c("Logistic Regression", "Random Forest"),
    Accuracy = c(log_accuracy, rf_accuracy))
  
}


accuracy_data <- do.call(rbind, args=accuracy_data_list)
accuracy_data$Covariates <- c(rep("All",2),rep("Area, RGB, HS",2),rep("RGB, HS",2),rep("RGB",2),rep("Area",2))
accuracy_data$Covariates <- factor(accuracy_data$Covariates,levels=c("Area","RGB","RGB, HS","Area, RGB, HS","All"))

# plot
g1 <- ggplot(accuracy_data, aes(x = Covariates, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "VTA v LC Model", x = "Covariates", y = "Accuracy") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  geom_text(aes(label = round(Accuracy, 2)), vjust = -0.5, size = 3.5,position = position_dodge(width=0.9))

arrange <- ggarrange(plotlist=list(g1), nrow=2, ncol=2, widths = c(2,2))
ggsave("bars_predmod_CTR.VTALC.pdf", arrange)


### apply best model to SNV
data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% c("SND","SNL","SNM","SNV") | df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% c("VTA","LC") & df_agg_merged$Diagnosis != "CTR" ,]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- rbind(data_table,test_data) # combine with test data not used for training
df <- data_table[complete.cases(data_table$ROI),]

# predict
rf_predictions <- predict(rf_model1, df)

# relabel predictions
df$predictions <- as.character(rf_predictions)
df$predictions[df$predictions == "LC"] <- "bad"
df$predictions[df$predictions == "VTA"] <- "good"


# linear mixed-effects model of iNM size between diagnosis in each Brain region; n/um2
data_table <- df
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, ROI, Brainregion,Age,Sex,PMD,predictions,ROI.Area..µm.., sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$predictions_per.mm2 <- (dplot$n / as.numeric(dplot$ROI.Area..µm..)) * 1000000
data_table <- dplot[complete.cases(dplot$ROI),]
colnames(data_table) <- make.names(colnames(data_table))

# select which species to plot
data_table <- data_table[data_table$predictions == "good",]

y_list <- c("predictions_per.mm2")
res <- list()
roi_uniq <- unique(data_table$ROI)
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab = "Predictions (No./mm²)"
    colour_palette = Diagnosis_col
    ymax <- max(data_table[,y_variable]) * 2
    yp <- max(data_table[,y_variable]) * 1.5
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(roi_cont) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(0,ymax)
    
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.05)){
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.05, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position =  yp, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))
#ggsave("Vln_good_predictions_per.mm2_ROI_vtalcmodel.pdf", arrange)



# ######################################################
# #### Part 8 - Evaluate the color channels vs sizes ###
# ######################################################
# density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# define region - LC
region1 <- c("LC")
# select data
dat <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% region1 & df_agg_merged$Age > 57),]
dat <- dat[complete.cases(dat$ROI),]
# format
colnames(dat) <- make.names(colnames(dat))

# Scale the data
data_scaled <- scale(dat[,c("log10_Area","Direction.Red","Direction.Blue","Direction.Green")])
# Run K-Means Clustering
set.seed(123) # For reproducibility of clustering
k <- 3 # Number of clusters
kmeans_result <- kmeans(data_scaled, centers = k)
dat$iNM_cluster <- kmeans_result$cluster
dat1 <- dat

# switch LC clusters to align with SNVd
dat1$iNM_cluster[dat1$iNM_cluster == "1"] <- "4" # flip clusters for LC
dat1$iNM_cluster[dat1$iNM_cluster == "3"] <- "1" # flip clusters for LC
dat1$iNM_cluster[dat1$iNM_cluster == "4"] <- "3" # flip clusters for LC

# define region - SNV
region1 <- c("SNV")
# select data
dat <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% region1 & df_agg_merged$Age > 57),]
dat <- dat[complete.cases(dat$ROI),]
# define young controls
dat_y <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM" & df_agg_merged$ROI %in% region1 & df_agg_merged$Age < 57),]
dat_y$iNM_cluster <- "young CTR"
# format
colnames(dat) <- make.names(colnames(dat))
colnames(dat_y) <- make.names(colnames(dat_y))
# Scale the data
data_scaled <- scale(dat[,c("log10_Area","Direction.Red","Direction.Blue","Direction.Green")])
# Run K-Means Clustering
set.seed(123) # For reproducibility of clustering
k <- 4 # Number of clusters
kmeans_result <- kmeans(data_scaled, centers = k)
dat$iNM_cluster<- kmeans_result$cluster
# add young ctr
dat2 <- rbind(dat,dat_y)

# switch SNV clusters to align with LC
dat2$iNM_cluster[dat2$iNM_cluster == "2"] <- "A" # flip clusters 
dat2$iNM_cluster[dat2$iNM_cluster == "4"] <- "2" # flip clusters 
dat2$iNM_cluster[dat2$iNM_cluster == "A"] <- "4" # flip clusters 

## combine and relabel clusters
dat3 <- rbind(dat1,dat2)
# dat3$iNM_cluster[dat3$iNM_cluster == "1"] <- "protective"
# dat3$iNM_cluster[dat3$iNM_cluster == "2"] <- "novel"
# dat3$iNM_cluster[dat3$iNM_cluster == "3"] <- "toxic"
# dat3$iNM_cluster[dat3$iNM_cluster == "4"] <- "neutral"
dat3$iNM_cluster[dat3$iNM_cluster == "1"] <- "non-toxic"
dat3$iNM_cluster[dat3$iNM_cluster == "2"] <- "transition"
dat3$iNM_cluster[dat3$iNM_cluster == "3"] <- "toxic"
dat3$iNM_cluster[dat3$iNM_cluster == "4"] <- "novel"

df_agg_iNM <- dat3

#save(df_agg_iNM,file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/nm_data_090224.Rdata")
save(df_agg_iNM,file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/nm_data_280224.Rdata")

## 3D plots
# create 3D html
dat <- dat1[order(dat1$iNM_cluster),] # assign relavent data frame
#cluster_colors <- c("purple","grey","hotpink","orange","cyan") # for SNV with 2 young controls
cluster_colors <- c("purple","grey", "hotpink") # for LC
dat$cluster <- as.character(dat$iNM_cluster)

# Make sure cluster is a factor
dat$cluster <- as.factor(dat$cluster)  # Make sure cluster is a factor

# create colors
names(cluster_colors) <- levels(dat$cluster)
RGB <- paste0("rgb(",round(dat$Mean..Red.,0),",",
              round(dat$Mean..Green.,0),",",
              round(dat$Mean..Blue.,0),")")

# Create your plotly plot
# fig <- plot_ly(data = dat, x = ~log10_Area, y = ~Direction.Red, z = ~Direction.Blue,
#                type = 'scatter3d', mode = 'markers',
#                marker = list(size = 4, color = ~RGB, line = list(width = 10)),
#                color = ~cluster, colors = cluster_colors)

fig <- plot_ly(data = dat, x = ~log10_Area, y = ~Direction.Red, z = ~Direction.Blue,
               type = 'scatter3d', mode = 'markers',
               marker = list(size = 4,line = list(width = 2, color = 'grey')),
               color = ~cluster, colors = cluster_colors)

# fig <- plot_ly(data = dat, x = ~Direction.Green, y = ~Direction.Red, z = ~Direction.Blue,
#                type = 'scatter3d', mode = 'markers',
#                marker = list(size = 4,line = list(width = 2, color = 'grey')),
#                color = ~cluster, colors = cluster_colors)

# fig <- plot_ly(data = dat, x = ~Direction.Green, y = ~Direction.Red, z = ~Direction.Blue,
#                type = 'scatter3d', mode = 'markers',
#                marker = list(size = ~Area..µm../70,  
#                              line = list(width = 2, color = 'grey')),
#                color = ~cluster, colors = cluster_colors)

# Customize the layout as needed
fig <- fig %>% layout(scene = list(xaxis = list(title = 'log10(Area)', titlefont = list(size = 24)),
                                   #xaxis = list(title = 'Direction Green'),
                                   yaxis = list(title = 'Direction Red', titlefont = list(size = 24)),
                                   zaxis = list(title = 'Direction Blue', titlefont = list(size = 24))))
# Show the plot
fig


## Make heatmaps for clustering
library(ComplexHeatmap)
library(circlize)

# convert labels to number
df_agg_iNM$iNM_cluster_number <- df_agg_iNM$iNM_cluster
df_agg_iNM$iNM_cluster_number[df_agg_iNM$iNM_cluster_number == "non-toxic"] <- 1
df_agg_iNM$iNM_cluster_number[df_agg_iNM$iNM_cluster_number == "transition"] <- 2
df_agg_iNM$iNM_cluster_number[df_agg_iNM$iNM_cluster_number == "toxic"] <- 3
df_agg_iNM$iNM_cluster_number[df_agg_iNM$iNM_cluster_number == "novel"] <- 4
df_agg_iNM$iNM_cluster_number[df_agg_iNM$iNM_cluster_number == "young CTR"] <- 5
df_agg_iNM$iNM_cluster_number <- as.numeric(df_agg_iNM$iNM_cluster_number)

# Define colors for the heatmap based on the kmeans cluster
cluster_colors <- get_cluster_colors(dat$iNM_cluster)

# scale data
mat <- scale(df_agg_iNM[,c("log10_Area", "Direction.Red", "Direction.Blue", "Direction.Green")])

# Change column names
new_colnames <- c("iNM Area", "Red", "Blue", "Green")

# Create right annotations
row_ha <- rowAnnotation(
  col = list(ROI = c("LC" = viridis(6)[6], "SNV" = "purple")),
  ROI = df_agg_iNM$ROI,
  show_legend = TRUE,  # Show legend for the annotation
  simple_anno_size = unit(1, "cm")
)

# Define a custom viridis color palette with modified alpha values
custom_viridis_palette <- viridis(n = 256, option = "A")  # Original viridis palette
custom_viridis_palette <- custom_viridis_palette[c(rep(1,5),seq(1,80,3),seq(81,121,1),seq(121,256,3))]

# Plot heatmap with custom column names
heatmap_object <- Heatmap(
  mat,
  name = "z-score",
  column_labels = new_colnames,
  row_split = df_agg_iNM$iNM_cluster,
  row_dend_reorder = TRUE,
  show_row_names = FALSE,
  row_dend_width = unit(5, "cm"),
  cluster_columns = FALSE,
  row_gap = unit(1, "mm"),
  #border = TRUE,
  right_annotation = row_ha,
  left_annotation = rowAnnotation(
    foo = anno_block(
      gp = gpar(fill = c("orange", "dodgerblue", "purple", "grey", "hotpink")),
      width = unit(3, "cm")
    )
  ),
  row_title = NULL,
  col = custom_viridis_palette
)

# print
pdf("heatmap_iNM.cluster.pdf",height = 3)
draw(heatmap_object)
dev.off()



### Logistic regression and Violin Plots
# calculate group-wise differences
data_table <- df_agg_iNM[!df_agg_iNM$iNM_cluster %in% "young CTR",]
colnames(data_table) <- make.names(colnames(data_table))
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, ROI, Brainregion,Age,Sex,PMD,iNM_cluster,ROI.Area..µm.., sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$case.region <- make.names(paste0(dplot$Brainbank_ID_recode,sep="_",dplot$ROI))

# select each clusetr
# d1 <- dplot[dplot$iNM_cluster == "protective",]
# d2 <- dplot[dplot$iNM_cluster == "novel",]
# d3 <- dplot[dplot$iNM_cluster == "toxic",]
# d4 <- dplot[dplot$iNM_cluster == "neutral",]
d1 <- dplot[dplot$iNM_cluster == "non-toxic",]
d2 <- dplot[dplot$iNM_cluster == "transition",]
d3 <- dplot[dplot$iNM_cluster == "toxic",]
d4 <- dplot[dplot$iNM_cluster == "novel",]

tmp <- merge(d1,d2[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n.x","n.y")] <- c("n.1","n.2")
tmp <- merge(tmp,d3[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n")] <- c("n.3")
tmp <- merge(tmp,d4[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n")] <- c("n.4")

tmp[is.na(tmp)] <- 0

# tmp$protective_prc.total <- tmp$n.1 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
# tmp$novel_prc.total <- tmp$n.2 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
# tmp$toxic_prc.total <- tmp$n.3 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
# tmp$neutral_prc.total <- tmp$n.4 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$nontoxic_prc.total <- tmp$n.1 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$transition_prc.total <- tmp$n.2 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$toxic_prc.total <- tmp$n.3 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$novel_prc.total <- tmp$n.4 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100

data_table <- tmp

#y_list <- c("protective_prc.total","novel_prc.total","toxic_prc.total","neutral_prc.total")
y_list <- c("nontoxic_prc.total","transition_prc.total","toxic_prc.total","novel_prc.total")
y_labs <- c("Non-toxic (%)","Transition (%)","Toxic (%)","Novel (%)")

# y_list <- c("nontoxic_prc.total","transition_prc.total","toxic_prc.total")
# y_labs <- c("Non-toxic (%)","Transition (%)","Toxic (%)")

res <- list()
#roi_uniq <- c("LC")
roi_uniq <- c("SNV")
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
  for(i in 1:length(y_list)){
    print(i)
    # define variables
    x_variable = "Diagnosis"
    y_variable = y_list[i]
    x_lab = "Diagnosis"
    y_lab =  y_variable
    colour_palette = Diagnosis_col
   # ymax <- max(data_table[,y_variable]) * 2
   # yp <- max(data_table[,y_variable]) * 1.5
    
    # Fit a linear mixed-effects model
    model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Age) + (1 | Sex) + (1 | PMD)")), data = data_table2)
    
    # perform post-hoc tests to compare different regions using Tukey method
    posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
    summary(posthoc)
    
    # format for plotting
    tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
    fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
    data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
    
    # make violin plot 
    bxp <- ggviolin(
      data_table2, x = x_variable, y = y_variable, 
      fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      ylab(y_labs[i]) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") + ylim(0,150)
    
    
    # set x and y in data_table2 for plotting
    data_table2$x <- data_table2[, x_variable]
    data_table2$y <- data_table2[, y_variable]
    
    # create stat.test df
    summary_posthoc <- summary(posthoc)
    group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
    group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
    p.adj = as.numeric(summary_posthoc$test$pvalues)
    stat.test <- as.data.frame(cbind(group1,group2,p.adj))
    stat.test$p.adj <- as.numeric(stat.test$p.adj)
    stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", " "))
    print(stat.test)
    # get sig values
    if(any(stat.test$p.adj < 0.1)){
      # add p-values to plot and save the result
      stat.test <- stat.test[stat.test$p.adj < 0.1, ]
      bxp <- bxp + stat_pvalue_manual(stat.test,
                                      y.position = 120, step.increase = 0.15,
                                      label = "p.adj.signif") 
    }
    
    res[[length(res) + 1]] <- bxp
    # save plot
    #ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  #res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=3, widths = c(2,2))

#ggsave("Vln_iNMcluster.prc.total.LC_v3.pdf", arrange,width = 8, height = 6)
#ggsave("Vln_iNMcluster.prc.total.SNV_v3.pdf", arrange,width = 8, height = 6)


## Stacked barcharts
# select clusters of interest
# df <- dat
# df <- df[!df$cluster %in% c("young CTR"),]

df <- df_agg_iNM
df <- df[!df$iNM_cluster %in% c("young CTR"),]
df <- df[df$ROI == "SNV",]

# Calculate mean RGB for each Diagnosis x CLuster
mean_rgb <- df %>%
  group_by(Diagnosis, iNM_cluster) %>%
  summarise(mean_R = mean(Mean..Red.), mean_G = mean(Mean..Green.), mean_B = mean(Mean..Blue.), .groups = 'drop')

# Calculate percentages
cluster_count <- df %>%
  group_by(Diagnosis, iNM_cluster) %>%
  summarise(count = n(), .groups = 'drop')

total_count <- df %>%
  group_by(Diagnosis) %>%
  summarise(total = n(), .groups = 'drop')

percentage_df <- merge(cluster_count, total_count, by = "Diagnosis")
percentage_df$percentage <- (percentage_df$count / percentage_df$total) * 100

# Merge percentage data with mean RGB values
plot_data <- merge(percentage_df, mean_rgb, by = c("Diagnosis", "iNM_cluster"))

# Define a color scheme for borders based on clusters
#border_colors <- c("1" = "purple", "2" = "grey", "3" = "hotpink", "4" = "orange") 
border_colors <- c("non-toxic" = "purple", "transition" = "grey", "toxic" = "hotpink", "novel" = "orange") 

# Convert mean RGB to hexadecimal colors for filling
plot_data$fill_color <- rgb(plot_data$mean_R/255, plot_data$mean_G/255, plot_data$mean_B/255)

# reorder
plot_data <- plot_data %>%
  arrange(Diagnosis, percentage) %>%
  mutate(iNM_cluster = factor(iNM_cluster, levels = c("novel","toxic","transition","non-toxic")))

# Now, update the plotting code with larger borders and ordered clusters
g1 <- ggplot(plot_data, aes(x = Diagnosis, y = percentage, fill = fill_color, group = iNM_cluster)) +
  geom_bar(stat = "identity", color = border_colors[as.character(plot_data$iNM_cluster)], 
           position = "stack", size = 0.5) + # Adjusted border size here
  scale_fill_identity() +
  labs(y = "Percentage", x = "Diagnosis", title = "iNM Clusters in SNV") +
  theme_minimal() +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), position = position_stack(vjust = 0.5), color = "black", size = 3.5)

g2 <- ggplot(plot_data, aes(x = Diagnosis, y = percentage, fill = iNM_cluster, group = iNM_cluster)) +
  geom_bar(stat = "identity", color = border_colors[as.character(plot_data$iNM_cluster)], 
           position = "stack", size = 1) + # Adjusted border size
  scale_fill_manual(values = border_colors) + # Use manual fill colors for clusters
  labs(y = "Percentage", x = "Diagnosis", title = "iNM Clusters in SNV") +
  theme_minimal() +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), position = position_stack(vjust = 0.5), color = "black", size = 3.5)

## save
res <- list(g1,g2)
arrange <- ggarrange(plotlist=res, nrow=2, ncol=2, widths = c(2,2))

ggsave("Bar_cluster.prc.class.SNV_k4_v3.pdf", arrange)

############################################################################################
#### Part X: tabulation of cohort statistics
############################################################################################
# Brightfield
dplot <- df_agg_merged[which(df_agg_merged$Age > 57),]

dplot$Sex[dplot$Sex %in% c("male","Male")] <- "M"
dplot$Sex[dplot$Sex %in% c("female","Female")] <- "F"
dplot$Diagnosis <- as.factor(dplot$Diagnosis)
tmp <- dplot

# Age
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_age = median(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_age = IQR(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(sd_age = sd(Age,na.rm=TRUE))

# Kruskal wallis test
tmp %>% kruskal_effsize(Age ~ Diagnosis)

# Dunns test
res <- tmp %>% 
  dunn_test(Age ~ Diagnosis, p.adjust.method = "BH") 
res

#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))

tmp2 <- tmp[!duplicated(tmp$Brainbank_ID),]

tmp3 <- tmp2[tmp2$Diagnosis %in% c("CTR","ILBD"),]
chisq.test(tmp3$Sex,tmp3$Diagnosis)

tmp3 <- tmp2[tmp2$Diagnosis %in% c("CTR","ePD"),]
chisq.test(tmp3$Sex,tmp3$Diagnosis)

tmp3 <- tmp2[tmp2$Diagnosis %in% c("CTR","lPD"),]
chisq.test(tmp3$Sex,tmp3$Diagnosis)


# postmortem interval
tmp$PMD <- as.numeric(tmp$PMD)
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = median(PMD,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = IQR(PMD,na.rm=TRUE))

# Kruskal wallis test
tmp %>% kruskal_effsize(PMD ~ Diagnosis)

# Dunns test
res <- tmp %>% 
  dunn_test(PMD ~ Diagnosis, p.adjust.method = "BH") 
res





# GeoMx
dplot <- df_agg_merged[which(df_agg_merged$Age > 57),]

dplot$Sex[dplot$Sex %in% c("male","Male")] <- "M"
dplot$Sex[dplot$Sex %in% c("female","Female")] <- "F"
dplot$Diagnosis <- as.factor(dplot$Diagnosis)
tmp <- dplot

# Age
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_age = median(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_age = IQR(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(sd_age = sd(Age,na.rm=TRUE))

# Kruskal wallis test
tmp %>% kruskal_effsize(Age ~ Diagnosis)

# Dunns test
res <- tmp %>% 
  dunn_test(Age ~ Diagnosis, p.adjust.method = "BH") 
res

#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))

tmp2 <- tmp[!duplicated(tmp$Brainbank_ID),]

tmp3 <- tmp2[tmp2$Diagnosis %in% c("CTR","ILBD"),]
chisq.test(tmp3$Sex,tmp3$Diagnosis)

tmp3 <- tmp2[tmp2$Diagnosis %in% c("CTR","ePD"),]
chisq.test(tmp3$Sex,tmp3$Diagnosis)

tmp3 <- tmp2[tmp2$Diagnosis %in% c("CTR","lPD"),]
chisq.test(tmp3$Sex,tmp3$Diagnosis)


# postmortem interval
tmp$PMD <- as.numeric(tmp$PMD)
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = median(PMD,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = IQR(PMD,na.rm=TRUE))

# Kruskal wallis test
tmp %>% kruskal_effsize(PMD ~ Diagnosis)

# Dunns test
res <- tmp %>% 
  dunn_test(PMD ~ Diagnosis, p.adjust.method = "BH") 
res











# ################################################################################################
# #### Part 5 - Evaluate the color channels relative to the expected by euNM and pheoNM colors ###
# ################################################################################################
# data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" & df_agg_merged$Diagnosis == "CTR" ,]
# colnames(data_table) <- make.names(colnames(data_table))
# data_table <- data_table[complete.cases(data_table$ROI),]
# #data_table_res <- data_table
# 
# # calculate euclidean distance to yellow and brown
# yellow = col2rgb('yellow')
# brown = col2rgb('brown')
# 
# # Function to calculate Euclidean distance between two colors
# color_distance <- function(color1, color2) {
#   sqrt(sum((color1 - color2) ^ 2))
# }
# 
# # Compare each sample to yellow and brown
# samples <- data_table[,c("Mean..Red.","Mean..Green.","Mean..Blue.")]
# res_y <- list()
# res_b <- list()
# res_class <- list()
# for (i in 1:nrow(samples)) {
#   sample <- samples[i, ]
#   distance_to_yellow <- color_distance(sample, yellow)
#   distance_to_brown <- color_distance(sample, brown)
#   
#   res_y[[i]] <- distance_to_yellow
#   res_b[[i]] <- distance_to_brown
#   
#   if (distance_to_yellow < distance_to_brown) {
#     cat("Sample", i, "is closer to Yellow\n")
#     res_class[[i]] <- "PH"
#   } else {
#     cat("Sample", i, "is closer to Brown\n")
#     res_class[[i]] <- "EU"
#   }
# }
# 
# # assign value
# data_table$distance_to_yellow <- unlist(res_y) 
# data_table$distance_to_brown <- unlist(res_b) 
# data_table$ratio_y.b <- (255 - unlist(res_y)) / (255 - unlist(res_b))
# hist(data_table$ratio_y.b)
# var_interest <- c("distance_to_yellow","distance_to_brown","ratio_y.b")
# y_list <- var_interest
# 
# # if eumelanin (EM) then will have high blue, with less green and less red
# for (i in 1:length(var_interest)) {
#   print(var_interest[i])
#   x <- data_table[,var_interest[i]]
#   
#   # Find the modes of a KDE
#   findmodes <- function(kde) {
#     kde$x[which(c(kde$y[-1],NA) < kde$y & kde$y > c(NA,kde$y[-length(kde$y)]))]
#   }
#   
#   # Compute mode trace, varying the bandwidth within a factor of 10
#   m <- mean(x)
#   id <- 1
#   bw <- density(x)$bw * 10^seq(1,-1, length.out=101) 
#   modes.lst <- lapply(bw, function(h) {
#     m.new <- sort(findmodes(density(x, bw=h)))
#     #  Associate each previous mode with a nearest new mode.
#     if (length(m.new)==1) delta <- Inf else delta <- min(diff(m.new))/2
#     d <- outer(m.new, m, function(x,y) abs(x-y))
#     i <- apply(d, 2, which.min)
#     g <- rep(NA_integer_, length(m.new))
#     g[i] <- id[1:ncol(d)]
#     # Create new ids for new modes that appear.
#     k <- is.na(g)
#     g[k] <- (sum(!k)+1):length(g)
#     id <<- g
#     m <<- m.new
#     data.frame(bw=h, Mode=m.new, id=g)
#   })
#   X <- do.call(rbind, args=modes.lst)
#   X$id <- factor(X$id)
#   
#   # Locate the modes at the most vertical portions of traces.
#   minslope <- function(x, y) {
#     f <- splinefun(x, y)
#     e <- diff(range(x)) * 1e-4
#     df2 <- function(x) ((f(x+e)-f(x-e)) / (2*e))^2 # Numerical derivative, squared
#     v <- optimize(df2, c(min(x),max(x)))
#     c(bw=v$minimum, slope=v$objective, Mode=f(v$minimum))
#   }
#   # Retain the desired modes.
#   n.modes <- 2 # USER SELECTED: Following visual assessment
#   bw.max <- max(subset(X, id==n.modes)$bw)
#   modes <- sapply(1:n.modes, function(i) {
#     Y <- subset(X, id==i & bw <= bw.max)
#     minslope(Y$bw, Y$Mode)
#   })
#   #
#   print(modes)
#   # Plot the results.
#   g1 <- ggplot(X, aes(bw, Mode)) +
#     geom_line(aes(col=id), size=1.2, show.legend=FALSE) +
#     geom_point(aes(bw, Mode), data=as.data.frame(t(modes)), size=3, col="Black", alpha=1/2) +
#     scale_x_log10() + ylab("Mode") + xlab("Bandwidth") +
#     coord_flip() + 
#     ggtitle(var_interest[i])
#   
#   
#   g2 <- ggplot(data.frame(x), aes(x, ..density..)) +
#     geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
#     geom_density(alpha=.2, fill="grey") + 
#     geom_vline(data=as.data.frame(t(modes)),
#                mapping=aes(xintercept=Mode), col="#D18A4e", size=1) +
#     xlab(paste0("(",var_interest[i],")")) + ylab("Density") 
#   
#   arrange <- ggarrange(plotlist=list(g1,g2), nrow=2, ncol=1, widths = c(2,2))
#   ggsave(paste0("Mode_Trace_test.",var_interest[i],".png"), arrange)
# }
# 
# # color v ROI
# for(i in 1:length(y_list)){
#   print(i)
#   # define variables
#   x_variable = "ROI"
#   y_variable = y_list[i]
#   x_lab = "ROI"
#   y_lab = y_variable
#   colour_palette = ROI_col 
#   
#   # Fit a linear mixed-effects model
#   model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | log10_Area)")), data = data_table)
#   
#   # perform post-hoc tests to compare different regions using Tukey method
#   posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
#   summary(posthoc)
#   
#   # format for plotting
#   tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
#   fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#   data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
#   
#   # make violin plot 
#   bxp <- ggviolin(
#     data_table, x = x_variable, y = y_variable, 
#     fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#   
#   # set x and y in data_table for plotting
#   data_table$x <- data_table[, x_variable]
#   data_table$y <- data_table[, y_variable]
#   
#   # create stat.test df
#   summary_posthoc <- summary(posthoc)
#   group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#   group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#   p.adj = as.numeric(summary_posthoc$test$pvalues)
#   stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#   stat.test$p.adj <- as.numeric(stat.test$p.adj)
#   stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                    symbols = c("***", "**", "*", ".", " "))
#   print(stat.test)
#   
#   # get sig values
#   if(length(stat.test$p.adj < 0.05) > 1){
#     # add p-values to plot and save the result
#     stat.test <- stat.test[stat.test$p.adj < 0.05, ]
#     bxp <- bxp + stat_pvalue_manual(stat.test,
#                                     y.position = max(data_table$y) * 1.6, step.increase = 0.1,
#                                     label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
#   }
#   
#   # save plot
#   ggsave(paste("violin_iNM.",y_variable,".png"), bxp)
#   
# }
# 
# 
# 
# ########################################################################
# #### Part 6.1 - Assessing distributions colour within iNM of Controls  ###
# ########################################################################
# # select iNM
# data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" & df_agg_merged$Diagnosis == "CTR",]
# colnames(data_table) <- make.names(colnames(data_table))
# data_table <- data_table[complete.cases(data_table$ROI),]
# 
# # color v region
# y_list <- c("Mean..Green.","Mean..Red.","Mean..Blue.",
#             "Direction.of.colour..G.","Direction.of.colour..R.","Direction.of.colour..B.")
# 
# for(i in 1:length(y_list)){
#   print(i)
#   # define variables
#   x_variable = "ROI"
#   y_variable = y_list[i]
#   x_lab = "ROI"
#   y_lab = y_variable
#   colour_palette = ROI_col 
#   
#   # Fit a linear mixed-effects model
#   model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Diagnosis) + (1 | Age) + (1 | Sex) + (1 | Sex) + (1 | log10_Area)")), data = data_table)
#   
#   # perform post-hoc tests to compare different regions using Tukey method
#   posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
#   summary(posthoc)
#   
#   # format for plotting
#   tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
#   fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#   data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
#   
#   # make violin plot 
#   bxp <- ggviolin(
#     data_table, x = x_variable, y = y_variable, 
#     fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#   
#   # set x and y in data_table for plotting
#   data_table$x <- data_table[, x_variable]
#   data_table$y <- data_table[, y_variable]
#   
#   # create stat.test df
#   summary_posthoc <- summary(posthoc)
#   group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#   group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#   p.adj = as.numeric(summary_posthoc$test$pvalues)
#   stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#   stat.test$p.adj <- as.numeric(stat.test$p.adj)
#   stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                    symbols = c("***", "**", "*", ".", " "))
#   print(stat.test)
#   
#   # get sig values
#   if(length(stat.test$p.adj < 0.05) > 1){
#     # add p-values to plot and save the result
#     stat.test <- stat.test[stat.test$p.adj < 0.05, ]
#     bxp <- bxp + stat_pvalue_manual(stat.test,
#                                     y.position = max(data_table$y) * 1.6, step.increase = 0.1,
#                                     label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
#   }
#   
#   # save plot
#   ggsave(paste("violin_iNM.",y_variable,".png"), bxp)
#   
# }
# 
# 
# # color v diagnosis
# y_list <- c("Mean..Green.","Mean..Red.","Mean..Blue.",
#             "Direction.of.colour..G.","Direction.of.colour..R.","Direction.of.colour..B.")
# 
# for(i in 1:length(y_list)){
#   print(i)
#   # define variables
#   x_variable = "Diagnosis"
#   y_variable = y_list[i]
#   x_lab = "Diagnosis"
#   y_lab = y_variable
#   colour_palette = Diagnosis_col
#   
#   # Fit a linear mixed-effects model
#   model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Diagnosis) + (1 | Age) + (1 | Sex) + (1 | log10_Area)")), data = data_table)
#   
#   # perform post-hoc tests to compare different regions using Tukey method
#   posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
#   summary(posthoc)
#   
#   # format for plotting
#   tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
#   fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#   data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
#   
#   # make violin plot 
#   bxp <- ggviolin(
#     data_table, x = x_variable, y = y_variable, 
#     fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#   
#   # set x and y in data_table for plotting
#   data_table$x <- data_table[, x_variable]
#   data_table$y <- data_table[, y_variable]
#   
#   # create stat.test df
#   summary_posthoc <- summary(posthoc)
#   group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#   group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#   p.adj = as.numeric(summary_posthoc$test$pvalues)
#   stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#   stat.test$p.adj <- as.numeric(stat.test$p.adj)
#   stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                    symbols = c("***", "**", "*", ".", " "))
#   print(stat.test)
#   # get sig values
#   if(any(stat.test$p.adj < 0.05)){
#     print("test")
#     # add p-values to plot and save the result
#     stat.test <- stat.test[stat.test$p.adj < 0.05, ]
#     bxp <- bxp + stat_pvalue_manual(stat.test,
#                                     y.position = max(data_table$y) * 1.6, step.increase = 0.1,
#                                     label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
#   }
#  
#   # save plot
#   ggsave(paste("violin_iNM.",y_variable,".",x_lab,".png"), bxp)
#   
# }
# 
# # color v diagnosis in each ROI
# y_list <- c("Mean..Green.","Mean..Red.","Mean..Blue.",
#             "Direction.of.colour..G.","Direction.of.colour..R.","Direction.of.colour..B.")
# 
# roi_uniq <- unique(data_table$ROI)
# for(z in 1:length(roi_uniq)){
#   print(z)
#   roi_cont <- roi_uniq[z]
#   data_table2 <- data_table[data_table$ROI == roi_cont,]
#   for(i in 1:length(y_list)){
#     print(i)
#     # define variables
#     x_variable = "Diagnosis"
#     y_variable = y_list[i]
#     x_lab = "Diagnosis"
#     y_lab = y_variable
#     colour_palette = Diagnosis_col
#     
#     # Fit a linear mixed-effects model
#     model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age) + (1 | Sex) + (1 | log10_Area)")), data = data_table2)
#     
#     # perform post-hoc tests to compare different regions using Tukey method
#     posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
#     summary(posthoc)
#     
#     # format for plotting
#     tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
#     fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#     data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
#     
#     # make violin plot 
#     bxp <- ggviolin(
#       data_table2, x = x_variable, y = y_variable, 
#       fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#     
#     # set x and y in data_table2 for plotting
#     data_table2$x <- data_table2[, x_variable]
#     data_table2$y <- data_table2[, y_variable]
#     
#     # create stat.test df
#     summary_posthoc <- summary(posthoc)
#     group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#     group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#     p.adj = as.numeric(summary_posthoc$test$pvalues)
#     stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#     stat.test$p.adj <- as.numeric(stat.test$p.adj)
#     stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                      symbols = c("***", "**", "*", ".", " "))
#     print(stat.test)
#     # get sig values
#     if(any(stat.test$p.adj < 0.05)){
#       print("test")
#       # add p-values to plot and save the result
#       stat.test <- stat.test[stat.test$p.adj < 0.05, ]
#       bxp <- bxp + stat_pvalue_manual(stat.test,
#                                       y.position = max(data_table2$y) * 1.6, step.increase = 0.1,
#                                       label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
#     }
#     
#     # save plot
#     ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
#     
#   }
# }
# 
# 
# ####################################################################
# #### Part 6.2 Assessing modalities of colour channels within iNM ### 
# ####################################################################
# #data_table_res <- data_table
# data_table <- data_table_res
# 
# # calculate euclidean distance to yellow and brown
# yellow = col2rgb('yellow')
# brown = col2rgb('brown')
# 
# 
# # Function to calculate Euclidean distance between two colors
# color_distance <- function(color1, color2) {
#   sqrt(sum((color1 - color2) ^ 2))
# }
# 
# # Compare each sample to yellow and brown
# samples <- data_table[,c("Mean..Red.","Mean..Green.","Mean..Blue.")]
# res_y <- list()
# res_b <- list()
# res_class <- list()
# for (i in 1:nrow(samples)) {
#   sample <- samples[i, ]
#   distance_to_yellow <- color_distance(sample, yellow)
#   distance_to_brown <- color_distance(sample, brown)
#   
#   res_y[[i]] <- distance_to_yellow
#   res_b[[i]] <- distance_to_brown
#   
#   if (distance_to_yellow < distance_to_brown) {
#     cat("Sample", i, "is closer to Yellow\n")
#     res_class[[i]] <- "PH"
#   } else {
#     cat("Sample", i, "is closer to Brown\n")
#     res_class[[i]] <- "EU"
#   }
# }
# 
# # assign value
# data_table$distance_to_yellow <- unlist(res_y) 
# data_table$distance_to_brown <- unlist(res_b) 
# data_table$ratio_y.b <- (255 - unlist(res_y)) / (255 - unlist(res_b))
# hist(data_table$ratio_y.b)
# var_interest <- c("distance_to_yellow","distance_to_brown","ratio_y.b")
# y_list <- var_interest
# 
# # if eumelanin (EM) then will have high blue, with less green and less red
# for (i in 1:length(var_interest)) {
#   print(var_interest[i])
#   x <- data_table[,var_interest[i]]
#   
#   # Find the modes of a KDE
#   findmodes <- function(kde) {
#     kde$x[which(c(kde$y[-1],NA) < kde$y & kde$y > c(NA,kde$y[-length(kde$y)]))]
#   }
#   
#   # Compute mode trace, varying the bandwidth within a factor of 10
#   m <- mean(x)
#   id <- 1
#   bw <- density(x)$bw * 10^seq(1,-1, length.out=101) 
#   modes.lst <- lapply(bw, function(h) {
#     m.new <- sort(findmodes(density(x, bw=h)))
#     #  Associate each previous mode with a nearest new mode.
#     if (length(m.new)==1) delta <- Inf else delta <- min(diff(m.new))/2
#     d <- outer(m.new, m, function(x,y) abs(x-y))
#     i <- apply(d, 2, which.min)
#     g <- rep(NA_integer_, length(m.new))
#     g[i] <- id[1:ncol(d)]
#     # Create new ids for new modes that appear.
#     k <- is.na(g)
#     g[k] <- (sum(!k)+1):length(g)
#     id <<- g
#     m <<- m.new
#     data.frame(bw=h, Mode=m.new, id=g)
#   })
#   X <- do.call(rbind, args=modes.lst)
#   X$id <- factor(X$id)
#   
#   # Locate the modes at the most vertical portions of traces.
#   minslope <- function(x, y) {
#     f <- splinefun(x, y)
#     e <- diff(range(x)) * 1e-4
#     df2 <- function(x) ((f(x+e)-f(x-e)) / (2*e))^2 # Numerical derivative, squared
#     v <- optimize(df2, c(min(x),max(x)))
#     c(bw=v$minimum, slope=v$objective, Mode=f(v$minimum))
#   }
#   # Retain the desired modes.
#   n.modes <- 2 # USER SELECTED: Following visual assessment
#   bw.max <- max(subset(X, id==n.modes)$bw)
#   modes <- sapply(1:n.modes, function(i) {
#     Y <- subset(X, id==i & bw <= bw.max)
#     minslope(Y$bw, Y$Mode)
#   })
#   #
#   print(modes)
#   # Plot the results.
#   g1 <- ggplot(X, aes(bw, Mode)) +
#     geom_line(aes(col=id), size=1.2, show.legend=FALSE) +
#     geom_point(aes(bw, Mode), data=as.data.frame(t(modes)), size=3, col="Black", alpha=1/2) +
#     scale_x_log10() + ylab("Mode") + xlab("Bandwidth") +
#     coord_flip() + 
#     ggtitle(var_interest[i])
#   
#   
#   g2 <- ggplot(data.frame(x), aes(x, ..density..)) +
#     geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
#     geom_density(alpha=.2, fill="grey") + 
#     geom_vline(data=as.data.frame(t(modes)),
#                mapping=aes(xintercept=Mode), col="#D18A4e", size=1) +
#     xlab(paste0("(",var_interest[i],")")) + ylab("Density") 
#   
#   arrange <- ggarrange(plotlist=list(g1,g2), nrow=2, ncol=1, widths = c(2,2))
#   ggsave(paste0("Mode_Trace_test.",var_interest[i],".png"), arrange)
# }
# 
# # color v diagnosis
# for(i in 1:length(y_list)){
#   print(i)
#   # define variables
#   x_variable = "Diagnosis"
#   y_variable = y_list[i]
#   x_lab = "Diagnosis"
#   y_lab = y_variable
#   colour_palette = Diagnosis_col
#   
#   # Fit a linear mixed-effects model
#   model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | ROI) + (1 | Age)")), data = data_table)
#   
#   # perform post-hoc tests to compare different regions using Tukey method
#   posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
#   summary(posthoc)
#   
#   # format for plotting
#   tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
#   fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#   data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
#   
#   # make violin plot 
#   bxp <- ggviolin(
#     data_table, x = x_variable, y = y_variable, 
#     fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#   
#   # set x and y in data_table for plotting
#   data_table$x <- data_table[, x_variable]
#   data_table$y <- data_table[, y_variable]
#   
#   # create stat.test df
#   summary_posthoc <- summary(posthoc)
#   group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#   group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#   p.adj = as.numeric(summary_posthoc$test$pvalues)
#   stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#   stat.test$p.adj <- as.numeric(stat.test$p.adj)
#   stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                    symbols = c("***", "**", "*", ".", " "))
#   print(stat.test)
#   # get sig values
#   if(any(stat.test$p.adj < 0.06)){
#     print("test")
#     # add p-values to plot and save the result
#     stat.test <- stat.test[stat.test$p.adj < 0.06, ]
#     bxp <- bxp + stat_pvalue_manual(stat.test,
#                                     y.position = max(data_table$y) * 1.6, step.increase = 0.1,
#                                     label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
#   }
#   
#   # save plot
#   ggsave(paste("violin_iNM.",y_variable,".",x_lab,".png"), bxp)
#   
# }
# 
# 
# # color v ROI
# for(i in 1:length(y_list)){
#   print(i)
#   # define variables
#   x_variable = "ROI"
#   y_variable = y_list[i]
#   x_lab = "ROI"
#   y_lab = y_variable
#   colour_palette = ROI_col 
#   
#   # Fit a linear mixed-effects model
#   model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Diagnosis) + (1 | Age)")), data = data_table)
#   
#   # perform post-hoc tests to compare different regions using Tukey method
#   posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
#   summary(posthoc)
#   
#   # format for plotting
#   tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
#   fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#   data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
#   
#   # make violin plot 
#   bxp <- ggviolin(
#     data_table, x = x_variable, y = y_variable, 
#     fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#   
#   # set x and y in data_table for plotting
#   data_table$x <- data_table[, x_variable]
#   data_table$y <- data_table[, y_variable]
#   
#   # create stat.test df
#   summary_posthoc <- summary(posthoc)
#   group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#   group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#   p.adj = as.numeric(summary_posthoc$test$pvalues)
#   stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#   stat.test$p.adj <- as.numeric(stat.test$p.adj)
#   stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                    symbols = c("***", "**", "*", ".", " "))
#   print(stat.test)
#   
#   # get sig values
#   if(length(stat.test$p.adj < 0.05) > 1){
#     # add p-values to plot and save the result
#     stat.test <- stat.test[stat.test$p.adj < 0.05, ]
#     bxp <- bxp + stat_pvalue_manual(stat.test,
#                                     y.position = max(data_table$y) * 1.6, step.increase = 0.1,
#                                     label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
#   }
#   
#   # save plot
#   ggsave(paste("violin_iNM.",y_variable,".png"), bxp)
#   
# }
# 
# 
# # color v Diagnosis per ROI
# roi_uniq <- unique(data_table$ROI)
# for(z in 1:length(roi_uniq)){
#   print(z)
#   roi_cont <- roi_uniq[z]
#   data_table2 <- data_table[data_table$ROI == roi_cont,]
#   for(i in 1:length(y_list)){
#     print(i)
#     # define variables
#     x_variable = "Diagnosis"
#     y_variable = y_list[i]
#     x_lab = "Diagnosis"
#     y_lab = y_variable
#     colour_palette = Diagnosis_col
#     
#     # Fit a linear mixed-effects model
#     model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID_recode) + (1 | Age)")), data = data_table2)
#     
#     # perform post-hoc tests to compare different regions using Tukey method
#     posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
#     summary(posthoc)
#     
#     # format for plotting
#     tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table2, mean)
#     fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
#     data_table2[,x_variable] <- factor(data_table2[,x_variable], levels = fact_lvls)
#     
#     # make violin plot 
#     bxp <- ggviolin(
#       data_table2, x = x_variable, y = y_variable, 
#       fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
#       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#     
#     # set x and y in data_table2 for plotting
#     data_table2$x <- data_table2[, x_variable]
#     data_table2$y <- data_table2[, y_variable]
#     
#     # create stat.test df
#     summary_posthoc <- summary(posthoc)
#     group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
#     group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
#     p.adj = as.numeric(summary_posthoc$test$pvalues)
#     stat.test <- as.data.frame(cbind(group1,group2,p.adj))
#     stat.test$p.adj <- as.numeric(stat.test$p.adj)
#     stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
#                                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                                      symbols = c("***", "**", "*", ".", " "))
#     print(stat.test)
#     # get sig values
#     if(any(stat.test$p.adj < 0.05)){
#       print("test")
#       # add p-values to plot and save the result
#       stat.test <- stat.test[stat.test$p.adj < 0.05, ]
#       bxp <- bxp + stat_pvalue_manual(stat.test,
#                                       y.position = max(data_table2$y) * 1.6, step.increase = 0.1,
#                                       label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
#     }
#     
#     # save plot
#     ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
#     
#   }
# }
