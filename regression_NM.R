# linear regression of GeoMx gene and cell-types x Neuromelanin quantification.

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

theme_set(theme_minimal())

## inputs
analysis_dir <- "/Users/zacc/USyd/NM_analysis"
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

# format age group
df_agg_merged$Age <- as.numeric(df_agg_merged$Age)
dplot <- df_agg_merged
dplot <- dplot %>% 
  mutate(
    # Create categories
    age_group = dplyr::case_when(
      Age <= 70            ~ "<70",
      Age > 70 & Age <= 80 ~ "70-80",
      Age > 80 & Age <= 90 ~ "80-90",
      Age > 90             ~ "> 90"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("<70", "70-80","80-90", "> 90")
    )
  )
df_agg_merged <- dplot

# colour palettes
Diagnosis_col = c("CTR"= "grey","ILBD" = "#00AFBB", "ePD" = "#E7B800","lPD" = "red")
age_group_col = magma(4)
ROI_col = c("SNV" = viridis(6)[1],
            "SNM" = viridis(6)[2],
            "SND" = viridis(6)[3],
            "SNL" = viridis(6)[4],
            "VTA" = viridis(6)[5],
            "LC" = viridis(6)[6])

# histograms of cohort
dplot <- df_agg_merged

pa <- ggplot(dplot, aes(x = log10(`Area [µm²]`))) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
  geom_density(alpha=.2, fill="grey") + xlim(0,4)

dplot <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM"),]

p1 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y = age_group)) +
  geom_density() + theme_minimal() +  theme(legend.position = "none") +
  geom_density_ridges(aes(fill = age_group))  + ylab("Age Group") + scale_fill_manual( values = age_group_col)

p2 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=Diagnosis)) +
  geom_density() + theme_minimal() + theme(legend.position = "none") +
  geom_density_ridges(aes(fill = Diagnosis)) + scale_fill_manual( values = Diagnosis_col)

p3 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=ROI)) +
  geom_density() + theme_minimal() + ylab("Brain Region") + theme(legend.position = "none") +
  geom_density_ridges(aes(fill = ROI)) + scale_fill_manual( values = ROI_col)


arrange <- ggarrange(plotlist=list(pa,p1,p2,p3), nrow=2, ncol=2, widths = c(2,2))
ggsave("distributions_iNMarea.png", arrange)


#####################################################
### Part 2.1: Evaluate the iNM in Controls by ROI ###
#####################################################
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" ,]
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
                                                                                     geom="pointrange", color="black")
  
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
  
}

g1 <- bxp

# evaluate counts/ um2
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" ,]
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, ROI,Age,Sex,PMD,`ROI Area [µm²]`, sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.µm2 <- dplot$n / as.numeric(dplot$`ROI Area [µm²]`)
data_table <- dplot[complete.cases(dplot$ROI),]
colnames(data_table) <- make.names(colnames(data_table))

y_list <- c("n_per.µm2")
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = "iNM (n/µm²)"
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
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                       geom="pointrange", color="black")
  
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
                                    label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  }
  
  # save plot
  ggsave(paste("violin_iNM.CTR",y_variable,".png"), bxp)
  
}

arrange <- ggarrange(plotlist=list(g1,bxp), nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_ctr_iNM.ROI.png", arrange)



#############################################################
### Part 2.2: Evaluate the iNM by PD group in Brainregion ###
#############################################################
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM"),]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

# linear mixed-effects model of iNM size between diagnosis in each Brain region; area
y_list <- c("log10_Area")
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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(roi_cont) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black")
    
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
                                      y.position = max(data_table2$y) * 1.6, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_iNM.Brainregion.png", arrange)


# linear mixed-effects model of iNM size between diagnosis in each Brain region; n/um2
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM"),]
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, Brainregion,Age,Sex,PMD,`ROI Area [µm²]`, sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.µm2 <- dplot$n / as.numeric(dplot$`ROI Area [µm²]`)
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
    y_lab = "iNM (n/µm²)"
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
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black")
  
    
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
    ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_iNM_n_per.µm2.Brainregion.png", arrange)


#############################################################
### Part 3.1: Evaluate color variables in Controls by ROI ###
#############################################################
data_table <- df_agg_merged[df_agg_merged$Diagnosis == "CTR" & df_agg_merged$intra.extra == "iNM" ,]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

# run shapiro - not normal
shapiro.test(data_table$log10_Area[1:5000])

# linear mixed-effects model of RGB v region
y_list <- c("Mean..Red.","Mean..Green.","Mean..Blue.")
res <- list()
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = y_variable
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
                                                                                     geom="pointrange", color="black")
  
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

arrange <- ggarrange(plotlist=res, nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_iNM_RGB.CTR.ROI.png", arrange)


# linear mixed-effects model of Hue/ Saturation v region
y_list <- c("Mean..Hue.","Mean..Saturation.")
res <- list()
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = y_variable
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
                                                                                     geom="pointrange", color="black")
  
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

arrange <- ggarrange(plotlist=res, nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_iNM_HueSat.CTR.ROI.png", arrange)


############################################################################
### Part 3.2: Evaluate color variables in Controls vs PD by Brain region ###
############################################################################
# linear mixed-effects model of color variables and Diagnosis within each ROI; area
data_table <- df_agg_merged[which(df_agg_merged$intra.extra == "iNM"),]
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[complete.cases(data_table$ROI),]

# linear mixed-effects model of RGB between Ctr and ILBD within each ROI; area
y_list <- c("Mean..Red.","Mean..Green.","Mean..Blue.","Mean..Hue.","Mean..Saturation.","Direction.Red","Direction.Green","Direction.Blue")
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
    y_lab = y_variable
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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(roi_cont) + 
      ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black")
    
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
                                      y.position = max(data_table2$y) * 1.6, step.increase = 0.1,
                                      label = "p.adj.signif") 
    }
    
    
    # save plot
    ggsave(paste("violin_iNM.",y_variable,".",x_lab,".",roi_cont,".png"), bxp)
  }
  res[[z]] <- bxp
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_RGB.Hue.Sat.Diagnosis_ROI.png", arrange)


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


###########################################################
### Part 4: Evaluate color variables in Controls by ROI ###
###########################################################
# data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" & df_agg_merged$Diagnosis == "CTR" ,]
# colnames(data_table) <- make.names(colnames(data_table))
# df <- data_table[complete.cases(data_table$ROI),]
# 
# # Convert categorical variables to factors
# df$Sex <- as.factor(df$Sex)
# df$Age <- as.factor(df$Age)
# df$ROI <- as.factor(df$ROI) 
# df$PMD <- as.factor(df$PMD) 
# df$Brainbank_ID_recode <- as.factor(df$Brainbank_ID_recode)
# 
# # Split data into training and testing sets
# set.seed(123) 
# training_indices <- sample(1:nrow(df), 0.7*nrow(df))
# train_data <- df[training_indices, ]
# test_data <- df[-training_indices, ]
# 
# # Multinomial logistic regression
# multinom_model <- multinom(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + Sex + Age + PMD + Brainbank_ID_recode, data = train_data)
# 
# # Random Forest
# rf_model <- randomForest(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + Sex + Age + Brainbank_ID_recode, data = train_data)
# 
# # Support Vector Machine
# svm_model <- svm(ROI ~ Mean..Red. + Mean..Green. + Mean..Blue. + Sex + Age + Brainbank_ID_recode, data = train_data)
# 
# # Predictions
# multinom_predictions <- predict(multinom_model, test_data)
# rf_predictions <- predict(rf_model, test_data)
# svm_predictions <- predict(svm_model, test_data)
# 
# 
# # Evaluation Metrics
# confusionMatrix(multinom_predictions, test_data$group)
# confusionMatrix(rf_predictions, test_data$group)
# confusionMatrix(svm_predictions, test_data$group)
# 
# # Assuming you have the confusion matrices stored as log_conf, rf_conf, svm_conf
# log_accuracy <- sum(diag(log_conf)) / sum(log_conf)
# rf_accuracy <- sum(diag(rf_conf)) / sum(rf_conf)
# svm_accuracy <- sum(diag(svm_conf)) / sum(svm_conf)
# 
# # Create a dataframe for plotting
# accuracy_data <- data.frame(
#   Model = c("Logistic Regression", "Random Forest", "SVM"),
#   Accuracy = c(log_accuracy, rf_accuracy, svm_accuracy)
# )
# 
# # plot
# ggplot(accuracy_data, aes(x = Model, y = Accuracy, fill = Model)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   labs(title = "Model Comparison", x = "Model", y = "Accuracy") +
#   theme_minimal() +
#   scale_fill_brewer(palette = "Set1") +
#   geom_text(aes(label = round(Accuracy, 2)), vjust = -0.5, size = 3.5)
# 


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
