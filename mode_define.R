# Scripts to determine the modalities of NM Area distribution.

library(readxl)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
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


###############################################################
### Part 2: Evalutate the number of modes in a distribution ###
###############################################################
var_interest = "total"
for (i in 1:length(var_interest)) {
  print(i)
  x <- log10(df_agg$`Area [µm²]`)
  
  # Find the modes of a KDE
  findmodes <- function(kde) {
    kde$x[which(c(kde$y[-1],NA) < kde$y & kde$y > c(NA,kde$y[-length(kde$y)]))]
  }
  
  # Compute mode trace, varying the bandwidth within a factor of 10
  m <- mean(x)
  id <- 1
  bw <- density(x)$bw * 10^seq(1,-1, length.out=101) 
  modes.lst <- lapply(bw, function(h) {
    m.new <- sort(findmodes(density(x, bw=h)))
    #  Associate each previous mode with a nearest new mode.
    if (length(m.new)==1) delta <- Inf else delta <- min(diff(m.new))/2
    d <- outer(m.new, m, function(x,y) abs(x-y))
    i <- apply(d, 2, which.min)
    g <- rep(NA_integer_, length(m.new))
    g[i] <- id[1:ncol(d)]
    # Create new ids for new modes that appear.
    k <- is.na(g)
    g[k] <- (sum(!k)+1):length(g)
    id <<- g
    m <<- m.new
    data.frame(bw=h, Mode=m.new, id=g)
  })
  X <- do.call(rbind, args=modes.lst)
  X$id <- factor(X$id)
  
  # Locate the modes at the most vertical portions of traces.
  minslope <- function(x, y) {
    f <- splinefun(x, y)
    e <- diff(range(x)) * 1e-4
    df2 <- function(x) ((f(x+e)-f(x-e)) / (2*e))^2 # Numerical derivative, squared
    v <- optimize(df2, c(min(x),max(x)))
    c(bw=v$minimum, slope=v$objective, Mode=f(v$minimum))
  }
  # Retain the desired modes.
  n.modes <- 2 # USER SELECTED: Following visual assessment
  bw.max <- max(subset(X, id==n.modes)$bw)
  modes <- sapply(1:n.modes, function(i) {
    Y <- subset(X, id==i & bw <= bw.max)
    minslope(Y$bw, Y$Mode)
  })
  #
  print(modes)
  # Plot the results.
  g1 <- ggplot(X, aes(bw, Mode)) +
    geom_line(aes(col=id), size=1.2, show.legend=FALSE) +
    geom_point(aes(bw, Mode), data=as.data.frame(t(modes)), size=3, col="Black", alpha=1/2) +
    scale_x_log10() + ylab("Mode") + xlab("Bandwidth") +
    coord_flip() + ylim(0.5,3.5) +
    ggtitle(var_interest[i])
  
  
  g2 <- ggplot(data.frame(x), aes(x, ..density..)) +
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
    geom_density(alpha=.2, fill="grey") + xlim(0.5,3.5) +
    geom_vline(data=as.data.frame(t(modes)),
               mapping=aes(xintercept=Mode), col="#D18A4e", size=1) +
    xlab("log10(`Area [µm²]`)") + ylab("Density") 
  
  arrange <- ggarrange(plotlist=list(g1,g1,g2,g2), nrow=2, ncol=2, widths = c(2,2))
  ggsave(paste0("Mode_Trace_",var_interest[i],".png"), arrange)
}

#######################################################
# Results; Identified 2 modes (k=2)
#######################################################

