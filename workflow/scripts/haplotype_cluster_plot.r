#!/usr/bin/Rscript

### Plot haplotype clustering stats
## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library (ggplot2)

options(scipen = 999)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

WIDTH <- as.numeric(WIDTH)
HEIGHT <- as.numeric(HEIGHT)
MIN_LEN <- as.numeric(MIN_LEN)*1000000

# Plot settings
ylimit <- 1.5 # Upper y axis limit for sex depth difference. Sometimes there are windows with a high value that makes the windows of interest difficult to discern.

# Import data
data <- read.delim(DATA, sep="\t", head=T)
colnames(data)[c(1:13,15,16,18:22)] <- c("contig", "start", "end", "Classification", "Number of variants", "Proportion of heterogamets in smallest cluster",
                                "Proportion of haplotypes in smallest cluster", "No. of individuals\nin smallest cluster",
                                "No. of individuals in largest cluster", "Total SS", "Largest cluster SS", "Smallest cluster SS", "Between SS", "Variants per bp", "midpos", "Border change", "Phase switch", "Unknown", "Sex depth\ndifference", "Sex heterozygosity\ndifference")
data <- data[,c(1:13,15,16,18:22)]
index <- read.delim(INDEX, sep="\t", head=F)
colnames(index)[1:2] <- c("contig", "length")

# Modify data
short_contigs <- index[which(index$length < MIN_LEN),1]
data <- data[!(data$contig %in% short_contigs),]
data$contig <-  factor(data$contig, levels = index$contig)
data$Classification <- factor(data$Classification, order=T, labels=rev(c("Autosomal: no variation", "Autosomal: homogamete\nin cluster", "Autosomal: both heterogametic\nhaplotypes in cluster", "Autosomal: heterogametes\nmissing from cluster", "Sex-linked")),
                                                            levels=rev(c("Autosomal: no variation in window", "Autosomal: homogamete in cluster", "Autosomal: both haplotypes of heterogamete in cluster", "Autosomal: too few heterogametes in cluster", "Sex-linked")))
data$`No. of haplotypes\nin smallest cluster` <- round(data$`Proportion of haplotypes in smallest cluster` * (data$`No. of individuals in largest cluster` +  data$`No. of individuals\nin smallest cluster`))
data <- data[,c(1:4,10:15,19,20,21),]

if(nrow(data) == 0) {

  png(paste(OUT, "Cluster_sums_of_squares.png", sep=""), width=WIDTH, height=HEIGHT, res=300)
  plot.new()
  text(0.5, 0.5, paste("No contigs larger than the min_len parameter: ", MIN_LEN/1000000, "Mb. Lower this parameter in the config file. ", sep=""), cex = 5)
  dev.off()
  
  png(paste(OUT, "Cluster_info.png", sep=""), width=WIDTH, height=HEIGHT, res=300)
  plot.new()
  text(0.5, 0.5, paste("No contigs larger than the min_len parameter: ", MIN_LEN/1000000, "Mb. Lower this parameter in the config file. ", sep=""), cex = 5)
  dev.off()
  
} else {

  # Get absolute positions
  x <- 0
  for(i in 1:nrow(index)) {
    indicies <- which(data$contig == index$contig[i])
    data$midpos[indicies] <-  (data$end[indicies] - data$start[indicies])/2 + data$start[indicies] + x
    x <- x + index$length[i]
  }
  
  # Label positions
  contig_labels <- data %>%
    group_by(contig) %>%
    summarise(midpos = mean(midpos, na.rm = TRUE)) %>%
    arrange(midpos) 
  
  # Sums of squares from haplotype clustering
  data1 <- data %>% pivot_longer(
    cols = c("Total SS", "Between SS", "Largest cluster SS", "Smallest cluster SS"), 
    names_to = "sumsquares", 
    values_to = "Sums of squares")
  data1$sumsquares <- factor(data1$sumsquares, order=T, levels=c("Total SS", "Between SS", "Largest cluster SS", "Smallest cluster SS"))
  
  plot_SS <- ggplot() +
    geom_point(data=subset(data1, Classification=="Autosomal: homogamete\nin cluster"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
    geom_point(data=subset(data1, Classification=="Autosomal: both heterogametic\nhaplotypes in cluster"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
    geom_point(data=subset(data1, Classification== "Autosomal: heterogametes\nmissing from cluster"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
    geom_point(data=subset(data1, Classification=="Sex-linked"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
    scale_color_manual(name="Smallest cluster", values = c("Autosomal: homogamete\nin cluster"="#2c7bb6",
                                                           "Autosomal: both heterogametic\nhaplotypes in cluster"="#abd9e9",
                                                           "Autosomal: heterogametes\nmissing from cluster"="#fdae61",
                                                           "Sex-linked"="#d7191c")) +
    facet_grid(`sumsquares`~., scales = "free", switch= "y") +
    scale_x_continuous(expand = c(0,0), breaks = contig_labels$midpos, labels = contig_labels$contig) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(2, 'cm'),
      legend.title = element_text(size=15),
      legend.text = element_text(size=13),
      strip.text.y = element_text(size = 13),
      strip.placement = "outside", 
      axis.title.y = element_text(size=15),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size=10),
      axis.text.x = element_text(size=10, angle=90))
  
  png(paste(OUT, "Cluster_sums_of_squares.png", sep=""), width=WIDTH, height=HEIGHT*6, res=300)
  print(plot_SS)
  dev.off()
  rm(data1)
  
  
  # Number of individuals and haplotypes in smallest cluster & variant density
  data2 <- data %>% pivot_longer(
    cols = c("No. of haplotypes\nin smallest cluster", "Variants per bp", "Sex heterozygosity\ndifference", "Sex depth\ndifference"),
    names_to = "Smallest cluster data",
    values_to = "values")
  
  if(length(which(data2$`Smallest cluster data` == "Sex depth\ndifference" & (data2$values > ylimit))) > 0) {
    data2 <- data2[-which(data2$`Smallest cluster data` == "Sex depth\ndifference" & (data2$values > ylimit)),]
  }
  
  plot_info <- ggplot() +
    geom_point(data=subset(data2, Classification=="Autosomal: homogamete\nin cluster"), aes(x=midpos, y=values, color=Classification)) +
    geom_point(data=subset(data2, Classification=="Autosomal: both heterogametic\nhaplotypes in cluster"), aes(x=midpos, y=values, color=Classification)) +
    geom_point(data=subset(data2, Classification== "Autosomal: heterogametes\nmissing from cluster"), aes(x=midpos, y=values, color=Classification)) +
    geom_point(data=subset(data2, Classification=="Sex-linked"), aes(x=midpos, y=values, color=Classification)) +
    scale_color_manual(name="Smallest cluster", values = c("Autosomal: homogamete\nin cluster"="#2c7bb6",
                                                           "Autosomal: both heterogametic\nhaplotypes in cluster"="#abd9e9",
                                                           "Autosomal: heterogametes\nmissing from cluster"="#fdae61",
                                                           "Sex-linked"="#d7191c")) +
    facet_grid(`Smallest cluster data`~., scales = "free", switch="y") +
    scale_x_continuous(expand = c(0,0), breaks = contig_labels$midpos, labels = contig_labels$contig) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(2, 'cm'),
      legend.title = element_text(size=15),
      legend.text = element_text(size=13),
      strip.text.y = element_text(size = 15),
      strip.placement = "outside", 
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size=10),
      axis.text.x = element_text(size=10, angle=90))
  
  png(paste(OUT, "Cluster_info.png", sep=""), width=WIDTH, height=HEIGHT*6, res=300)
  print(plot_info)
  dev.off()
  
}
q()
