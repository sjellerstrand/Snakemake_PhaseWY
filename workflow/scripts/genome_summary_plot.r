#!/usr/bin/Rscript

### Plot genome summary
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

# Import data
data <- read.delim(DATA, sep="\t", head=F)
index <- read.delim(INDEX, sep="\t", head=F)
colnames(data) <- c("contig", "start" , "end", "data_type")
colnames(index)[1:2] <- c("contig", "length")

# Modify data
data$contig <-  factor(data$contig, levels = index[,1])
short_contigs <- index$contig[which(index$length < MIN_LEN)]
data <- data[!(data$contig %in% short_contigs),]

if(nrow(data) == 0) {
  
    png(paste(OUT, ".png", sep=""), width=WIDTH, height=HEIGHT, res=300)
    plot.new()
    text(0.5, 0.5, paste("No contigs larger than the min_len parameter: ", MIN_LEN/1000000, "Mb. Lower this parameter in the config file. ", sep=""), cex = 5)
    dev.off()

} else {

  if(rev(str_split(DATA, '_', simplify=T))[3] == "liftover") {
  data$data_type <- factor(data$data_type, labels=c("No alignment", "No data", "Sex sequencing\ndepth difference", "Sex haplotype clustering\n&\ndepth difference", "Sex haplotype clustering", "Autosomal"),
                                           levels=c("No alignment", "Missing data", "Sex depth difference", "Sex haplotype clustering & depth difference", "Sex haplotype clustering", "Autosomal"))
  } else {
    data$data_type <- factor(data$data_type, labels=c("No data", "Sex sequencing\ndepth difference", "Sex haplotype clustering\n&\ndepth difference", "Sex haplotype clustering", "Autosomal"),
                                             levels=c("Missing data", "Sex depth difference", "Sex haplotype clustering & depth difference", "Sex haplotype clustering", "Autosomal"))
  }
  data$start <- as.numeric(data$start)
  data$end <- as.numeric(data$end)
  
  # Get absolute positions
  x <- 0
  for(i in 1:nrow(index)) {
    indicies <- which(data$contig == index$contig[i])
    data$start[indicies] <-  data$start[indicies] + x + 1
    data$end[indicies] <-  data$end[indicies] + x
    x <- x + index$length[i]
  }
  
  # Label positions
  contig_labels <- data %>%
    group_by(contig) %>%
    summarise(midpos = mean((start+end)/2, na.rm = TRUE)) %>%
    arrange(midpos)
  
  # Plot genome summary
  if(rev(str_split(DATA, '_', simplify=T))[3] == "liftover") {
  
    plot <- ggplot() +
      geom_rect(data=data, aes(xmin=start, ymin=0, xmax = end, ymax = 1, fill=data_type)) +
      scale_fill_manual(name="Classification", values = c("No alignment"="white", "No data"="#404040", "Autosomal"="#E4EAF0", "Sex haplotype clustering\n&\ndepth difference"="#f03b20", "Sex sequencing\ndepth difference"="#b30000", "Sex haplotype clustering"="#fecc5c")) +
      scale_x_continuous(expand = c(0,0), breaks = contig_labels$midpos, labels = contig_labels$contig) +
      scale_y_continuous(limits=c(0,1), expand = c(0,0)) +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=15), ##
        legend.text = element_text(size=13),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90))
        
  } else {
  
    plot <- ggplot() +
      geom_rect(data=data, aes(xmin=start, ymin=0, xmax = end, ymax = 1, fill=data_type)) +
      scale_fill_manual(name="Classification", values = c("No data"="#404040", "Autosomal"="#E4EAF0", "Sex haplotype clustering\n&\ndepth difference"="#f03b20", "Sex sequencing\ndepth difference"="#b30000", "Sex haplotype clustering"="#fecc5c")) +
      scale_x_continuous(expand = c(0,0), breaks = contig_labels$midpos, labels = contig_labels$contig) +
      scale_y_continuous(limits=c(0,1), expand = c(0,0)) +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=15), ##
        legend.text = element_text(size=13),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90))
  }
  
  
  
  png(paste(OUT, ".png", sep=""), width=WIDTH, height=HEIGHT*2, res=300)
  print(plot)
  dev.off()
  
}
q()
