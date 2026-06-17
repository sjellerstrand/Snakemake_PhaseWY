#!/usr/bin/Rscript

### Plot sex depth distribution
## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library (ggplot2)
library(ggridges)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

THRESHOLD <- as.numeric(THRESHOLD)

# Import data
data <- read.delim(DATA, sep="\t", head=F)
colnames(data) <- c("depth","data_type")
data$data_type <- factor(data$data_type, labels=c("All", "Sex sequencing\ndepth difference", "Sex haplotype clustering\n&\ndepth difference", "Sex haplotype clustering", "Autosomal"),
                                          levels=c("All", "Sex depth difference", "Sex haplotype clustering & depth difference", "Sex haplotype clustering", "Autosomal"))
data <- data[!is.na(data$data_type),]

# Plot settings
bin_width <- 0.015
xlimit <- 1.5
bin_n <- xlimit/bin_width

data$data_type <- factor(data$data_type, order=T, levels=c("All", "Sex sequencing\ndepth difference", "Sex haplotype clustering\n&\ndepth difference", "Autosomal", "Sex haplotype clustering"))

# Plot stacked distributions
plot_stacked <- ggplot(data, aes(x = depth, fill = data_type)) +
  geom_histogram(binwidth = bin_width, color = "black", position = "stack") +
  geom_vline(aes(xintercept = THRESHOLD), color = "black", linetype = 2, linewidth = 1) +
  scale_x_continuous(limits = c(0, xlimit), expand = c(0, 0.01), breaks=seq(0, xlimit, 0.25)) +
  scale_fill_manual(values = c("Autosomal"="#E4EAF0", "Sex haplotype clustering\n&\ndepth difference"="#f03b20", "Sex sequencing\ndepth difference"="#b30000", "Sex haplotype clustering"="#fecc5c")) +
  labs(x = "Depth score", y = "Count", fill = "Classification") +
  theme_bw() +
  theme(legend.position = "right",
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"))

png(paste(OUT, "_stacked.png", sep=""), width=3000, height=2000, res=300)
print(plot_stacked)
dev.off()


# Plot distributions separate
data <- rbind(data, data)
data$data_type[(1+(nrow(data)/2)):nrow(data)] <- "All"

data$data_type <- factor(data$data_type, order=T, levels=c("All", "Sex sequencing\ndepth difference", "Sex haplotype clustering\n&\ndepth difference", "Sex haplotype clustering", "Autosomal"))

plot_separate <- ggplot() +
  geom_density_ridges(data=data, aes(x=depth, y=data_type, fill=data_type), stat="binline", bins=bin_n, alpha=1, scale=1, draw_baseline = FALSE) +
  geom_vline(aes(xintercept=THRESHOLD), color="black", linetype=2, linewidth=1) +
  scale_x_continuous(limits=c(0,xlimit), expand = c(0.01, 0.01), breaks=seq(0, xlimit, by=0.25)) +
  scale_fill_manual(values = c("All"="#404040", "Autosomal"="#E4EAF0", "Sex haplotype clustering\n&\ndepth difference"="#f03b20", "Sex sequencing\ndepth difference"="#b30000", "Sex haplotype clustering"="#fecc5c")) +
  labs(x="Depth score", y="Category", fill="Classification") +
  theme_bw() +
  theme(legend.position= "none",
        plot.margin = margin(t = 1, r = 2, b = 1, l = 1, unit = "cm"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=15, colour="black"))

png(paste(OUT, "_separate.png", sep=""), width=3000, height=2000, res=300)
print(plot_separate)
dev.off()
q()
