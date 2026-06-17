## Load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(ggridges)

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/PhaseWY/Results/")


PhaseWY_SW <- read.delim("Larks/Genome_summary_Skylark_West.bed", sep="\t", head=F)
PhaseWY_SW$Species <- "Skylark"
PhaseWY_R <- read.delim("Larks/Genome_summary_Rasolark.bed", sep="\t", head=F)
PhaseWY_R$Species <- "Raso lark"

PhaseWY <- rbind(PhaseWY_SW, PhaseWY_R)
colnames(PhaseWY) <- c("scaffold", "start", "end", "data_type", "Species")
PhaseWY$Species <- factor(PhaseWY$Species, order=T, level=c("Skylark", "Raso lark"))
PhaseWY$data_type <- factor(PhaseWY$data_type, labels=c("No data", "Sex sequencing\ndepth difference", "Sex haplotype clustering\n&\ndepth difference", "Sex haplotype clustering", "Autosomal"),
                         levels=c("Missing data", "Sex depth difference", "Sex haplotype clustering & depth difference", "Sex haplotype clustering", "Autosomal"))


# Sex chrom structure
struct <- as.data.frame(matrix(NA, 9, 5))
colnames(struct) <- c("category", "size", "label", "label_height", "label_pos")
struct[1,] <- c("prop_PAR5", 15.9, "PAR 5", 2, 15.9/2)
struct[2,] <- c("prop_5", 36.2, "5", 2, 15.9  + 36.2/2)
struct[3,] <- c("prop_ancZ", 76.35, "Ancestral (S0, S1, S2, S3)", 2, 15.9 + 36.2 + 76.35/2)
struct[4,] <- c("prop_4A", 9.3, "4A", 2, 15.9 + 36.2 + 76.35 + 9.3/2)
struct[5,] <- c("prop_3a", 8.45, "3-a", 2, 15.9 + 36.2 + 76.35 + 9.3 + 8.45/2)
struct[6,] <- c("prop_3b", 3.55, "3-b", 2, 15.9 + 36.2 + 76.35 + 9.3 + 8.45 + 3.55/2)
struct[7,] <- c("prop_3c", 64.2, "3-c", 2, 15.9 + 36.2 + 76.35 + 9.3 + 8.45 + 3.55 + 64.2/2)
struct[8,] <- c("prop_PAR3", 26.6, "PAR 3", 2, 15.9 + 36.2 + 76.35 + 9.3 + 8.45 + 3.55 + 64.2 + 26.6/2)
struct[9,] <- c("end", NA, "", 2, max(PhaseWY$end)/1000000)
struct$size <- as.numeric(struct$size)
struct$label_pos <- as.numeric(struct$label_pos) * 1000000


struct$category <- factor(struct$category, order=T, levels = struct$category)  # Ensure factor levels
struct$label_height <- as.numeric(0.75)
struct$label_height[5] <- struct$label_height[6] - 0.5
struct$tick <- struct$label_height + 0.15
struct$tick[9] <- NA
struct$label2 <- paste(as.character(struct$size), "Mb", sep=" ")
struct$label2[5] <- "8.45 Mb"
struct$label2[9] <- NA
struct$age <- c(NA, "10.2 MY", "40.3 - 135.3 MY", "22 MY", "22 MY", "19.7 MY", "≤ 6.3 MY", NA, NA)

plot_PhaseWY <- ggplot() +
  geom_rect(data=PhaseWY, aes(xmin=start, ymin=0, xmax = end, ymax = 1, fill=data_type)) +
  scale_fill_manual(name="PhaseWY classification", values = c("No data"="#404040", "Autosomal"="#E4EAF0", "Sex haplotype clustering\n&\ndepth difference"="#f03b20", "Sex sequencing\ndepth difference"="#b30000", "Sex haplotype clustering"="#fecc5c")) +
  facet_wrap(~Species, nrow=2) +
  geom_vline(xintercept=15926250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=52100000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=128450000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=137750000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=146200000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=149736250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=213950000, color="black", linetype=2, linewidth=0.5) +
  scale_x_continuous(expand = c(0,0), limits = c(min(PhaseWY$start), max(PhaseWY$end))) +
  scale_y_continuous(limits=c(0,1), expand = c(0.0,0)) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box.margin = margin(l = 40),
    legend.key.size = unit(1.5, 'cm'),
    legend.title = element_text(size=23),
    legend.text = element_text(size=18),
    panel.spacing = unit(3, 'cm'),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank())


label_plot <- ggplot(struct) +
  geom_text(aes(x = label_pos, y = label_height - 0.00, label = label), size = 6) +
  geom_text(aes(x = label_pos, y = label_height - 0.15, label = label2), size = 5) +
  geom_text(aes(x = label_pos, y = label_height - 0.30, label = age), size = 5) +
  geom_segment(aes(y=tick, yend=1, x=label_pos, xend=label_pos)) +
  scale_x_continuous(expand = c(0,0), limits = c(min(PhaseWY$start), max(PhaseWY$end))) +
  scale_y_continuous(limits=c(-0.2,1), expand=c(0.0,0.0)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"))

plot_PhaseWY2 <- plot_PhaseWY / label_plot +
  plot_layout(heights = c(2.0, 1.3)) +   plot_annotation(theme = theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10)))

png("Figures/classification.png", width=4800, height=1500, res=300)
plot_PhaseWY2
dev.off()


######### Pipeline stats
struct2 <- struct
struct2$label_height[c(1,8)] <- struct2$label_height[2] + 0.075

plot_PhaseWY_Skylark <- ggplot() +
  geom_rect(data=subset(PhaseWY, Species=="Skylark"), aes(xmin=start, ymin=0, xmax = end, ymax = 0.5, fill=data_type)) +
  scale_fill_manual(name="PhaseWY classification", values = c("No data"="#404040", "Autosomal"="#E4EAF0", "Sex haplotype clustering\n&\ndepth difference"="#f03b20", "Sex sequencing\ndepth difference"="#b30000", "Sex haplotype clustering"="#fecc5c")) +
  geom_vline(xintercept=15926250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=52100000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=128450000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=137750000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=146200000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=149736250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=213950000, color="black", linetype=2, linewidth=0.5) +
  geom_rect(aes(xmin=0, ymin=0.5, xmax = max(PhaseWY$end), ymax = 1.5), fill="white") +
  geom_text(data=struct2, aes(x = label_pos, y = -label_height + 1.65, label = label), size = 6) +
  geom_text(data=struct2, aes(x = label_pos, y = -label_height + 1.5, label = age), size = 5) +
  geom_segment(data=struct2, aes(y=-tick+1.5, yend=0.5, x=label_pos, xend=label_pos)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,NA)) +
  scale_y_continuous(limits=c(0,1.5), expand = c(0.0,0.0)) +
  labs(y=expression(atop("Skylark", "classification"))) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box.margin = margin(l = 40),
    legend.key.size = unit(1.5, 'cm'),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15),
    panel.border = element_blank(),
    panel.spacing = unit(3, 'cm'),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    axis.title.y = element_text(size=15, angle=90, hjust=0.0),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank())

# Import data
data <- read.delim("Larks/Skylark_phase_windows.bed", head=T)
colnames(data)[c(1:13,15,16,18:22)] <- c("scaffold", "start", "end", "Classification", "Number of variants", "Proportion of heterogamets in smallest cluster",
                                         "Proportion of haplotypes in smallest cluster", "No. of individuals\nin smallest cluster",
                                         "No. of individuals in largest cluster", "Total SS", "Largest cluster SS", "Smallest cluster SS", "Between SS", "Variants per bp", "midpos", "Border change", "Phase switch", "Unknown", "Sex depth\ndifference", "Sex heterozygosity\ndifference")
data <- data[,c(1:13,15,16,18:22)]
index <- read.delim("Larks/subsetted_genome.fasta.fai", sep="\t", head=F)
colnames(index)[1:2] <- c("scaffold", "length")

# Modify data
data$scaffold <-  factor(data$scaffold, levels = index$scaffold)
data$Classification <- factor(data$Classification, order=T, labels=rev(c("Autosomal: no variation", "Autosomal:\nhomogamete\nin cluster", "Autosomal:\nboth heterogametic\nhaplotypes in cluster", "Autosomal:\nheterogametes\nmissing from cluster", "Sex-linked")),
                              levels=rev(c("Autosomal: no variation in window", "Autosomal: homogamete in cluster", "Autosomal: both haplotypes of heterogamete in cluster", "Autosomal: too few heterogametes in cluster", "Sex-linked")))
data$`No. of haplotypes\nin smallest cluster` <- round(data$`Proportion of haplotypes in smallest cluster` * (data$`No. of individuals in largest cluster` +  data$`No. of individuals\nin smallest cluster`))
data <- data[,c(1:4,10:15,19,20,21),]

# Get absolute positions
x <- 0
for(i in 1:nrow(index)) {
  indicies <- which(data$scaffold == index$scaffold[i])
  data$midpos[indicies] <-  (data$end[indicies] - data$start[indicies])/2 + data$start[indicies] + x
  x <- x + index$length[i]
}

# Label positions
scaffold_labels <- data %>%
  group_by(scaffold) %>%
  summarise(midpos = mean(midpos, na.rm = TRUE)) %>%
  arrange(midpos) 

# Sums of squares from haplotype clustering
data1 <- data %>% pivot_longer(
  cols = c("Total SS", "Between SS", "Largest cluster SS", "Smallest cluster SS"), 
  names_to = "sumsquares", 
  values_to = "Sums of squares")
data1$sumsquares <- factor(data1$sumsquares, order=T, levels=c("Total SS", "Between SS", "Largest cluster SS", "Smallest cluster SS"))


plot_SS <- ggplot() +
  geom_point(data=subset(data1, Classification=="Autosomal:\nhomogamete\nin cluster" & sumsquares == "Between SS"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
  geom_point(data=subset(data1, Classification=="Autosomal:\nboth heterogametic\nhaplotypes in cluster" & sumsquares == "Between SS"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
  geom_point(data=subset(data1, Classification== "Autosomal:\nheterogametes\nmissing from cluster" & sumsquares == "Between SSr"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
  geom_point(data=subset(data1, Classification=="Sex-linked" & sumsquares == "Between SS"), aes(x=midpos, y=`Sums of squares`, color=Classification)) +
  geom_vline(xintercept=15926250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=52100000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=128450000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=137750000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=146200000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=149736250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=213950000, color="black", linetype=2, linewidth=0.5) +
  scale_color_manual(name="Classification", values = c("Autosomal:\nhomogamete\nin cluster"="#2c7bb6",
                                                       "Autosomal:\nboth heterogametic\nhaplotypes in cluster"="#abd9e9",
                                                       "Autosomal:\nheterogametes\nmissing from cluster"="#fdae61",
                                                       "Sex-linked"="#d7191c")) +
  scale_x_continuous(expand = c(0,0), breaks = scaffold_labels$midpos, labels = scaffold_labels$scaffold) +
  guides(color = "none") +
  labs(y=expression(atop("Between clusters", "sums of squares"))) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.size = unit(2, 'cm'),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15),
    strip.text.y = element_text(size = 13),
    strip.placement = "outside", 
    axis.title.y = element_text(size=15),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12, color="black"),
    axis.text.x = element_blank())


# Number of individuals and haplotypes in smallest cluster & variant density
data2 <- data %>% pivot_longer(
  cols = c("No. of haplotypes\nin smallest cluster", "Variants per bp", "Sex heterozygosity\ndifference", "Sex depth\ndifference"),
  names_to = "Smallest cluster data",
  values_to = "values")

plot_haplotypes <- ggplot() +
  geom_point(data=subset(data2, Classification=="Autosomal:\nhomogamete\nin cluster" & `Smallest cluster data` == "No. of haplotypes\nin smallest cluster"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification=="Autosomal:\nboth heterogametic\nhaplotypes in cluster" & `Smallest cluster data` == "No. of haplotypes\nin smallest cluster"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification== "Autosomal:\nheterogametes\nmissing from cluster" & `Smallest cluster data` == "No. of haplotypes\nin smallest cluster"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification=="Sex-linked" & `Smallest cluster data` == "No. of haplotypes\nin smallest cluster"), aes(x=midpos, y=values, color=Classification)) +
  geom_vline(xintercept=15926250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=52100000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=128450000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=137750000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=146200000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=149736250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=213950000, color="black", linetype=2, linewidth=0.5) +
  scale_color_manual(name="Classification criteria", values = c("Autosomal:\nhomogamete\nin cluster"="#2c7bb6",
                                                                "Autosomal:\nboth heterogametic\nhaplotypes in cluster"="#abd9e9",
                                                                "Autosomal:\nheterogametes\nmissing from cluster"="#fdae61",
                                                                "Sex-linked"="#d7191c")) +
  scale_x_continuous(expand = c(0,0), breaks = scaffold_labels$midpos, labels = scaffold_labels$scaffold) +
  scale_y_continuous(breaks= c(2, 6, 10, 14,  18)) +
  guides(color=guide_legend(override.aes = list(size = 4))) +
  labs(y=expression(atop("No. of haplotypes", "in smallest cluster"))) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.size = unit(1.5, 'cm'),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15),
    strip.text.y = element_text(size = 13),
    strip.placement = "outside", 
    axis.title.y = element_text(size=15),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12, color="black"),
    axis.text.x = element_blank())

plot_depth <- ggplot() +
  geom_point(data=subset(data2, Classification=="Autosomal:\nhomogamete\nin cluster" & `Smallest cluster data` == "Sex depth\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification=="Autosomal:\nboth heterogametic\nhaplotypes in cluster" & `Smallest cluster data` == "Sex depth\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification== "Autosomal:\nheterogametes\nmissing from cluster" & `Smallest cluster data` == "Sex depth\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification=="Sex-linked" & `Smallest cluster data` == "Sex depth\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_vline(xintercept=15926250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=52100000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=128450000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=137750000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=146200000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=149736250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=213950000, color="black", linetype=2, linewidth=0.5) +
  scale_color_manual(name="Smallest cluster", values = c("Autosomal:\nhomogamete\nin cluster"="#2c7bb6",
                                                         "Autosomal:\nboth heterogametic\nhaplotypes in cluster"="#abd9e9",
                                                         "Autosomal:\nheterogametes\nmissing from cluster"="#fdae61",
                                                         "Sex-linked"="#d7191c")) +
  scale_x_continuous(expand = c(0,0), breaks = scaffold_labels$midpos, labels = scaffold_labels$scaffold) +
  scale_y_continuous(limits=c(0.25,1.5)) +
  guides(color = "none") +
  labs(y=expression(atop("Sex depth", "difference"))) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.size = unit(2, 'cm'),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15),
    strip.text.y = element_text(size = 13),
    strip.placement = "outside", 
    axis.title.y = element_text(size=15),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12, color="black"),
    axis.text.x = element_blank())

plot_heterozygosity <- ggplot() +
  geom_point(data=subset(data2, Classification=="Autosomal:\nhomogamete\nin cluster" & `Smallest cluster data` == "Sex heterozygosity\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification=="Autosomal:\nboth heterogametic\nhaplotypes in cluster" & `Smallest cluster data` == "Sex heterozygosity\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification== "Autosomal:\nheterogametes\nmissing from cluster" & `Smallest cluster data` == "Sex heterozygosity\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_point(data=subset(data2, Classification=="Sex-linked" & `Smallest cluster data` == "Sex heterozygosity\ndifference"), aes(x=midpos, y=values, color=Classification)) +
  geom_vline(xintercept=15926250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=52100000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=128450000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=137750000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=146200000, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=149736250, color="black", linetype=2, linewidth=0.5) +
  geom_vline(xintercept=213950000, color="black", linetype=2, linewidth=0.5) +
  scale_color_manual(name="Smallest cluster", values = c("Autosomal:\nhomogamete\nin cluster"="#2c7bb6",
                                                         "Autosomal:\nboth heterogametic\nhaplotypes in cluster"="#abd9e9",
                                                         "Autosomal:\nheterogametes\nmissing from cluster"="#fdae61",
                                                         "Sex-linked"="#d7191c")) +
  scale_x_continuous(expand = c(0,0), breaks = scaffold_labels$midpos, labels = scaffold_labels$scaffold) +
  guides(color = "none") +
  labs(y=expression(atop("Sex heterozygosity", "difference"))) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.size = unit(2, 'cm'),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12, color="black"),
    axis.text.x = element_blank())


# Import data
data <- read.delim("Larks/Skylark_sexdiff_genome_summary_subset_0.01_percent.bed", sep="\t", head=F)
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
  geom_vline(aes(xintercept = 0.75), color = "black", linetype = 2, linewidth = 1) +
  scale_x_continuous(limits = c(0.25, xlimit), expand = c(0, 0.01), breaks=seq(0, xlimit, 0.25)) +
  scale_fill_manual(values = c("Autosomal"="#E4EAF0", "Sex haplotype clustering\n&\ndepth difference"="#f03b20", "Sex sequencing\ndepth difference"="#b30000", "Sex haplotype clustering"="#fecc5c")) +
  labs(x = "Sex depth difference score", y = "Count", fill = "Classification") +
  guides(fill="none") +
  theme_bw() +
  theme(legend.position = "right",
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"))


## Get sample data
Skylark_inds <- read.table("Larks/INDS_Skylark_Europe.tsv", sep="\t", head=T)
samples <- Skylark_inds

### PCA
# Import data
dataSApc <- read.table("Larks/Skylark_phased_all_variants.eigenvec")
dataSApc$Data <- "Before: All variants"
dataSAeig <- read.table("Larks/Skylark_phased_all_variants.eigenval")
dataSZpc <- read.table("Larks/Skylark_homogametic_filtered.eigenvec")
dataSZpc$Data <- "After: Z-linked variants"
dataSZeig <- read.table("Larks/Skylark_homogametic_filtered.eigenval")

eigenvals <- as.data.frame(rbind(t(dataSAeig), t(dataSZeig)))
eigenvals <- cbind(eigenvals, c("All variants", "Z-linked variants"))
colnames(eigenvals) <- c(paste0("PC", 1:(ncol(eigenvals)-2)), "Data")

for(i in 1:nrow(eigenvals)) {
  eigenvals[i,1:(ncol(eigenvals)-2)] <- eigenvals[i,1:(ncol(eigenvals)-2)]/sum(eigenvals[i,1:(ncol(eigenvals)-2)])*100
}

PC_data <- rbind(dataSApc, dataSZpc)
PC_data <- PC_data[,-1]
colnames(PC_data) <- c("ID", paste0("PC", 1:(ncol(PC_data)-2)), "Data")
PC_data$Sex <- rep(NA, nrow(PC_data))

for(i in 1:nrow(PC_data)) {
  PC_data$Sex[i] <- samples$sex[which(samples$sample_name == PC_data$ID[i])]
}
PC_data$Sex[which(PC_data$Sex == "HOMGAM")] <- "Male"
PC_data$Sex[which(PC_data$Sex == "HETGAM")] <- "Female"
PC_data$Data <- factor(PC_data$Data, order=T, levels=c("Before: All variants", "After: Z-linked variants"))



PCA <- ggplot() +
  geom_point(data=PC_data, aes(x=PC1, y=PC2, color=Sex), position="identity", size=5) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("Female", "Male")) +
  facet_wrap(~Data, nrow=2, scales="free") +
  scale_y_continuous(expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=18, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x =  element_text(size=12, color="black"))


### Folded SFS
SFS_data1 <- read.table("Larks/Skylark_phased_all_variants.frq", head=T, row.names=NULL)
SFS_data1$Data <- "Before: All variants"
SFS_data2 <- read.table("Larks/Skylark_homogametic_filtered.frq", head=T, row.names=NULL)
SFS_data2$Data <- "After: Z-linked variants"


SFS <- rbind(SFS_data1, SFS_data2)
colnames(SFS) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQA", "FREQB", "Data")
SFS$UFSFS <- SFS |> dplyr::select(FREQA, FREQB) |> apply(1, function(z) min(z))
SFS$Data <- factor(SFS$Data, order=T, levels=c("Before: All variants", "After: Z-linked variants"))
SFS <- SFS[which(SFS$UFSFS > 0),]

W_frequency <- as.data.frame(c(10/(18*2), 10/(10 + 8*2)))
colnames(W_frequency) <- "freq"
W_frequency$Data <- c("Before: All variants", "After: Z-linked variants")
W_frequency$Data <- factor(W_frequency$Data, order=T, levels=c("Before: All variants", "After: Z-linked variants"))


SFS_plot <- ggplot() +
  geom_density(data=SFS, aes(x=UFSFS, y=after_stat(count/sum(count))), adjust=1.5, linewidth=1.5) +
  geom_vline(data=W_frequency, aes(xintercept=freq), color="#d7191c", linewidth=1, linetype=2) +
  facet_wrap(~Data, nrow=2, scales="free_y") +
  labs(x ="Minor allele frequency", y = "Proportion") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limit = c(0,0.5)) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=18, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, color="black"),
        axis.text.x =  element_text(size=12, color="black"))



### Inbreeding coefficient
samples <- read.table("Larks/INDS_Skylark_Europe.tsv", head=T)

data_Fis1 <- read.table("Larks/Skylark_phased_all_variants.het", head=T)
data_Fis1$Data <- "Before: All variants"
data_Fis2 <- read.table("Larks/Skylark_autosomal_filtered.het", head=T)
data_Fis2$Data <- "After: Autosomal variants"

InbCoef <- rbind(data_Fis1, data_Fis2)
InbCoef$Data <- factor(InbCoef$Data, order=T, levels=c("Before: All variants", "After: Autosomal variants"))
InbCoef$Sex <- rep(NA, nrow(InbCoef))

for(i in 1:nrow(InbCoef)) {
  InbCoef$Sex[i] <- samples$sex[which(samples$sample_name == InbCoef$INDV[i])]
}
InbCoef$Sex[which(InbCoef$Sex == "HOMGAM")] <- "Male"
InbCoef$Sex[which(InbCoef$Sex == "HETGAM")] <- "Female"

Fis <- ggplot() +
    geom_histogram(data=InbCoef, aes(x=F, fill=Sex), position="identity") +
    scale_fill_manual(values = c("#d7191c", "#2c7bb6"), limits=c("Female", "Male")) +
    facet_wrap(~Data, nrow=2) +
    labs(x=expression("Inbreeding Coefficient (F"[IS]*")"), y = "Count") +
    scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(limits=c(-0.1, 0.2)) +
    theme_bw() +
    theme(legend.key.size = unit(1, 'line'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=20), #change legend text font size
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.spacing = unit(1, "lines"),
          strip.background = element_rect(color="black", fill="white", linewidth=1),
          strip.text = element_text(size=20, color="black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=15, color="black"),
          axis.text.x =  element_text(size=15, color="black"))
          

row_last <- (plot_stacked | PCA | SFS_plot) + plot_layout(width=c(2,1,1))

fig <- plot_PhaseWY_Skylark / plot_haplotypes / plot_SS / plot_depth / plot_heterozygosity / row_last +
  plot_annotation(tag_levels = "A") + plot_layout(heights=c(1.5,1,1,1,1,2), guides = "collect") & theme(plot.margin = margin(10, 10, 10, 10), plot.tag = element_text(size = 20, face="bold"), legend.position = "right")
  
png("Figures/classification_info.png", width=6000, height=7200, res=300)
fig
dev.off()


png("Figures/Supplement/inbreeding_coef.png", width=2000, height=2000, res=300)
Fis
dev.off()


