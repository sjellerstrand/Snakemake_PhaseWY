## Load libraries
rm(list=ls())
library(tidyverse)
library(patchwork)
library(glmmTMB)
library(MuMIn)
library(ggeffects)
library(ape)
library(phytools)
library(gridExtra)
library(gridGraphics)
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/PhaseWY/Results/")


# Degeneration

age_gen_model <- readRDS("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/age_gen_model.RDS")
age_seq <- as.data.frame(seq(0, max(135.25)*100)/100)
colnames(age_seq) <- "cumAge"
age_seq$Strata_Age_Generations <- predict(age_gen_model, newdata=age_seq)
agefit <- lm(cumAge ~ -1 + Strata_Age_Generations + I(Strata_Age_Generations^2) + offset(rep(0, length(Strata_Age_Generations))), data=age_seq)

# Set up data
data <- read.delim("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/Skylark_2021/Results/Genes/Skylark_2021_Rasolark_2021_organised_data1.tsv", sep="\t", head=T)
data_deg <- data[which(data$Filter3=="OK" & data$Filter4=="OK" & data$Filter5=="OK"),]
data_deg$Strata2 <- data_deg$Strata
data_deg$Strata[which(data_deg$Strata == "PAR3" | data_deg$Strata == "PAR5")] <- "PAR"
data_deg$Strata <- factor(data_deg$Strata, order=T, labels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")), levels=rev(c("S0", "S1", "S2", "S3", "4A", "3a", "3b", "5", "3c", "PAR", "autosomal")))
data_deg$Species <- factor(data_deg$Species, order=T, labels=c("Skylark", "Raso lark"), levels=c("Skylark", "Rasolark"))

deg_prop <- data_deg |> count(Species, Strata, Wdegeneration)
deg_prop$prop <- rep(NA, nrow(deg_prop))

for(j in unique(deg_prop$Strata)) {
  deg_prop$prop[which(deg_prop$Strata == j & deg_prop$Species == "Raso lark")] <- deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Raso lark")] / sum(deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Raso lark")])
  deg_prop$prop[which(deg_prop$Strata == j & deg_prop$Species == "Skylark")] <- deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Skylark")] / sum(deg_prop$n[which(deg_prop$Strata == j & deg_prop$Species == "Skylark")])
}

deg_prop$Wdegeneration <- factor(deg_prop$Wdegeneration, order=T, labels=c("Functional", "Loss-of-function mutation", "Partial exon loss", "Full exon loss"), levels=c("W functional", "W loss of function", "W partially degenerated", "W degenerated"))
data_deg$bin_function <- rep(NA, nrow(data_deg))
data_deg$bin_function[which(data_deg$Wdegeneration == "W functional")] <- 0
data_deg$bin_function[which(data_deg$Wdegeneration != "W functional")] <- 1
data_deg$bin_function <- factor(data_deg$bin_function, levels=c(0,1))
data_deg3 <- data_deg[which(data_deg$Strata != "Autosomal"),]
data_deg3 <- data_deg3[which(!is.na(data_deg3$pHaplo)),]
data_deg3$logGeneLen <- log10(data_deg3$geneLengthDataBase)

# Fit model
globalmodel1 <- glmmTMB(bin_function ~
                          Strata_Age_Generations + pHaplo + logGeneLen + Species +
                          Strata_Age_Generations:pHaplo + Strata_Age_Generations:logGeneLen +
                          Strata_Age_Generations:Species +
                          pHaplo:logGeneLen + pHaplo:Species +
                          logGeneLen:Species,
                        family=binomial, data = data_deg3, na.action = "na.fail", REML=F)

options(na.action = "na.omit")
combinations1 <- dredge(global.model=globalmodel1, rank="AIC")

# Get model average
FinalModel <- model.avg(get.models(combinations1, subset = delta <= 2))

#Plot as function of age
pred_seq1a <- as.data.frame(seq(0, max(data_deg3$Strata_Age_Generations)*100)/100)
colnames(pred_seq1a) <- "Strata_Age_Generations"
pred_seq1a$pHaplo <- median(data_deg3$pHaplo)
pred_seq1a$logGeneLen <- median(data_deg3$logGeneLen)
pred_seq1a$Species <- "Skylark"
pred_seq1b <- pred_seq1a
pred_seq1b$Species <- "Raso lark"
pred_se <- predict(FinalModel, newdata=pred_seq1a, type="response", se.fit=T)
pred_seq1a$LOF <- pred_se$fit
pred_seq1a$SE <- pred_se$se.fit
pred_se <- predict(FinalModel, newdata=pred_seq1b, type="response", se.fit=T)
pred_seq1b$LOF <- pred_se$fit
pred_seq1b$SE <- pred_se$se.fit
pred_seq1_age <- pred_seq1a
pred_seq1_age$LOF <- rowMeans(cbind(pred_seq1a$LOF, pred_seq1b$LOF))
pred_seq1_age$SE <- rowMeans(cbind(pred_seq1a$SE, pred_seq1b$SE))
pred_seq1_age$lower <- pred_seq1_age$LOF - 1.96 * pred_seq1_age$SE
pred_seq1_age$upper <- pred_seq1_age$LOF + 1.96 * pred_seq1_age$SE
pred_seq1_age$lower[pred_seq1_age$lower < 0] <- 0
pred_seq1_age$upper[pred_seq1_age$upper > 1] <- 1
pred_seq1_age$Species <- "NA"

pred_seq2a <- as.data.frame(unique(data_deg3$Strata_Age_Generations))
colnames(pred_seq2a) <- c("Strata_Age_Generations")
pred_seq2a$pHaplo <- median(data_deg3$pHaplo)
pred_seq2a$logGeneLen <- median(data_deg3$logGeneLen)
pred_seq2a$Species <- "Skylark"
pred_seq2b <- pred_seq2a
pred_seq2b$Species <- "Raso lark"
pred_seq2a$LOF <- predict(FinalModel, newdata=pred_seq2a, type="response")
pred_seq2b$LOF <- predict(FinalModel, newdata=pred_seq2b, type="response")
pred_seq2_age <- pred_seq2a
pred_seq2_age$LOF <- rowMeans(cbind(pred_seq2a$LOF, pred_seq2b$LOF))
pred_seq2_age$Species <- "NA"

pred_seq2_age$Strata <- rep(NA, nrow(pred_seq2_age))
pred_seq2_age$prop <- rep(NA, nrow(pred_seq2_age))
pred_seq2_age$Strata[1] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[1])]))
pred_seq2_age$prop[1] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[1] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[2] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[2])]))
pred_seq2_age$prop[2] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[2] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[3] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[3])]))
pred_seq2_age$prop[3] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[3] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[4] <- "3-a"

pred_seq2_age$prop[4] <-sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[4] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age <- rbind(pred_seq2_age, pred_seq2_age[4,])
pred_seq2_age$Strata[5] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[5])]))
pred_seq2_age$prop[5] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[5] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[6] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[6])]))
pred_seq2_age$prop[6] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[6] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[7] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[7])]))
pred_seq2_age$prop[7] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[7] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[8] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[8])]))
pred_seq2_age$prop[8] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[8] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[9] <- unique(as.character(data_deg3$Strata[which(data_deg3$Strata_Age_Generations == pred_seq2_age$Strata_Age_Generations[9])]))
pred_seq2_age$prop[9] <- sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[9] & deg_prop$Wdegeneration != "Functional")])/2
pred_seq2_age$Strata[10] <- "4A"
pred_seq2_age$prop[10] <-sum(deg_prop$prop[which(deg_prop$Strata == pred_seq2_age$Strata[10] & deg_prop$Wdegeneration != "Functional")])/2

intercept95_age <- pred_seq1_age$Strata_Age_Generations[which(abs(pred_seq1_age$LOF-0.95) == min(abs(pred_seq1_age$LOF-0.95)))]
intercept_year <- (intercept95_age * coef(agefit)[1]) + (coef(agefit)[2] * (intercept95_age ^2))


deg_prop$Strata <- factor(deg_prop$Strata, order=T, labels=rev(c("W-S0", "W-S1", "W-S2", "W-S3", "W-4A", "W-3-a", "W-3-b", "W-5", "W-3-c", "PAR", "Autosomal")),
                          levels=rev(c("S0", "S1", "S2", "S3", "4A", "3-a", "3-b", "5", "3-c", "PAR", "Autosomal")))

plot_curve_age <- ggplot() +
  geom_ribbon(data=pred_seq1_age, aes(x=Strata_Age_Generations, ymin=LOF-SE*1.96, ymax=LOF+SE*1.96), alpha=0.1) +
  geom_line(data=pred_seq1_age, aes(x=Strata_Age_Generations, y=LOF-SE*1.96), alpha=0.1, linewidth=1) +
  geom_line(data=pred_seq1_age, aes(x=Strata_Age_Generations, y=LOF+SE*1.96), alpha=0.1, linewidth=1) +
  geom_line(data=pred_seq1_age, aes(x=Strata_Age_Generations, y=LOF, color="Probability of\nnon-functionality"), linewidth=1) +
  geom_point(data=pred_seq2_age, aes(x=Strata_Age_Generations, y=prop, color="Strata (Proportions)"), size=3.5) +
  geom_vline(aes(xintercept=intercept95_age, color="95% probability of\nnon-functionality"), linewidth=1, linetype=2) +
  scale_color_manual(values = c("black", "#2c7bb6", "#d7191c"), limits=c("Probability of\nnon-functionality", "Strata (Proportions)", "95% probability of\nnon-functionality")) +
  scale_x_continuous(limits=c(0,35), sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name="Age (million years in larks)", breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30,35)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_text(data = pred_seq2_age, aes(x=Strata_Age_Generations, y=prop, label = Strata, hjust = c(-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5), vjust = c(1.5,1.5,1.5,0.8,1.5,1.5,-0.5,1.5,1.5,1.5)), color = "Black", size=6) +
  labs(x ="Age (million generations)", y = expression(atop("Probability of", "non-functionality")), title = NULL) +
  guides(color="none") +
  theme_bw() +
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"))

plot_curve_age



### Stratum divergence

# Set up data
files <- dir("Larks/divergence_CDS_inv/")
files <- files[grep("ufboot", files)]

get_clade_distance <- function(tree, clade1, clade2) {
  d <- cophenetic.phylo(tree)
  sub_d <- d[clade1, clade2]
  mean(sub_d)
}

div_data <- as.data.frame(matrix(NA, length(files), 1002))
colnames(div_data)[1:2] <- c("Species", "Strata")

for(i in 1:length(files)) {
  div_data$Species[i] <- "Skylark"
  div_data$Strata[i] <- str_split(files[i], "\\.", simplify=T)[1]
  tree <- read.tree(paste("Larks/divergence_CDS_inv/", files[i], sep=""))
  tips <- tree[[1]]$tip.label
  clade_W  <- tips[grepl("_W$", tips)]
  clade_Z <- tips[!grepl("_W$", tips)]
  div_data[i,3:1002] <- sapply(tree, get_clade_distance, clade1 = clade_W, clade2 = clade_Z)
}

div_data <- div_data |>
  pivot_longer(cols = 3:1002, names_to = "bootstrap", values_to = "distance")
div_data$distance <- div_data$distance/2
div_data$distance_scaled <- div_data$distance/(7.16*10^-9)
div_data$Age <- rep(NA, nrow(div_data))
div_data$Age2 <- rep(NA, nrow(div_data))

div_data$Age[which(div_data$Strata == "5")] <- "3.7 MG"
div_data$Age[which(div_data$Strata == "S3")] <- "13.5 MG"
div_data$Age[which(div_data$Strata == "S2")] <- "23.1 MG"
div_data$Age[which(div_data$Strata == "S1")] <- "26.0 MG"
div_data$Age[which(div_data$Strata == "S0")] <- "32.3 MG"
div_data$Age[which(div_data$Strata == "4A")] <- "7.8 MG"
div_data$Age[which(div_data$Strata == "3a")] <- "7.8 MG"
div_data$Age[which(div_data$Strata == "3b")] <- "7.0 MG"
div_data$Age[which(div_data$Strata == "3c")] <- "< 2.3 MG"

div_data$Age2[which(div_data$Strata == "5")] <- 3.7
div_data$Age2[which(div_data$Strata == "S3")] <- 13.5
div_data$Age2[which(div_data$Strata == "S2")] <- 23.1
div_data$Age2[which(div_data$Strata == "S1")] <- 26.0
div_data$Age2[which(div_data$Strata == "S0")] <- 32.3
div_data$Age2[which(div_data$Strata == "4A")] <- 7.8
div_data$Age2[which(div_data$Strata == "3a")] <- 7.8
div_data$Age2[which(div_data$Strata == "3b")] <- 7.0
div_data$Age2[which(div_data$Strata == "3c")] <- 2.3
div_data$Strata_age <- factor(div_data$Strata, order=T, labels=c("3-c: < 2.3 MG", "5: 3.7 MG", "3-b: 7.0 MG", "3-a: 7.8 MG", "4A: 7.8 MG", "S3: 13.5 MG", "S2: 23.1 MG", "S1: 26.0 MG", "S0: 32.3 MG"), levels=c("3c", "5", "3b", "3a", "4A", "S3", "S2", "S1", "S0"))
div_data$Strata <- factor(div_data$Strata, order=T, labels=c("3-c", "5", "3-b", "3-a", "4A", "S3", "S2", "S1", "S0"), levels=c("3c", "5", "3b", "3a", "4A", "S3", "S2", "S1", "S0"))
div_data$Age <- factor(div_data$Age, order=T, levels=c("< 2.3 MG", "3.7 MG", "7.0 MG", "7.8 MG", "13.5 MG", "23.1 MG", "26.0 MG", "32.3 MG"))

pred_seq2_age2 <- pred_seq2_age[-which(pred_seq2_age$Strata == "PAR"),]
pred_seq2_age2$prop <- log10(pred_seq2_age2$prop)

pred_seq2_age2$Strata[which(pred_seq2_age$Strata == "S0")] <- "S0?"
pred_seq2_age2$Strata[which(pred_seq2_age$Strata == "S1")] <- "S1?"
div_data <- div_data[-which(div_data$Strata == "S1"),]

stratum_div_plot <- ggplot() +
  geom_boxplot(data=div_data, aes(x=Age2, y=log10(distance_scaled), group=Strata_age)) +
  geom_smooth(data=div_data, aes(x=Age2, y=log10(distance_scaled)), formula= y~poly(x,2), method="lm", color="black", linewidth=1) +
  scale_x_continuous(limits=c(0,35), sec.axis=sec_axis(LOF~ (. * coef(agefit)[1]) + (coef(agefit)[2] * (. ^2)), name="Age (million years in larks)", breaks=seq(0,140,20)), breaks=c(0,5,10,15,20,25,30,35)) +
  geom_text(data = pred_seq2_age2, aes(x=Strata_Age_Generations, y=c(5.4,6.1,5.6,5.3,6.1,5.9,6.1,4.8,5.4), label = Strata), color = "black", size=6) +
  labs(x ="Age (million generations)", y=expression(log[10]*"(Scaled divergence (generations))")) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"))
stratum_div_plot

factor <- div_data |>
  group_by(Strata, Age2) |>
  summarise(mean_distance = mean(distance_scaled, na.rm = TRUE), .groups = "drop")

factor$Age2 <- factor$Age2*1000000

factor$fact <- factor$Age2/factor$mean_distance


### Co-inheritance

### Skylark
tree_SW <- read.tree("Larks/coinheritance/Skylark_W_filtered_snps.min4.phy.treefile")
tree_SZ <- read.tree("Larks/coinheritance/Skylark_Z_filtered_snps.min4.phy.varsites.phy.treefile")
tree_SM <- read.tree("Larks/coinheritance/Skylark_M_filtered_snps.min4.phy.treefile")

# Make ultrametric
tree_SW <- chronos(tree_SW)
tree_SM <- chronos(tree_SM)
tree_SZ <- chronos(tree_SZ)

# Make W-M cophylo
tree_SWM <- cophylo(tree_SW, tree_SM)
plot(tree_SWM)
grid.newpage()
grid.echo(function() {
plot(tree_SWM, type="phylogram", fsize=1.0 ,part=0.38, link.lty="dashed", link.lwd=2, pts=FALSE, lwd=2, mar=c(0.1,0.1,2.1,0.1))

nodelabels.cophylo(tree_SWM$trees[[1]]$node.label[2:Nnode(tree_SWM$trees[[1]])],
                   2:Nnode(tree_SWM$trees[[1]])+Ntip(tree_SWM$trees[[1]]),frame="none",
                   cex=1.0,adj=c(1.3,-0.4),which="left")

nodelabels.cophylo(tree_SWM$trees[[2]]$node.label[2:Nnode(tree_SWM$trees[[2]])],
                   2:Nnode(tree_SWM$trees[[2]])+Ntip(tree_SWM$trees[[2]]),frame="none",
                   cex=1.0,adj=c(-0.3,-0.4),which="right")
mtext("W-chromosome", at=-0.5, adj=0, cex=1.5)
mtext("Mitochondrion", at=0.5, adj=1, cex=1.5)
})
grob_WM <- grid.grab()

# Make W-Z cophylo
tree_SWZ <- cophylo(tree_SW, tree_SZ)
plot(tree_SWZ)

grid.newpage()
grid.echo(function() {
plot(tree_SWZ, type="phylogram", fsize=1.0 ,part=0.38, link.lty="dashed", link.lwd=2, pts=FALSE, lwd=2, mar=c(0.1,0.1,2.1,0.1))

nodelabels.cophylo(tree_SWZ$trees[[1]]$node.label[2:Nnode(tree_SWZ$trees[[1]])],
                   2:Nnode(tree_SWZ$trees[[1]])+Ntip(tree_SWZ$trees[[1]]),frame="none",
                   cex=1.0,adj=c(1.3,-0.4),which="left")

nodelabels.cophylo(tree_SWZ$trees[[2]]$node.label[2:Nnode(tree_SWZ$trees[[2]])],
                   2:Nnode(tree_SWZ$trees[[2]])+Ntip(tree_SWZ$trees[[2]]),frame="none",
                   cex=1.0,adj=c(-0.3,-0.4),which="right")
mtext("W-chromosome", at=-0.5, adj=0, cex=1.5)
mtext("Z-chromosome", at=0.5, adj=1, cex=1.5)
})
grob_WZ <- grid.grab()

cophylo_plot <- ggplot() + 
  xlim(0, 2) + ylim(0, 1) + 
  labs(y="Co-segregation") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),   
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_blank(),
        axis.text.x = element_blank())

cophylo_plot <- cophylo_plot + annotation_custom(grob = grob_WM, xmin = -0.143, xmax = 1, ymin = 0, ymax = 1)
cophylo_plot <- cophylo_plot + annotation_custom(grob = grob_WZ, xmin = 1, xmax = 2.143, ymin = 0, ymax = 1)
cophylo_plot


### A:Z:W diversity
data_RA <- read.delim("Larks/scans/autosomal_windows_10000_steps_10000_exon_dist_20000/Rasolark_autosomal.pi_tajimas_D", sep=",", head=T)
data_RA_windows <- read.delim("Larks/scans/autosomal_windows_10000_steps_10000_exon_dist_20000/Rasolark_autosomal_windows_10000_steps_10000_exon_dist_20000.txt", sep="\t", head=T)
data_RA <- left_join(data_RA, data_RA_windows, by = "start")
data_RA$pi_abs <- data_RA$pi_all*(data_RA$sites/data_RA$N_callable_sites)
data_RA$pi_abs[which(is.na(data_RA$pi_abs))] <- 0
data_RA$Species <- "Raso lark"
data_RA$Data <- "A"

data_RZ <- read.delim("Larks/scans/Z_windows_10000_steps_10000_exon_dist_20000/Rasolark_Z.pi_tajimas_D", sep=",", head=T)
data_RZ_windows <- read.delim("Larks/scans/Z_windows_10000_steps_10000_exon_dist_20000/Rasolark_Z_windows_10000_steps_10000_exon_dist_20000.txt", sep="\t", head=T)
data_RZ <- left_join(data_RZ, data_RZ_windows, by = "start")
data_RZ$pi_abs <- data_RZ$pi_all*(data_RZ$sites/data_RZ$N_callable_sites)
data_RZ$pi_abs[which(is.na(data_RZ$pi_abs))] <- 0
data_RZ$Species <- "Raso lark"
data_RZ$Data <- "Z"

data_RW <- read.delim("Larks/scans/W_windows_10000_steps_10000_exon_dist_20000/Rasolark_W.pi_tajimas_D", sep=",", head=T)
data_RW_windows <- read.delim("Larks/scans/W_windows_10000_steps_10000_exon_dist_20000/Rasolark_W_windows_10000_steps_10000_exon_dist_20000.txt", sep="\t", head=T)
data_RW <- left_join(data_RW, data_RW_windows, by = "start")
data_RW$pi_abs <- data_RW$pi_all*(data_RW$sites/data_RW$N_callable_sites)
data_RW$pi_abs[which(is.na(data_RW$pi_abs))] <- 0
data_RW$Species <- "Raso lark"
data_RW$Data <- "W"

data_SA <- read.delim("Larks/scans/autosomal_windows_10000_steps_10000_exon_dist_20000/Skylark_autosomal.pi_tajimas_D", sep=",", head=T)
data_SA_windows <- read.delim("Larks/scans/autosomal_windows_10000_steps_10000_exon_dist_20000/Skylark_autosomal_windows_10000_steps_10000_exon_dist_20000.txt", sep="\t", head=T)
data_SA <- left_join(data_SA, data_SA_windows, by = "start")
data_SA$pi_abs <- data_SA$pi_all*(data_SA$sites/data_SA$N_callable_sites)
data_SA$pi_abs[which(is.na(data_SA$pi_abs))] <- 0
data_SA$Species <- "Skylark"
data_SA$Data <- "A"

data_SZ <- read.delim("Larks/scans/Z_windows_10000_steps_10000_exon_dist_20000/Skylark_Z.pi_tajimas_D", sep=",", head=T)
data_SZ_windows <- read.delim("Larks/scans/Z_windows_10000_steps_10000_exon_dist_20000/Skylark_Z_windows_10000_steps_10000_exon_dist_20000.txt", sep="\t", head=T)
data_SZ <- left_join(data_SZ, data_SZ_windows, by = "start")
data_SZ$pi_abs <- data_SZ$pi_all*(data_SZ$sites/data_SZ$N_callable_sites)
data_SZ$pi_abs[which(is.na(data_SZ$pi_abs))] <- 0
data_SZ$Species <- "Skylark"
data_SZ$Data <- "Z"

data_SW <- read.delim("Larks/scans/W_windows_10000_steps_10000_exon_dist_20000/Skylark_W.pi_tajimas_D", sep=",", head=T)
data_SW_windows <- read.delim("Larks/scans/W_windows_10000_steps_10000_exon_dist_20000/Skylark_W_windows_10000_steps_10000_exon_dist_20000.txt", sep="\t", head=T)
data_SW <- left_join(data_SW, data_SW_windows, by = "start")
data_SW$pi_abs <- data_SW$pi_all*(data_SW$sites/data_SW$N_callable_sites)
data_SW$pi_abs[which(is.na(data_SW$pi_abs))] <- 0
data_SW$Species <- "Skylark"
data_SW$Data <- "W"


# Perform Bootstraps
Nboot <- 1000

Ratios <- as.data.frame(matrix(NA, Nboot*4, 3))
colnames(Ratios) <- c("Ratio", "Data", "Species")

set.seed(123)
for(i in seq(1, Nboot*4, 4)) {
  # Resample
  resample_RA <- sample(data_RA$pi_abs, nrow(data_RA), replace = TRUE)
  resample_RZ <- sample(data_RZ$pi_abs, nrow(data_RZ), replace = TRUE)
  resample_RW <- sample(data_RW$pi_abs, nrow(data_RW), replace = TRUE)
  resample_SA <- sample(data_SA$pi_abs, nrow(data_SA), replace = TRUE)
  resample_SZ <- sample(data_SZ$pi_abs, nrow(data_SZ), replace = TRUE)
  resample_SW <- sample(data_SW$pi_abs, nrow(data_SW), replace = TRUE)
  
  # Record data
  Ratios$Ratio[i] <- mean(resample_RZ)/mean(resample_RA)
  Ratios$Data[i] <- "Z:A"
  Ratios$Species[i] <- "Raso lark"
  Ratios$Ratio[i+1] <- mean(resample_RW)/mean(resample_RA)
  Ratios$Data[i+1] <- "W:A"
  Ratios$Species[i+1] <- "Raso lark"
  Ratios$Ratio[i+2] <- mean(resample_SZ)/mean(resample_SA)
  Ratios$Data[i+2] <- "Z:A"
  Ratios$Species[i+2] <- "Skylark"
  Ratios$Ratio[i+3] <- mean(resample_SW)/mean(resample_SA)
  Ratios$Data[i+3] <- "W:A"
  Ratios$Species[i+3] <- "Skylark"
}


mean(data_RA$pi_abs)
mean(data_RZ$pi_abs)
mean(data_RW$pi_abs)
mean(data_SA$pi_abs)
mean(data_SZ$pi_abs)
mean(data_SW$pi_abs)

mean(Ratios$Ratio[which(Ratios$Species == "Raso lark" & Ratios$Data == "Z:A")])
mean(Ratios$Ratio[which(Ratios$Species == "Skylark" & Ratios$Data == "Z:A")])
mean(Ratios$Ratio[which(Ratios$Species == "Raso lark" & Ratios$Data == "W:A")])
mean(Ratios$Ratio[which(Ratios$Species == "Skylark" & Ratios$Data == "W:A")])

intercepts <- data.frame(Data = c("Z:A", "W:A"), intercept = c(3/4, 1/4))
intercepts$Data <- factor(intercepts$Data, order=T, levels=c("Z:A", "W:A"))
yscales <- data.frame(Data = c("Z:A", "W:A"), ymin = c(0, 0), ymax = c(1, 0.25/(0.75)))
Ratios <- Ratios |> left_join(yscales, by = "Data")
Ratios$Species <- factor(Ratios$Species, order=T, levels=c("Skylark", "Raso lark"))
Ratios$Data <- factor(Ratios$Data, order=T, levels=c("Z:A", "W:A"))

breaks_fun <- function(x) {
  if (max(x) > 0.5) {
    c(0.00, 0.25, 0.50, 0.75, 1.00)
  } else {
    c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
  }
}

ratio2 <- ggplot() +
  geom_violin(data=Ratios, aes(x=Species, y=Ratio, fill=Species), color="black", width=1, position=position_dodge(1)) +
  geom_boxplot(data=Ratios, aes(x=Species, y=Ratio, fill=Species), color="black", width=0.05, position=position_dodge(1)) +
  geom_hline(data=intercepts, aes(yintercept = intercept), linewidth=1, linetype = 2, color = "#d7191c") +
  scale_fill_manual(values = c("Raso lark" = "#92c5de", "Skylark" = "#f4a582")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(breaks = breaks_fun, limits = c(0, NA), expand = c(0,0)) + 
  geom_blank(data = Ratios, aes(y = ymin)) +
  geom_blank(data = Ratios, aes(y = ymax)) + 
  facet_wrap(~Data, scales = "free_y", nrow=1, labeller = as_labeller(c("Z:A"=" Z:PAR", "W:A"=" W:PAR"))) +
  guides(fill="none") +
  ylab("Bootstrapped π ratio") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font sizestrip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=16, color="black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=20, color="black"))

ratio2


### Sexual antagonism
genome <- read.delim("Larks/subsetted_genome.fasta.fai", sep="\t", head=F)
genome <- genome[,c(1,3)]
colnames(genome) <- c("scaffold", "CumSize")

data_S <- read.delim("Larks/scans/Z_windows_10000_steps_2500_exon_dist_NA/Skylark_Z.pairwise_populations", sep=",", head=T)
data_S_windows <- read.delim("Larks/scans/Z_windows_10000_steps_2500_exon_dist_NA/Skylark_Z_windows_10000_steps_2500_exon_dist_NA.txt", sep="\t", head=T)
data_S <- left_join(data_S, data_S_windows, by = "start")
data_S$Species <- "Skylark"
data_S$dxy_Female_Male <- data_S$dxy_Female_Male*(data_S$sites/data_S$N_callable_sites)

data_Slong <- data_S |>
  pivot_longer(cols = c(dxy_Female_Male, Fst_Female_Male), names_to="Data", values_to="value")
data_Slong$Data <- factor(data_Slong$Data, order=T, labels=c("Fst", "Dxy", "GWAS: -log10(P)"), levels=c("Fst_Female_Male", "dxy_Female_Male", "GWAS"))
data_Slong <- data_Slong[,c("mid", "value", "Data", "Species")]

gwas_S <- read.delim("Larks/Skylark_GWAS_out.assoc.txt", sep="\t", head=T)
gwas_S$Data <- "GWAS: -log10(P)"
gwas_S$Species <- "Skylark"
gwas_S$Data <- factor(gwas_S$Data, order=T, levels=c("Fst", "Dxy", "GWAS: -log10(P)"),)
gwas_S$mid <- gwas_S$ps
gwas_S$value <- -log10(gwas_S$p_wald)
colnames(gwas_S)
gwas_S <- gwas_S[,c("mid", "value", "Data", "Species")]

sex_ant_S <- rbind(data_Slong, gwas_S)
sex_ant_S <- sex_ant_S[which(sex_ant_S$mid > 15900000 & sex_ant_S$mid < 213950000),]
bonf_corr_p <- as.data.frame(matrix(c(NA, "GWAS: -log10(P)"), 1, 2))
colnames(bonf_corr_p) <- c("val", "Data")
bonf_corr_p$val <- -log10(0.05/nrow(gwas_S))

sex_ant_S[which(sex_ant_S$Data == "Fst" & sex_ant_S$value > 0.05),]
data_S[which(data_S$Fst_Female_Male > 0.05),]

sex_ant_S2 <- sex_ant_S %>%
  slice_sample(prop = 0.01)

sex_ant <- ggplot() +
  geom_hline(data=bonf_corr_p, aes(yintercept=val+2), color="white", alpha=0, linewidth=1, linetype=2) +
  geom_point(data=sex_ant_S, aes(x=mid, y=value), color="black") +
  geom_hline(data=bonf_corr_p, aes(yintercept=val), color="#d7191c", linewidth=1, linetype=2) +
  scale_x_continuous(expand = c(0,0), limits=c(0, 240507424), labels = scales::label_scientific(digits = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~Data, nrow=3, scales="free_y") +
  labs(x="Position", y="Female-Male") +
  theme_bw() +
  theme(strip.placement = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x =  element_text(size=20, color="black"))
#sex_ant

fig <- stratum_div_plot / plot_curve_age / cophylo_plot / ratio2 / sex_ant +
  plot_annotation(tag_levels = "A") + plot_layout(heights=c(1,1,1,1,1), guides = "collect") & theme(plot.margin = margin(10, 10, 10, 10), plot.tag = element_text(size = 20, face="bold"), legend.position = "right")


png("Figures/applications.png", width=6000, height=9000, res=300)
fig
dev.off()
