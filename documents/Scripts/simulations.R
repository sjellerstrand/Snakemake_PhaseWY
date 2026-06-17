## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(viridis)
library(patchwork)
library(performance)
library(gamlss)
sessionInfo()

options(scipen=999)
setwd("C:/Users/Simon JE/OneDrive - Lund University/Dokument/Simon/PhD/Projects/PhaseWY/Results")

data <- read.delim("Simulations/simulation_info.tsv", head=T)
data$Population.size2 <- factor(data$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(100, 1000, 10000, 100000, 1000000))
data <- data[-which(data$Generation == 100),]
data$Generation2 <- factor(data$Generation, order=T, labels=c("1 000", "10 000", "100 000", "1 000 000", "10 000 000"), levels=c(1000, 10000, 100000, 1000000, 10000000))
data$WhatsHap <- factor(data$WhatsHap)
data$WhatsHap  <- relevel(data$WhatsHap, ref = "OFF")
data$Simulation <- factor(data$Simulation)
data$Method <- data$Pairs


data_2 <- read.delim("Simulations/simulation_info_method.tsv", head=T)
data_2$Population.size2 <- factor(data_2$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(100, 1000, 10000, 100000, 1000000))
data_2 <- data_2[-which(data_2$Generation == 100),]
data_2$Generation2 <- factor(data_2$Generation, order=T, labels=c("1 000", "10 000", "100 000", "1 000 000", "10 000 000"), levels=c(1000, 10000, 100000, 1000000, 10000000))
data_2$Method <- factor(data_2$Method)
data_2$Method  <- relevel(data_2$Method, ref = "Python script")
data_2$Simulation <- factor(data_2$Simulation)

ne_labels <- c("100" = "bold(N[e]: '100')", "1 000" = "bold(N[e]: '1 000')", "10 000" = "bold(N[e]: '10 000')", "100 000" = "bold(N[e]: '100 000')", "1 000 000" = "bold(N[e]: '1 000 000')", "10 000 000" = "bold(N[e]: '10 000 000')")


### Test 1 of sex heterozygosity differences
data_test1 <- data[which(data$Pairs == 10 & data$WhatsHap == "ON" & data$MAC == 1 & data$Model == "Hamming" & data$Window.size == 10000),c("Heterozygosity.A", "Heterozygosity.sex", "Population.size", "Generation", "Simulation", "Population.size2", "Generation2")]

test_grid <- as.data.frame(unique(cbind(as.character(data_test1$Population.size2), as.character(data_test1$Generation2))))
colnames(test_grid) <- c("Population.size2","Generation2")
test_grid$Pval_t <- rep(NA, nrow(test_grid))
test_grid$Pval_par <- rep(NA, nrow(test_grid))

for(i in 1:nrow(test_grid)) {
  ttest_data <- data_test1[which(data_test1$Population.size2 == test_grid$Population.size2[i] & data_test1$Generation2 == test_grid$Generation2[i]),]
  ttest_result1 <- t.test(ttest_data$Heterozygosity.A, ttest_data$Heterozygosity.sex, paired=T, alternative = "less")
  test_grid$Pval_t[i] <- ttest_result1$p.value
  ttest_result2 <- wilcox.test(ttest_data$Heterozygosity.A, ttest_data$Heterozygosity.sex, paired=T, alternative = "less")
  test_grid$Pval_par[i] <- ttest_result2$p.value
}


test_grid$Pval_t_corr <- test_grid$Pval_t*nrow(test_grid)
test_grid$Pval_t_corr[which(test_grid$Pval_t_corr > 1)] <- 1
test_grid$Pval_par_corr <- test_grid$Pval_par*nrow(test_grid)
test_grid$Pval_par_corr[which(test_grid$Pval_par_corr > 1)] <- 1


test_grid$sig <- rep("No", nrow(test_grid))
test_grid$sig[which(test_grid$Pval_par_corr < 0.05)] <- "Yes"
test_grid$sig2 <- rep("", nrow(test_grid))
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.1)] <- "."
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.05)] <- "*"
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.01)] <- "**"
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.001)] <- "***"


test_grid$Population.size2 <- factor(test_grid$Population.size2, order=T, levels=c("100", "1 000", "10 000", "100 000", "1 000 000", "10 000 000"))

test_grid$Generation2 <- factor(test_grid$Generation, order=T, levels=c("1 000", "10 000", "100 000", "1 000 000", "10 000 000"))

data_test1_plot <- data_test1 %>% pivot_longer(cols = c("Heterozygosity.A", "Heterozygosity.sex"), names_to = "Region", values_to = "Heterozygosity")

pred_df1 <- data_test1_plot %>%
  group_by(Population.size2, Region) %>%
  do({model <- lm(log10(Heterozygosity) ~ log10(Generation), data = .)
  pred <- data.frame(Generation = .$Generation)
  cbind(pred, predict(model, newdata = pred, interval = "confidence"))}) %>%
  ungroup() %>%
  rename(fit = fit, lwr = lwr, upr = upr)

heterozygosity_plot <- ggplot() +
  geom_ribbon(data=pred_df1, aes(x = log10(Generation), ymin = lwr, ymax = upr, fill=Region), colour=1, alpha=0.5) +
  geom_line(data=pred_df1, aes(x = log10(Generation), y=fit, group = Region), size = 1.5, colour="black") +
  geom_line(data=pred_df1, aes(x = log10(Generation), y = fit, colour=Region), size=0.8) +
  facet_wrap(~Population.size2, nrow=1,  labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  scale_colour_manual(labels=c("Heterozygosity.A"="Pseudoautosomes","Heterozygosity.sex"="Sex chromosomes"), values=c("#E4EAF0", "#b30000"), name="Region\n") +
  scale_fill_manual(labels=c("Heterozygosity.A"="Pseudoautosomes","Heterozygosity.sex"="Sex chromosomes"), values=c("#E4EAF0", "#b30000"), name="Region\n") +
  scale_x_continuous(expand = c(0.1,0.1)) +
  labs(x = expression(log[10]*"(Generations)"), y = expression(atop("Sex heterozygosity difference", log[10]*"(Heterogametes/Homogametes)"))) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))


supp1_raw <- ggplot(data_test1_plot, aes(x = Generation2, y = log10(Heterozygosity), fill=Region)) +
  geom_boxplot() +
  facet_wrap(~Population.size2, nrow=1,  labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  scale_fill_manual(labels=c("Heterozygosity.A"="Pseudoautosomes","Heterozygosity.sex"="Sex chromosomes"), values=c("#E4EAF0", "#b30000"), name="Region\n") +
  labs(x = "Generations", y = expression(atop("Sex heterozygosity difference", log[10]*"(Heterogametes/Homogametes)"))) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black", angle=90, hjust = 1))

supp1_sig <- ggplot(test_grid, aes(x = Generation2, y = Population.size2, fill=Pval_t_corr)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") +
  geom_text(aes(label = sig2), size = 20, color = "white") +
  labs(x = "Generations", y = "Effective popualtion size", fill = "\nCorrected P-value\n") +
  theme_void() +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        title = element_text(size=25),
        axis.title.y = element_text(size=30, angle=90),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25),
        axis.text.x =  element_text(size=25))



supp1 <- (supp1_raw / supp1_sig) + plot_layout(guides = "collect") +
  plot_annotation(title="Female-Male pairs=10, Window size=10 000\n", tag_levels = 'A') & theme(plot.title = element_text(size = 30, face="bold"), plot.tag = element_text(size = 25, face="bold"), plot.tag.position = c(0.05, 1.00), legend.position = "right", legend.box.spacing = unit(1, "cm"))

png("Figures/Supplement/FigureS1.png", width=5500, height=8000, res=300)
supp1
dev.off()


### Test 2 of sex depth differences
data_test2 <- data[which(data$MAC==1 & data$Pairs==10 & data$Model == "Hamming" & data$WhatsHap == "ON" & data$Window.size == 1000),c("Dropout.A", "Dropout.sex", "Population.size", "Generation", "Simulation","Population.size2", "Generation2")]

test_grid <- as.data.frame(unique(cbind(as.character(data_test2$Population.size2), as.character(data_test2$Generation2))))
colnames(test_grid) <- c("Population.size2","Generation2")
test_grid$Pval_t <- rep(NA, nrow(test_grid))
test_grid$Pval_par <- rep(NA, nrow(test_grid))

for(i in 1:nrow(test_grid)) {
  ttest_data <- data_test2[which(data_test2$Population.size2 == test_grid$Population.size2[i] & data_test2$Generation2 == test_grid$Generation2[i]),]
  ttest_result1 <- t.test(ttest_data$Dropout.A, ttest_data$Dropout.sex, paired=T, alternative = "less")
  test_grid$Pval_t[i] <- ttest_result1$p.value
  ttest_result2 <- wilcox.test(ttest_data$Dropout.A, ttest_data$Dropout.sex, paired=T, alternative = "less")
  test_grid$Pval_par[i] <- ttest_result2$p.value
}


test_grid$Pval_t_corr <- test_grid$Pval_t*nrow(test_grid)
test_grid$Pval_t_corr[which(test_grid$Pval_t_corr > 1)] <- 1
test_grid$Pval_par_corr <- test_grid$Pval_par*nrow(test_grid)
test_grid$Pval_par_corr[which(test_grid$Pval_par_corr > 1)] <- 1
test_grid$sig2 <- rep("", nrow(test_grid))
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.1)] <- "."
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.05)] <- "*"
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.01)] <- "**"
test_grid$sig2[which(test_grid$Pval_par_corr  < 0.001)] <- "***"

test_grid$sig <- rep("No", nrow(test_grid))
test_grid$sig[which(test_grid$Pval_par_corr < 0.05)] <- "Yes"
test_grid$Population.size2 <- factor(test_grid$Population.size2, order=T, levels=c("100", "1 000", "10 000", "100 000", "1 000 000", "10 000 000"))
test_grid$Generation2 <- factor(test_grid$Generation, order=T, levels=c("1 000", "10 000", "100 000", "1 000 000", "10 000 000"))

test_grid$Generation2 <- factor(test_grid$Generation, order=T, levels=c("1 000", "10 000", "100 000", "1 000 000", "10 000 000"))


data_test2_plot <- data_test2 %>% pivot_longer(cols = c("Dropout.A", "Dropout.sex"), names_to = "Region", values_to = "Dropout")


supp2_raw <- ggplot(data_test2_plot, aes(x = Generation2, y = Dropout, fill=Region)) +
  geom_boxplot() +
  facet_wrap(~Population.size2, nrow=1, labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  scale_fill_manual(labels=c("Heterozygosity.A"="Pseudoautosomes","Heterozygosity.sex"="Sex chromosomes"), values=c("#E4EAF0", "#b30000"), name="Region\n") +
  labs(x = "Generations", y = expression(atop("Proportion sex depth difference", "(Heterogametes/Homogametes < 0.75)"))) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black", angle=90, hjust = 1))


supp2_sig <- ggplot(test_grid, aes(x = Generation2, y = Population.size2, fill=Pval_t_corr)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") +
  geom_text(aes(label = sig2), size = 20, color = "white") +
  labs(x = "Generations", y = "Effective popualtion size", fill = "\nCorrected P-value\n") +
  theme_void() +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        title = element_text(size=25),
        axis.title.y = element_text(size=30, angle=90),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25),
        axis.text.x =  element_text(size=25))



supp2 <- (supp2_raw / supp2_sig) + plot_layout(guides = "collect") +
  plot_annotation(title="Female-Male pairs=10, Window size=10 000\n", tag_levels = 'A') & theme(plot.title = element_text(size = 30, face="bold"), plot.tag = element_text(size = 25, face="bold"), plot.tag.position = c(0.05, 1.00), legend.position = "right", legend.box.spacing = unit(1, "cm"))

png("Figures/Supplement/FigureS2.png", width=5500, height=8000, res=300)
supp2
dev.off()





### Test 3 of predictors on accuracy: WhatsHap, MAC, and Window size
data_test3 <- data[which(data$Pairs==10 & data$Model == "Hamming"),c("Sensitivity.Y2", "Sensitivity.X", "Sensitivity.A", "Precision.Y2", "Precision.X", "Precision.A", "Population.size", "Generation", "WhatsHap", "MAC", "Window.size", "Simulation")]

# Transform data
hist(data_test3$Population.size, breaks=10)
data_test3$Population.size <- log10(data_test3$Population.size)
hist(data_test3$Population.size, breaks=10)

hist(data_test3$Generation, breaks=10)
data_test3$Generation <- log10(data_test3$Generation)
hist(data_test3$Generation, breaks=10)

hist(data_test3$Window.size, breaks=10)
data_test3$Window.size <- log10(data_test3$Window.size)
hist(data_test3$Window.size, breaks=10)

hist(data_test3$Sensitivity.Y2, breaks=10)
hist(data_test3$Sensitivity.X, breaks=10)
hist(data_test3$Sensitivity.A, breaks=10)
hist(data_test3$Precision.Y2, breaks=10)
hist(data_test3$Precision.X, breaks=10)
hist(data_test3$Precision.A, breaks=10)

data_test3long <- data_test3 |>
  pivot_longer(
    cols = c(Sensitivity.Y2, Sensitivity.X, Sensitivity.A, Precision.Y2, Precision.X, Precision.A),
    names_to = "Region",
    values_to = "Accuracy") |>
  filter(!is.na(Accuracy))

data_test3long$Region <- factor(data_test3long$Region)

if(length(which(dir() == "gamlss_models3.rds")) == 0) {
### Test for colinearity of all factors
  colin_model3 <- gamlss(Accuracy ~
                           Region + Population.size + Generation + WhatsHap + MAC + Window.size,
                         family = BEINF, data = data_test3long)
  check_collinearity(colin_model3)
  
  ### Test model
  
  data_test3long$Region  <- relevel(data_test3long$Region, ref = "Sensitivity.A")
  model3_SA <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + WhatsHap + MAC + Window.size +
      Population.size:Generation + Population.size:WhatsHap + Population.size:MAC + Population.size:Window.size +
      Generation:WhatsHap + Generation:MAC + Generation:Window.size +
      WhatsHap:MAC + WhatsHap:Window.size +
      MAC:Window.size),
    family = BEINF, data = data_test3long)
  sum_model3_SA <- as.data.frame(summary(model3_SA))
  
  data_test3long$Region  <- relevel(data_test3long$Region, ref = "Precision.A")
  model3_PA <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + WhatsHap + MAC + Window.size +
      Population.size:Generation + Population.size:WhatsHap + Population.size:MAC + Population.size:Window.size +
      Generation:WhatsHap + Generation:MAC + Generation:Window.size +
      WhatsHap:MAC + WhatsHap:Window.size +
      MAC:Window.size),
    family = BEINF, data = data_test3long)
  sum_model3_PA <- as.data.frame(summary(model3_PA))
  
  
  data_test3long$Region  <- relevel(data_test3long$Region, ref = "Sensitivity.X")
  model3_SX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + WhatsHap + MAC + Window.size +
      Population.size:Generation + Population.size:WhatsHap + Population.size:MAC + Population.size:Window.size +
      Generation:WhatsHap + Generation:MAC + Generation:Window.size +
      WhatsHap:MAC + WhatsHap:Window.size +
      MAC:Window.size),
    family = BEINF, data = data_test3long)
  sum_model3_SX <- as.data.frame(summary(model3_SX))
  
  data_test3long$Region  <- relevel(data_test3long$Region, ref = "Precision.X")
  model3_PX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + WhatsHap + MAC + Window.size +
      Population.size:Generation + Population.size:WhatsHap + Population.size:MAC + Population.size:Window.size +
      Generation:WhatsHap + Generation:MAC + Generation:Window.size +
      WhatsHap:MAC + WhatsHap:Window.size +
      MAC:Window.size),
    family = BEINF, data = data_test3long)
  sum_model3_PX <- as.data.frame(summary(model3_PX))
  
  data_test3long$Region  <- relevel(data_test3long$Region, ref = "Sensitivity.Y2")
  model3_SY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + WhatsHap + MAC + Window.size +
      Population.size:Generation + Population.size:WhatsHap + Population.size:MAC + Population.size:Window.size +
      Generation:WhatsHap + Generation:MAC + Generation:Window.size +
      WhatsHap:MAC + WhatsHap:Window.size +
      MAC:Window.size),
    family = BEINF, data = data_test3long)
  sum_model3_SY <- as.data.frame(summary(model3_SY))
  
  data_test3long$Region  <- relevel(data_test3long$Region, ref = "Precision.Y2")
  model3_PY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + WhatsHap + MAC + Window.size +
      Population.size:Generation + Population.size:WhatsHap + Population.size:MAC + Population.size:Window.size +
      Generation:WhatsHap + Generation:MAC + Generation:Window.size +
      WhatsHap:MAC + WhatsHap:Window.size +
      MAC:Window.size),
    family = BEINF, data = data_test3long)
  sum_model3_PY <- as.data.frame(summary(model3_PY))
  
  models_list <- list(
    model3_SA = model3_SA,
    model3_PA = model3_PA,
    model3_SX = model3_SX,
    model3_PX = model3_PX,
    model3_SY = model3_SY,
    model3_PY = model3_PY,
    sum_model3_SA = sum_model3_SA,
    sum_model3_PA = sum_model3_PA,
    sum_model3_SX = sum_model3_SX,
    sum_model3_PX = sum_model3_PX,
    sum_model3_SY = sum_model3_SY,
    sum_model3_PY = sum_model3_PY
  )

  saveRDS(models_list, file = "gamlss_models3.rds")
} else {
  models_list <- readRDS("gamlss_models3.rds")
  model3_SA <- models_list$model3_SA
  model3_PA <- models_list$model3_PA
  model3_SX <- models_list$model3_SX
  model3_PX <- models_list$model3_PX
  model3_SY <- models_list$model3_SY
  model3_PY <- models_list$model3_PY
  sum_model3_SA <- models_list$sum_model3_SA
  sum_model3_PA <- models_list$sum_model3_PA
  sum_model3_SX <- models_list$sum_model3_SX
  sum_model3_PX <- models_list$sum_model3_PX
  sum_model3_SY <- models_list$sum_model3_SY
  sum_model3_PY <- models_list$sum_model3_PY
}

plot(model3_SA)          # residuals & fitted
wp(model3_SA, ylim.all=2, xlim.all=5)

region_tables3 <- rbind(
  as.data.frame(sum_model3_SA) |> mutate(Region = "Autosomal") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model3_PA) |> mutate(Region = "Autosomal") |> mutate(Measure = "Precision"),
  as.data.frame(sum_model3_SX) |> mutate(Region = "X") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model3_PX) |> mutate(Region = "X") |> mutate(Measure = "Precision"),
  as.data.frame(sum_model3_SY) |> mutate(Region = "Y") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model3_PY) |> mutate(Region = "Y") |> mutate(Measure = "Precision"))

region_tables3 <- region_tables3[grep("Region", rownames(region_tables3), invert=T),]
region_tables3 <- region_tables3[!grepl("^X\\.Intercept\\.\\.", rownames(region_tables3)), ]

region_tables3$logOdds_Ratio <- log10(exp(region_tables3$Estimate))
region_tables3$logSE <- region_tables3$`Std. Error` * log10(exp(1))
region_tables3$logME <- 1.96 * region_tables3$logSE

region_tables3$Variable <- rep(NA, nrow(region_tables3))


region_tables3[which(region_tables3$`Pr(>|t|)` <= 0.1 & region_tables3$Region == "Y" & region_tables3$Measure == "Sensitivity"),]

region_tables3$Variable2[which(region_tables3$Region == "Autosomal" & region_tables3$Measure == "Sensitivity")] <-
  c("Intercept", "log10(Population size)", "log10(Generations)", "WhatsHap (ON)", "MAC filter", "log10(Window size)",
    "log10(Population size) × log10(Generations)","log10(Population size) × WhatsHap (ON)", "log10(Population size) × MAC filter",
    "log10(Population size) × log10(Window size)",
    "log10(Generations) x WhatsHap (ON)", "log10(Generations) x MAC filter", "log10(Generations) x log10(Window size)",
    "WhatsHap (ON) x MAC filter", "WhatsHap (ON) x log10(Window size)", "MAC filter x log10(Window size)")

region_tables3$Variable2[which(region_tables3$Region == "Autosomal" & region_tables3$Measure == "Precision")] <- region_tables3$Variable2[which(region_tables3$Region == "Autosomal" & region_tables3$Measure == "Sensitivity")]
region_tables3$Variable2[which(region_tables3$Region == "X" & region_tables3$Measure == "Sensitivity")] <- region_tables3$Variable2[which(region_tables3$Region == "Autosomal" & region_tables3$Measure == "Sensitivity")]
region_tables3$Variable2[which(region_tables3$Region == "X" & region_tables3$Measure == "Precision")] <- region_tables3$Variable2[which(region_tables3$Region == "Autosomal" & region_tables3$Measure == "Sensitivity")]
region_tables3$Variable2[which(region_tables3$Region == "Y" & region_tables3$Measure == "Sensitivity")] <- region_tables3$Variable2[which(region_tables3$Region == "Autosomal" & region_tables3$Measure == "Sensitivity")]
region_tables3$Variable2[which(region_tables3$Region == "Y" & region_tables3$Measure == "Precision")] <- region_tables3$Variable2[which(region_tables3$Region == "Autosomal" & region_tables3$Measure == "Sensitivity")]
region_tables3$sig <- rep(NA, nrow(region_tables3))
region_tables3$sig[which(region_tables3$`Pr(>|t|)` <= 0.1)] <- "."
region_tables3$sig[which(region_tables3$`Pr(>|t|)` <= 0.05)] <- "*"
region_tables3$sig[which(region_tables3$`Pr(>|t|)` <= 0.01)] <- "**"
region_tables3$sig[which(region_tables3$`Pr(>|t|)` <= 0.001)] <- "***"
region_tables3 <- region_tables3[which(region_tables3$Variable2 != "Intercept"),]
region_tables3$Region[which(region_tables3$Region == "Autosomal")] <- "Pseudoautosomal"
region_tables3$Region <- factor(region_tables3$Region, order=T, c("Pseudoautosomal", "X", "Y"))
region_tables3$Measure <- factor(region_tables3$Measure, order=T, c("Sensitivity", "Precision"))


# plot logodds
labels3 <- as.data.frame(region_tables3$Variable2[which(region_tables3$Region == "Pseudoautosomal" & region_tables3$Measure == "Sensitivity")])
colnames(labels3) <- "term"
labels3$contVar <- rev(1:nrow(labels3))
region_tables3$contVar <- labels3$contVar

plot_OR3 <- ggplot(region_tables3, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logOdds_Ratio-logME, ymax = logOdds_Ratio+logME), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=labels3$term, breaks=labels3$contVar, sec.axis = dup_axis(name=NULL)) +
  scale_y_continuous(limit=c(-2.5,2.5), expand=c(0,0)) +
  geom_text(aes(label = sig), y= -2, vjust=0.75, size = 5,) +
  facet_wrap(Region~Measure, nrow=3) +
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(x = expression(atop("","Predictors")), y = expression("log"[10]*"(Odds Ratio)")) +
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y.left = element_text(size=1, color="white"),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=10, color="black"),
        axis.ticks.y.left = element_blank())


plot_OR3


png("Figures/Supplement/FigureS3.png", width=3000, height=4500, res=300)
plot_OR3
dev.off()



### Test 4 of predictors on accuracy: Window size and Model
data_test4 <- data[which(data$Pairs==10 & data$WhatsHap == "ON" & data$MAC == 1), c("Sensitivity.Y2", "Sensitivity.X", "Sensitivity.A", "Precision.Y2", "Precision.X", "Precision.A", "Population.size", "Generation", "Window.size", "Model", "Simulation")]



# Transform data
hist(data_test4$Population.size, breaks=10)
data_test4$Population.size <- log10(data_test4$Population.size)
hist(data_test4$Population.size, breaks=10)

hist(data_test4$Generation, breaks=10)
data_test4$Generation <- log10(data_test4$Generation)
hist(data_test4$Generation, breaks=10)

hist(data_test4$Window.size, breaks=10)
data_test4$Window.size <- log10(data_test4$Window.size)
hist(data_test4$Window.size, breaks=10)


hist(data_test4$Sensitivity.Y2, breaks=10)
hist(data_test4$Sensitivity.X, breaks=10)
hist(data_test4$Sensitivity.A, breaks=10)
hist(data_test4$Precision.Y2, breaks=10)
hist(data_test4$Precision.X, breaks=10)
hist(data_test4$Precision.A, breaks=10)

data_test4long <- data_test4 |>
  pivot_longer(
    cols = c(Sensitivity.Y2, Sensitivity.X, Sensitivity.A, Precision.Y2, Precision.X, Precision.A),
    names_to = "Region",
    values_to = "Accuracy") |>
  filter(!is.na(Accuracy))

data_test4long$Region <- factor(data_test4long$Region)


if(length(which(dir() == "gamlss_models4.rds")) == 0) {
  ### Test for colinearity of all factors
  colin_model4 <- gamlss(Accuracy ~
                           Region + Population.size + Generation + Window.size + Model,
                         family = BEINF, data = data_test4long)
  check_collinearity(colin_model4)
  
  
  ### Test model
  data_test4long$Region  <- relevel(data_test4long$Region, ref = "Sensitivity.A")
  model4_SA <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Window.size + Model +
      Population.size:Generation + Population.size:Window.size + Population.size:Model +
      Generation:Window.size + Generation:Model +
      Window.size:Model),
    family = BEINF, data = data_test4long)
  sum_model4_SA <- as.data.frame(summary(model4_SA))
  
  data_test4long$Region  <- relevel(data_test4long$Region, ref = "Precision.A")
  model4_PA <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Window.size + Model +
      Population.size:Generation + Population.size:Window.size + Population.size:Model +
      Generation:Window.size + Generation:Model +
      Window.size:Model),
    family = BEINF, data = data_test4long)
  sum_model4_PA <- as.data.frame(summary(model4_PA))
  
  data_test4long$Region  <- relevel(data_test4long$Region, ref = "Sensitivity.X")
  model4_SX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Window.size + Model +
      Population.size:Generation + Population.size:Window.size + Population.size:Model +
      Generation:Window.size + Generation:Model +
      Window.size:Model),
    family = BEINF, data = data_test4long)
  sum_model4_SX <- as.data.frame(summary(model4_SX))
  
  data_test4long$Region  <- relevel(data_test4long$Region, ref = "Precision.X")
  model4_PX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Window.size + Model +
      Population.size:Generation + Population.size:Window.size + Population.size:Model +
      Generation:Window.size + Generation:Model +
      Window.size:Model),
    family = BEINF, data = data_test4long)
  sum_model4_PX <- as.data.frame(summary(model4_PX)) 
  
  data_test4long$Region  <- relevel(data_test4long$Region, ref = "Sensitivity.Y2")
  model4_SY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Window.size + Model +
      Population.size:Generation + Population.size:Window.size + Population.size:Model +
      Generation:Window.size + Generation:Model +
      Window.size:Model),
    family = BEINF, data = data_test4long)
  sum_model4_SY <- as.data.frame(summary(model4_SY))
  
  data_test4long$Region  <- relevel(data_test4long$Region, ref = "Precision.Y2")
  model4_PY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Window.size + Model +
      Population.size:Generation + Population.size:Window.size + Population.size:Model +
      Generation:Window.size + Generation:Model +
      Window.size:Model),
    family = BEINF, data = data_test4long)
  sum_model4_PY <- as.data.frame(summary(model4_PY)) 
  
  models_list <- list(
    model4_SA = model4_SA,
    model4_PA = model4_PA,
    model4_SX = model4_SX,
    model4_PX = model4_PX,
    model4_SY = model4_SY,
    model4_PY = model4_PY,
    sum_model4_SA = sum_model4_SA,
    sum_model4_PA = sum_model4_PA,
    sum_model4_SX = sum_model4_SX,
    sum_model4_PX = sum_model4_PX,
    sum_model4_SY = sum_model4_SY,
    sum_model4_PY = sum_model4_PY
  )

  saveRDS(models_list, file = "gamlss_models4.rds")
} else {
  models_list <- readRDS("gamlss_models4.rds")
  model4_SA <- models_list$model4_SA
  model4_PA <- models_list$model4_PA
  model4_SX <- models_list$model4_SX
  model4_PX <- models_list$model4_PX
  model4_SY <- models_list$model4_SY
  model4_PY <- models_list$model4_PY
  sum_model4_SA <- models_list$sum_model4_SA
  sum_model4_PA <- models_list$sum_model4_PA
  sum_model4_SX <- models_list$sum_model4_SX
  sum_model4_PX <- models_list$sum_model4_PX
  sum_model4_SY <- models_list$sum_model4_SY
  sum_model4_PY <- models_list$sum_model4_PY

}

plot(model4_SA)
wp(model4_SA, ylim.all=2)
#term.plot(model4_SA)

region_tables4 <- rbind(
  as.data.frame(sum_model4_SA) |> mutate(Region = "Autosomal") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model4_PA) |> mutate(Region = "Autosomal") |> mutate(Measure = "Precision"),
  as.data.frame(sum_model4_SX) |> mutate(Region = "X") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model4_PX) |> mutate(Region = "X") |> mutate(Measure = "Precision"),
  as.data.frame(sum_model4_SY) |> mutate(Region = "Y") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model4_PY) |> mutate(Region = "Y") |> mutate(Measure = "Precision"))

region_tables4 <- region_tables4[grep("Region", rownames(region_tables4), invert=T),]
region_tables4 <- region_tables4[!grepl("^X\\.Intercept\\.\\.", rownames(region_tables4)), ]

region_tables4$logOdds_Ratio <- log10(exp(region_tables4$Estimate))
region_tables4$logSE <- region_tables4$`Std. Error` * log10(exp(1))
region_tables4$logME <- 1.96 * region_tables4$logSE

region_tables4$Variable <- rep(NA, nrow(region_tables4))
region_tables4[which(region_tables4$`Pr(>|t|)` <= 0.1 & region_tables4$Region == "Y" & region_tables4$Measure == "Sensitivity"),]

region_tables4$Variable2 <- rep(NA, nrow(region_tables4))
region_tables4$Variable2[which(region_tables4$Region == "Autosomal" & region_tables4$Measure == "Sensitivity")]  <- c("Intercept", "log10(Population size)", "log10(Generations)", "log10(Window size)", "Model (Inverse MAF)",
                                                                                                                      "log10(Population size) × log10(Generations)", "log10(Population size) × log10(Window size)", "log10(Population size) × Model (Inverse MAF)",
                                                                                                                      "log10(Generations) x log10(Window size)", "log10(Generations) x Model (Inverse MAF)",
                                                                                                                      "log10(Window size) x Model (Inverse MAF)")

region_tables4$Variable2[which(region_tables4$Region == "Autosomal" & region_tables4$Measure == "Precision")] <- region_tables4$Variable2[which(region_tables4$Region == "Autosomal" & region_tables4$Measure == "Sensitivity")]
region_tables4$Variable2[which(region_tables4$Region == "X" & region_tables4$Measure == "Sensitivity")] <- region_tables4$Variable2[which(region_tables4$Region == "Autosomal" & region_tables4$Measure == "Sensitivity")]
region_tables4$Variable2[which(region_tables4$Region == "X" & region_tables4$Measure == "Precision")] <- region_tables4$Variable2[which(region_tables4$Region == "Autosomal" & region_tables4$Measure == "Sensitivity")]
region_tables4$Variable2[which(region_tables4$Region == "Y" & region_tables4$Measure == "Sensitivity")] <- region_tables4$Variable2[which(region_tables4$Region == "Autosomal" & region_tables4$Measure == "Sensitivity")]
region_tables4$Variable2[which(region_tables4$Region == "Y" & region_tables4$Measure == "Precision")] <- region_tables4$Variable2[which(region_tables4$Region == "Autosomal" & region_tables4$Measure == "Sensitivity")]
region_tables4$sig <- rep(NA, nrow(region_tables4))
region_tables4$sig[which(region_tables4$`Pr(>|t|)` <= 0.1)] <- "."
region_tables4$sig[which(region_tables4$`Pr(>|t|)` <= 0.05)] <- "*"
region_tables4$sig[which(region_tables4$`Pr(>|t|)` <= 0.01)] <- "**"
region_tables4$sig[which(region_tables4$`Pr(>|t|)` <= 0.001)] <- "***"
region_tables4 <- region_tables4[which(region_tables4$Variable2 != "Intercept"),]
region_tables4$Region[which(region_tables4$Region == "Autosomal")] <- "Pseudoautosomal"
region_tables4$Region <- factor(region_tables4$Region, order=T, c("Pseudoautosomal", "X", "Y"))
region_tables4$Measure <- factor(region_tables4$Measure, order=T, c("Sensitivity", "Precision"))

# plot logodds
labels4 <- as.data.frame(region_tables4$Variable2[which(region_tables4$Region == "Pseudoautosomal" & region_tables4$Measure == "Sensitivity")])
colnames(labels4) <- "term"
labels4$contVar <- rev(1:nrow(labels4))
region_tables4$contVar <- labels4$contVar

plot_OR4 <- ggplot(region_tables4, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logOdds_Ratio-logME, ymax = logOdds_Ratio+logME), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=labels4$term, breaks=labels4$contVar, sec.axis = dup_axis(name=NULL)) +
  scale_y_continuous(limit=c(-2.5,2.5), expand=c(0,0)) +
  geom_text(aes(label = sig), y= -2, vjust=0.75, size = 5) +
  facet_wrap(Region~Measure, nrow=3) +
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(x = expression(atop("","Predictors")), y = expression("log"[10]*"(Odds Ratio)")) +
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y.left = element_text(size=1, color="white"),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=10, color="black"),
        axis.ticks.y.left = element_blank())


plot_OR4


png("Figures/Supplement/FigureS4.png", width=3000, height=4500, res=300)
plot_OR4
dev.off()




# Test 5 of predictors on accuracy: Pairs and Window size
data_test5 <- data[which(data$WhatsHap == "ON" & data$MAC == 1 & data$Model == "Hamming"),c("Sensitivity.Y2", "Sensitivity.X", "Sensitivity.A", "Precision.Y2", "Precision.X", "Precision.A", "Population.size", "Generation", "Pairs", "Window.size", "Simulation")]

# Transform data
hist(data_test5$Population.size, breaks=10)
data_test5$Population.size <- log10(data_test5$Population.size)
hist(data_test5$Population.size, breaks=10)

hist(data_test5$Generation, breaks=10)
data_test5$Generation <- log10(data_test5$Generation)
hist(data_test5$Generation, breaks=10)

hist(data_test5$Window.size, breaks=10)
data_test5$Window.size <- log10(data_test5$Window.size)
hist(data_test5$Window.size, breaks=10)


hist(data_test5$Sensitivity.Y2, breaks=10)
hist(data_test5$Sensitivity.X, breaks=10)
hist(data_test5$Sensitivity.A, breaks=10)
hist(data_test5$Precision.Y2, breaks=10)
hist(data_test5$Precision.X, breaks=10)
hist(data_test5$Precision.A, breaks=10)

data_test5long <- data_test5 |>
  pivot_longer(
    cols = c(Sensitivity.Y2, Sensitivity.X, Sensitivity.A, Precision.Y2, Precision.X, Precision.A),
    names_to = "Region",
    values_to = "Accuracy") |>
  filter(!is.na(Accuracy))

data_test5long$Region <- factor(data_test5long$Region)

if(length(which(dir() == "gamlss_models5.rds")) == 0) {
  ### Test for colinearity of all factors
  colin_model5 <- gamlss(Accuracy ~
                           Region + Population.size + Generation + Pairs + Window.size,
                         family = BEINF, data = data_test5long)
  check_collinearity(colin_model5)
  
  ### Test model
  data_test5long$Region  <- relevel(data_test5long$Region, ref = "Sensitivity.A")
  model5_SA <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Pairs + Window.size +
      Population.size:Generation + Population.size:Pairs + Population.size:Window.size +
      Generation:Pairs + Generation:Window.size +
      Pairs:Window.size),
    family = BEINF, data = data_test5long)
  sum_model5_SA <- as.data.frame(summary(model5_SA))
  
  data_test5long$Region  <- relevel(data_test5long$Region, ref = "Precision.A")
  model5_PA <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Pairs + Window.size +
      Population.size:Generation + Population.size:Pairs + Population.size:Window.size +
      Generation:Pairs + Generation:Window.size +
      Pairs:Window.size),
    family = BEINF, data = data_test5long)
  sum_model5_PA <- as.data.frame(summary(model5_PA))
  
  data_test5long$Region  <- relevel(data_test5long$Region, ref = "Sensitivity.X")
  model5_SX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Pairs + Window.size +
      Population.size:Generation + Population.size:Pairs + Population.size:Window.size +
      Generation:Pairs + Generation:Window.size +
      Pairs:Window.size),
    family = BEINF, data = data_test5long)
  sum_model5_SX <- as.data.frame(summary(model5_SX))
  
  data_test5long$Region  <- relevel(data_test5long$Region, ref = "Precision.X")
  model5_PX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Pairs + Window.size +
      Population.size:Generation + Population.size:Pairs + Population.size:Window.size +
      Generation:Pairs + Generation:Window.size +
      Pairs:Window.size),
    family = BEINF, data = data_test5long)
  sum_model5_PX <- as.data.frame(summary(model5_PX))
  
  data_test5long$Region  <- relevel(data_test5long$Region, ref = "Sensitivity.Y2")
  model5_SY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Pairs + Window.size +
      Population.size:Generation + Population.size:Pairs + Population.size:Window.size +
      Generation:Pairs + Generation:Window.size +
      Pairs:Window.size),
    family = BEINF, data = data_test5long)
  sum_model5_SY <- as.data.frame(summary(model5_SY))
  
  data_test5long$Region  <- relevel(data_test5long$Region, ref = "Precision.Y2")
  model5_PY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Pairs + Window.size +
      Population.size:Generation + Population.size:Pairs + Population.size:Window.size +
      Generation:Pairs + Generation:Window.size +
      Pairs:Window.size),
    family = BEINF, data = data_test5long)
  sum_model5_PY <- as.data.frame(summary(model5_PY))
  
  models_list <- list(
    model5_SA = model5_SA,
    model5_PA = model5_PA,
    model5_SX = model5_SX,
    model5_PX = model5_PX,
    model5_SY = model5_SY,
    model5_PY = model5_PY,
    sum_model5_SA = sum_model5_SA,
    sum_model5_PA = sum_model5_PA,
    sum_model5_SX = sum_model5_SX,
    sum_model5_PX = sum_model5_PX,
    sum_model5_SY = sum_model5_SY,
    sum_model5_PY = sum_model5_PY
  )
  
  saveRDS(models_list, file = "gamlss_models5.rds")
} else {
  models_list <- readRDS("gamlss_models5.rds")
  model5_SA <- models_list$model5_SA
  model5_PA <- models_list$model5_PA
  model5_SX <- models_list$model5_SX
  model5_PX <- models_list$model5_PX
  model5_SY <- models_list$model5_SY
  model5_PY <- models_list$model5_PY
  sum_model5_SA <- models_list$sum_model5_SA
  sum_model5_PA <- models_list$sum_model5_PA
  sum_model5_SX <- models_list$sum_model5_SX
  sum_model5_PX <- models_list$sum_model5_PX
  sum_model5_SY <- models_list$sum_model5_SY
  sum_model5_PY <- models_list$sum_model5_PY
}

plot(model5_SA)
wp(model5_SA, ylim.all=2, xlim.all = 5)
#term.plot(model5_SA)


region_tables5 <- rbind(
  as.data.frame(sum_model5_SA) |> mutate(Region = "Autosomal") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model5_PA) |> mutate(Region = "Autosomal") |> mutate(Measure = "Precision"),
  as.data.frame(sum_model5_SX) |> mutate(Region = "X") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model5_PX) |> mutate(Region = "X") |> mutate(Measure = "Precision"),
  as.data.frame(sum_model5_SY) |> mutate(Region = "Y") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model5_PY) |> mutate(Region = "Y") |> mutate(Measure = "Precision"))

region_tables5 <- region_tables5[grep("Region", rownames(region_tables5), invert=T),]
region_tables5 <- region_tables5[!grepl("^X\\.Intercept\\.\\.", rownames(region_tables5)), ]

region_tables5$logOdds_Ratio <- log10(exp(region_tables5$Estimate))
region_tables5$logSE <- region_tables5$`Std. Error` * log10(exp(1))
region_tables5$logME <- 1.96 * region_tables5$logSE

region_tables5$Variable <- rep(NA, nrow(region_tables5))
region_tables5[which(region_tables5$`Pr(>|t|)` <= 0.1 & region_tables5$Region == "Y" & region_tables5$Measure == "Sensitivity"),]

region_tables5$Variable2 <- rep(NA, nrow(region_tables5))


region_tables5$Variable2 <- rep(NA, nrow(region_tables5))
region_tables5$Variable2[which(region_tables5$Region == "Autosomal" & region_tables5$Measure == "Sensitivity")] <- c("Intercept", "log10(Population size)", "log10(Generations)", "Female-Male pairs", "log10(Window size)",
                                                                                                                     "log10(Population size) × log10(Generations)", "log10(Population size) × Female-Male pairs", "log10(Population size) × log10(Window size)",
                                                                                                                     "log10(Generations) x Female-Male pairs", "log10(Generations) x log10(Window size)",
                                                                                                                     "Female-Male pairs x log10(Window size)")
region_tables5$Variable2[which(region_tables5$Region == "Autosomal" & region_tables5$Measure == "Precision")] <- region_tables5$Variable2[which(region_tables5$Region == "Autosomal" & region_tables5$Measure == "Sensitivity")]
region_tables5$Variable2[which(region_tables5$Region == "X" & region_tables5$Measure == "Sensitivity")] <- region_tables5$Variable2[which(region_tables5$Region == "Autosomal" & region_tables5$Measure == "Sensitivity")]
region_tables5$Variable2[which(region_tables5$Region == "X" & region_tables5$Measure == "Precision")] <- region_tables5$Variable2[which(region_tables5$Region == "Autosomal" & region_tables5$Measure == "Sensitivity")]
region_tables5$Variable2[which(region_tables5$Region == "Y" & region_tables5$Measure == "Sensitivity")] <- region_tables5$Variable2[which(region_tables5$Region == "Autosomal" & region_tables5$Measure == "Sensitivity")]
region_tables5$Variable2[which(region_tables5$Region == "Y" & region_tables5$Measure == "Precision")] <- region_tables5$Variable2[which(region_tables5$Region == "Autosomal" & region_tables5$Measure == "Sensitivity")]
region_tables5$sig <- rep(NA, nrow(region_tables5))
region_tables5$sig[which(region_tables5$`Pr(>|t|)` <= 0.1)] <- "."
region_tables5$sig[which(region_tables5$`Pr(>|t|)` <= 0.05)] <- "*"
region_tables5$sig[which(region_tables5$`Pr(>|t|)` <= 0.01)] <- "**"
region_tables5$sig[which(region_tables5$`Pr(>|t|)` <= 0.001)] <- "***"
region_tables5 <- region_tables5[which(region_tables5$Variable2 != "Intercept"),]
region_tables5$Region[which(region_tables5$Region == "Autosomal")] <- "Pseudoautosomal"
region_tables5$Region <- factor(region_tables5$Region, order=T, c("Pseudoautosomal", "X", "Y"))
region_tables5$Measure <- factor(region_tables5$Measure, order=T, c("Sensitivity", "Precision"))

# plot logodds
labels5 <- as.data.frame(region_tables5$Variable2[which(region_tables5$Region == "Pseudoautosomal" & region_tables5$Measure == "Sensitivity")])
colnames(labels5) <- "term"
labels5$contVar <- rev(1:nrow(labels5))
region_tables5$contVar <- labels5$contVar

plot_OR5 <- ggplot(region_tables5, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logOdds_Ratio-logME, ymax = logOdds_Ratio+logME), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=labels5$term, breaks=labels5$contVar, sec.axis = dup_axis(name=NULL)) +
  scale_y_continuous(limit=c(-2.5, 2.5), expand=c(0,0)) +
  geom_text(aes(label = sig), y= -2, vjust=0.75, size = 5) +
  facet_wrap(Region~Measure, nrow=3) +
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(x = expression(atop("","Predictors")), y = expression("log"[10]*"(Odds Ratio)")) +
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y.left = element_text(size=1, color="white"),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=10, color="black"),
        axis.ticks.y.left = element_blank())


plot_OR5


png("Figures/Supplement/FigureS5.png", width=3000, height=4500, res=300)
plot_OR5
dev.off()




### Test 6 of predictors on accuracy: Method and study design
data_test6 <- data_2[, c("Sensitivity.Y2", "Sensitivity.X", "Precision.Y2", "Precision.X", "Population.size", "Generation", "Method", "Simulation")]

# Transform data
hist(data_test6$Population.size, breaks=10)
data_test6$Population.size <- log10(data_test6$Population.size)
hist(data_test6$Population.size, breaks=10)

hist(data_test6$Generation, breaks=10)
data_test6$Generation <- log10(data_test6$Generation)
hist(data_test6$Generation, breaks=10)

hist(data_test6$Sensitivity.Y2, breaks=10)
hist(data_test6$Sensitivity.X, breaks=10)
hist(data_test6$Precision.Y2, breaks=10)
hist(data_test6$Precision.X, breaks=10)

data_test6long <- data_test6 |>
  pivot_longer(
    cols = c(Sensitivity.Y2, Sensitivity.X, Precision.Y2, Precision.X),
    names_to = "Region",
    values_to = "Accuracy") |>
  filter(!is.na(Accuracy))


if(length(which(dir() == "gamlss_models6.rds")) == 0) {
  ### Test for colinearity of all factors
  colin_model6 <- gamlss(Accuracy ~
                           Region + Population.size + Generation + Method,
                         family = BEINF, data = data_test6long)
  check_collinearity(colin_model6)
  
  
  
  ### Test model
  data_test6long$Region <- factor(data_test6long$Region, levels=c("Sensitivity.X", "Sensitivity.Y2", "Precision.Y2", "Precision.X"))
  model6_SX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Method +
      Population.size:Generation + Population.size:Method +
      Generation:Method),
    family = BEINF, data = data_test6long)
  sum_model6_SX <- as.data.frame(summary(model6_SX))
  
  data_test6long$Region <- factor(data_test6long$Region, levels=c("Precision.X", "Sensitivity.X", "Sensitivity.Y2", "Precision.Y2"))
  model6_PX <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Method +
      Population.size:Generation + Population.size:Method +
      Generation:Method),
    family = BEINF, data = data_test6long)
  sum_model6_PX <- as.data.frame(summary(model6_PX))
  
  data_test6long$Region <- factor(data_test6long$Region, levels=c("Sensitivity.Y2", "Sensitivity.X", "Precision.Y2", "Precision.X"))
  model6_SY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Method +
      Population.size:Generation + Population.size:Method +
      Generation:Method),
    family = BEINF, data = data_test6long)
  sum_model6_SY <- as.data.frame(summary(model6_SY))
  
  data_test6long$Region <- factor(data_test6long$Region, levels=c("Precision.Y2", "Sensitivity.X", "Sensitivity.Y2", "Precision.X"))
  model6_PY <- gamlss(Accuracy ~ Region * (
    Population.size + Generation + Method +
      Population.size:Generation + Population.size:Method +
      Generation:Method),
    family = BEINF, data = data_test6long)
  sum_model6_PY <- as.data.frame(summary(model6_PY))
  
  models_list <- list(
    model6_SX = model6_SX,
    model6_PX = model6_PX,
    model6_SY = model6_SY,
    model6_PY = model6_PY,
    sum_model6_SX = sum_model6_SX,
    sum_model6_PX = sum_model6_PX,
    sum_model6_SY = sum_model6_SY,
    sum_model6_PY = sum_model6_PY
  )
  
  saveRDS(models_list, file = "gamlss_models6.rds")
} else {
  models_list <- readRDS("gamlss_models6.rds")
  model6_SX <- models_list$model6_SX
  model6_PX <- models_list$model6_PX
  model6_SY <- models_list$model6_SY
  model6_PY <- models_list$model6_PY
  sum_model6_SX <- models_list$sum_model6_SX
  sum_model6_PX <- models_list$sum_model6_PX
  sum_model6_SY <- models_list$sum_model6_SY
  sum_model6_PY <- models_list$sum_model6_PY
}

plot(model6_SX)          # residuals & fitted
wp(model6_SX, ylim.all=2)
#term.plot(model6_X)


region_tables6 <- rbind(
  as.data.frame(sum_model6_SX) |> mutate(Region = "X") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model6_PX) |> mutate(Region = "X") |> mutate(Measure = "Precision"),
  as.data.frame(sum_model6_SY) |> mutate(Region = "Y") |> mutate(Measure = "Sensitivity"),
  as.data.frame(sum_model6_PY) |> mutate(Region = "Y") |> mutate(Measure = "Precision"))

region_tables6 <- region_tables6[grep("Region", rownames(region_tables6), invert=T),]
region_tables6 <- region_tables6[!grepl("^X\\.Intercept\\.\\.", rownames(region_tables6)), ]

region_tables6$logOdds_Ratio <- log10(exp(region_tables6$Estimate))
region_tables6$logSE <- region_tables6$`Std. Error` * log10(exp(1))
region_tables6$logME <- 1.96 * region_tables6$logSE

region_tables6$Variable <- rep(NA, nrow(region_tables6))
region_tables6[which(region_tables6$`Pr(>|t|)` <= 0.1 & region_tables6$Region == "Y"),]

region_tables6$Variable2 <- rep(NA, nrow(region_tables6))
region_tables6$Variable2[which(region_tables6$Region == "X" & region_tables6$Measure == "Sensitivity")] <- c("Intercept", "log10(Population size)", "log10(Generations)", "Method (1 female and 2 males)", "Method (2 females and 1 male)",
                                                                                                             "log10(Population size) × log10(Generations)", "log10(Population size) × Method (1 female and 2 males)", "log10(Population size) × Method (2 females and 1 male)",
                                                                                                             "log10(Generations) x Method (1 female and 2 males)", "log10(Generations) x Method (2 females and 1 male)")
region_tables6$Variable2[which(region_tables6$Region == "X" & region_tables6$Measure == "Precision")] <- region_tables6$Variable2[which(region_tables6$Region == "X" & region_tables6$Measure == "Sensitivity")]
region_tables6$Variable2[which(region_tables6$Region == "Y" & region_tables6$Measure == "Sensitivity")] <- region_tables6$Variable2[which(region_tables6$Region == "X" & region_tables6$Measure == "Sensitivity")]
region_tables6$Variable2[which(region_tables6$Region == "Y" & region_tables6$Measure == "Precision")] <- region_tables6$Variable2[which(region_tables6$Region == "X" & region_tables6$Measure == "Sensitivity")]
region_tables6$sig <- rep(NA, nrow(region_tables6))
region_tables6$sig[which(region_tables6$`Pr(>|t|)` <= 0.1)] <- "."
region_tables6$sig[which(region_tables6$`Pr(>|t|)` <= 0.05)] <- "*"
region_tables6$sig[which(region_tables6$`Pr(>|t|)` <= 0.01)] <- "**"
region_tables6$sig[which(region_tables6$`Pr(>|t|)` <= 0.001)] <- "***"
region_tables6 <- region_tables6[which(region_tables6$Variable2 != "Intercept"),]
region_tables6$Region <- factor(region_tables6$Region, order=T, c("X", "Y"))
region_tables6$Measure <- factor(region_tables6$Measure, order=T, c("Sensitivity", "Precision"))

# plot logodds
labels6 <- as.data.frame(region_tables6$Variable2[which(region_tables6$Region == "X" & region_tables6$Measure == "Sensitivity")])
colnames(labels6) <- "term"
labels6$contVar <- rev(1:nrow(labels6))
region_tables6$contVar <- labels6$contVar

plot_OR6 <- ggplot(region_tables6, aes(x = contVar, y = logOdds_Ratio)) +
  geom_point(size = 3, color = "Black") +  # Plot odds ratios as points
  geom_errorbar(aes(ymin = logOdds_Ratio-logME, ymax = logOdds_Ratio+logME), width = 0.2, color = "black") +  # Add CI bars
  geom_hline(yintercept = log10(1), linewidth=1, linetype = 2, color = "#d7191c") +  # Reference line at OR = 0
  scale_x_continuous(labels=labels6$term, breaks=labels6$contVar, sec.axis = dup_axis(name=NULL)) +
  scale_y_continuous(limit=c(-3.75,3.75), expand=c(0,0), breaks=c(-3, -2, -1, 0, 1, 2, 3)) +
  geom_text(aes(label = sig), y= -3.25, vjust=0.75, size = 5) +
  facet_wrap(Region~Measure) +
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(x = expression(atop("","Predictors")), y = expression("log"[10]*"(Odds Ratio)")) +
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y.left = element_text(size=1, color="white"),
        axis.text.y.right = element_text(size=12, color="black"),
        axis.text.x = element_text(size=10, color="black"),
        axis.ticks.y.left = element_blank())


plot_OR6


png("Figures/Supplement/FigureS6.png", width=3000, height=2500, res=300)
plot_OR6
dev.off()


#################################


newdata1 <- expand.grid(
  Window.size = unique(data_test3long$Window.size),
  WhatsHap = c("ON","OFF"),
  Region = c("Sensitivity.Y2", "Precision.Y2"),
  #Population.size = seq(min(data_test3long$Population.size), max(data_test3long$Population.size), length.out=100),
  Population.size = median(data_test3long$Population.size),
  Generation = seq(min(data_test3long$Generation), max(data_test3long$Generation), length.out=100),
  #Generation = unique(data_test3long$Generation),
  MAC = 1)
newdata1$pred <- predict(model3_SY, newdata = newdata1, type = "response")
newdata1$Region <- factor(newdata1$Region, labels=c("Y Sensitivity", "Y Precision"), levels=c("Sensitivity.Y2", "Precision.Y2"))

window_labels <- c("2" = "Window size: 100", "3" = "Window size: 1 000", "4" = "Window size: 10 000")

test3_plot_whats_gen_windY <- ggplot() +
  geom_line(data=newdata1, aes(x=Generation, y=pred, group=WhatsHap), color="black", size=2.8) +
  geom_line(data=newdata1, aes(x=Generation, y=pred, group=WhatsHap, colour = WhatsHap), size=2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Window.size, nrow=2, labeller = labeller(Window.size = as_labeller(window_labels))) +
  labs(x = expression(log[10]*"(Generations)"), y = expression(atop("Predicted accuracy"))) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))


newdata2 <- expand.grid(
  Window.size = max(data_test3long$Window.size),
  WhatsHap = c("ON","OFF"),
  Region = c("Sensitivity.Y2", "Precision.Y2"),
  #Population.size = seq(min(data_test3long$Population.size), max(data_test3long$Population.size), length.out=100),
  Population.size = unique(data_test3long$Population.size),
  Generation = seq(min(data_test3long$Generation), max(data_test3long$Generation), length.out=100),
  #Generation = unique(data_test3long$Generation),
  MAC = 1)
newdata2$pred <- predict(model3_SY, newdata = newdata2, type = "response")

newdata2$Region <- factor(newdata2$Region, labels=c("Y Sensitivity", "Y Precision"), levels=c("Sensitivity.Y2", "Precision.Y2"))
newdata2$Population.size2 <- factor(newdata2$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(2, 3, 4, 5, 6))

test3_plot_whats_popzise_genY <- ggplot() +
  geom_line(data=newdata2, aes(x=Generation, y=pred, group=WhatsHap), color="black", size=2.8) +
  geom_line(data=newdata2, aes(x=Generation, y=pred, group=WhatsHap, colour = WhatsHap), size=2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Population.size2, nrow=2, labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  labs(x = expression(log[10]*"(Generations)"), y = expression(atop("Predicted accuracy"))) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))


png("Figures/Supplement/test3Y_window.png", width=5500, height=3000, res=300)
test3_plot_whats_gen_windY
dev.off()


newdata2.5 <- expand.grid(
  Window.size = unique(data_test3long$Window.size),
  WhatsHap = c("ON","OFF"),
  Region = c("Sensitivity.Y2", "Precision.Y2"),
  Population.size = seq(min(data_test3long$Population.size), max(data_test3long$Population.size), length.out=100),
  Generation = median(data_test3long$Generation),
  MAC = 1)
newdata2.5$pred <- predict(model3_SY, newdata = newdata2.5, type = "response")

newdata2.5$Region <- factor(newdata2.5$Region, labels=c("Y Sensitivity", "Y Precision"), levels=c("Sensitivity.Y2", "Precision.Y2"))

test3_plot_whats_popzise_windY <- ggplot() +
  geom_line(data=newdata2.5, aes(x=Population.size, y=pred, group=WhatsHap), color="black", size=2.8) +
  geom_line(data=newdata2.5, aes(x=Population.size, y=pred, group=WhatsHap, colour = WhatsHap), size=2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Window.size, nrow=2, labeller = labeller(Window.size = as_labeller(window_labels))) +
  labs(x = expression(log[10]*"(Population Size)"), y = expression(atop("Predicted accuracy"))) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18, face="bold"), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=15, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black"))

newdata3 <- expand.grid(
  Window.size = unique(data_test3long$Window.size),
  WhatsHap = c("ON","OFF"),
  Region = c("Sensitivity.X", "Precision.X"),
  #Population.size = seq(min(data_test3long$Population.size), max(data_test3long$Population.size), length.out=100),
  Population.size = median(data_test3long$Population.size),
  Generation = seq(min(data_test3long$Generation), max(data_test3long$Generation), length.out=100),
  #Generation = unique(data_test3long$Generation),
  MAC = 1)
newdata3$pred <- predict(model3_SY, newdata = newdata3, type = "response")
newdata3$Region <- factor(newdata3$Region, labels=c("X Sensitivity", "X Precision"), levels=c("Sensitivity.X", "Precision.X"))


test3_plot_whats_gen_windX <- ggplot() +
  geom_line(data=newdata3, aes(x=Generation, y=pred, group=WhatsHap), color="black", size=2.8) +
  geom_line(data=newdata3, aes(x=Generation, y=pred, group=WhatsHap, colour = WhatsHap), size=2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Window.size, nrow=2, labeller = labeller(Window.size = as_labeller(window_labels))) +
  labs(x = expression(log[10]*"(Generations)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))


png("Figures/Supplement/test3X_window.png", width=5500, height=3000, res=300)
test3_plot_whats_gen_windX
dev.off()


newdata4 <- expand.grid(
  Window.size = max(data_test3long$Window.size),
  WhatsHap = c("ON","OFF"),
  Region = c("Sensitivity.X", "Precision.X"),
  #Population.size = seq(min(data_test3long$Population.size), max(data_test3long$Population.size), length.out=100),
  Population.size = unique(data_test3long$Population.size),
  Generation = seq(min(data_test3long$Generation), max(data_test3long$Generation), length.out=100),
  #Generation = unique(data_test3long$Generation),
  MAC = 1)
newdata4$pred <- predict(model3_SY, newdata = newdata4, type = "response")

newdata4$Region <- factor(newdata4$Region, labels=c("X Sensitivity", "X Precision"), levels=c("Sensitivity.X", "Precision.X"))
newdata4$Population.size2 <- factor(newdata4$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(2, 3, 4, 5, 6))

test3_plot_whats_popzise_genX <- ggplot() +
  geom_line(data=newdata4, aes(x=Generation, y=pred, group=WhatsHap), color="black", size=2.8) +
  geom_line(data=newdata4, aes(x=Generation, y=pred, group=WhatsHap, colour = WhatsHap), size=2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Population.size2, nrow=2, labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  labs(x = expression(log[10]*"(Generations)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))


png("Figures/Supplement/test3X_popsize.png", width=5500, height=3000, res=300)
test3_plot_whats_popzise_genX
dev.off()


newdata4.5 <- expand.grid(
  Window.size = unique(data_test3long$Window.size),
  WhatsHap = c("ON","OFF"),
  Region = c("Sensitivity.X", "Precision.X"),
  Population.size = seq(min(data_test3long$Population.size), max(data_test3long$Population.size), length.out=100),
  Generation = median(data_test3long$Generation),
  MAC = 1)
newdata4.5$pred <- predict(model3_SY, newdata = newdata4.5, type = "response")

newdata4.5$Region <- factor(newdata4.5$Region, labels=c("X Sensitivity", "X Precision"), levels=c("Sensitivity.X", "Precision.X"))

test3_plot_whats_popzise_windX <- ggplot() +
  geom_line(data=newdata4.5, aes(x=Population.size, y=pred, group=WhatsHap), color="black", size=2.8) +
  geom_line(data=newdata4.5, aes(x=Population.size, y=pred, group=WhatsHap, colour = WhatsHap), size=2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Window.size, nrow=2, labeller = labeller(Window.size = as_labeller(window_labels))) +
  labs(x = expression(log[10]*"(Population Size)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))


png("Figures/Supplement/test3X_popsize_wind.png", width=5500, height=3000, res=300)
test3_plot_whats_popzise_windX
dev.off()


newdata5 <- expand.grid(
  Window.size = unique(data_test5long$Window.size),
  Pairs = unique(data_test5long$Pairs),
  Region = c("Sensitivity.Y2", "Precision.Y2"),
  #Population.size = seq(min(data_test5long$Population.size), max(data_test5long$Population.size), length.out=100),
  Population.size = median(data_test5long$Population.size),
  Generation = seq(min(data_test5long$Generation), max(data_test5long$Generation), length.out=100))
#Generation = unique(data_test5long$Generation))
newdata5$pred <- predict(model5_SY, newdata = newdata5, type = "response")
newdata5$Region <- factor(newdata5$Region, labels=c("Y Sensitivity", "Y Precision"), levels=c("Sensitivity.Y2", "Precision.Y2"))


test5_plot_gen_wind_pairsY <- ggplot() +
  geom_line(data=newdata5, aes(x=Generation, y=pred, group=Pairs), color="black", size=2.8) +
  geom_line(data=newdata5, aes(x=Generation, y=pred, group=Pairs, colour = Pairs), size=2) +
  # scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Window.size, nrow=2, labeller = labeller(Window.size = as_labeller(window_labels))) +
  labs(x = expression(log[10]*"(Generations)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))

png("Figures/Supplement/test5Y_window.png", width=5500, height=3000, res=300)
test5_plot_gen_wind_pairsY
dev.off()


newdata6 <- expand.grid(
  Window.size = max(data_test3long$Window.size),
  Pairs = unique(data_test5long$Pairs),
  Region = c("Sensitivity.Y2", "Precision.Y2"),
  #Population.size = seq(min(data_test5long$Population.size), max(data_test5long$Population.size), length.out=100),
  Population.size = unique(data_test5long$Population.size),
  Generation = seq(min(data_test5long$Generation), max(data_test5long$Generation), length.out=100))
#Generation = unique(data_test5long$Generation))
newdata6$pred <- predict(model5_SY, newdata = newdata6, type = "response")

newdata6$Region <- factor(newdata6$Region, labels=c("Y Sensitivity", "Y Precision"), levels=c("Sensitivity.Y2", "Precision.Y2"))
newdata6$Population.size2 <- factor(newdata6$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(2, 3, 4, 5, 6))

test5_plot_popzise_gen_pairsY <- ggplot() +
  geom_line(data=newdata6, aes(x=Generation, y=pred, group=Pairs), color="black", size=2.8) +
  geom_line(data=newdata6, aes(x=Generation, y=pred, group=Pairs, colour = Pairs), size=2) +
  #  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Population.size2, nrow=2, labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  labs(x = expression(log[10]*"(Generations)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))


newdata7 <- expand.grid(
  Window.size = unique(data_test5long$Window.size),
  Pairs = unique(data_test5long$Pairs),
  Region = c("Sensitivity.X", "Precision.X"),
  #Population.size = seq(min(data_test5long$Population.size), max(data_test5long$Population.size), length.out=100),
  Population.size = median(data_test5long$Population.size),
  Generation = seq(min(data_test5long$Generation), max(data_test5long$Generation), length.out=100))
#Generation = unique(data_test5long$Generation))
newdata7$pred <- predict(model5_SY, newdata = newdata7, type = "response")
newdata7$Region <- factor(newdata7$Region, labels=c("X Sensitivity", "X Precision"), levels=c("Sensitivity.X", "Precision.X"))


test5_plot_gen_wind_pairsX <- ggplot() +
  geom_line(data=newdata7, aes(x=Generation, y=pred, group=Pairs), color="black", size=2.8) +
  geom_line(data=newdata7, aes(x=Generation, y=pred, group=Pairs, colour = Pairs), size=2) +
  # scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Window.size, nrow=2, labeller = labeller(Window.size = as_labeller(window_labels))) +
  labs(x = expression(log[10]*"(Generations)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))

png("Figures/Supplement/test5X_window.png", width=5500, height=3000, res=300)
test5_plot_gen_wind_pairsX
dev.off()

newdata8 <- expand.grid(
  Window.size = max(data_test3long$Window.size),
  Pairs = unique(data_test5long$Pairs),
  Region = c("Sensitivity.X", "Precision.X"),
  #Population.size = seq(min(data_test5long$Population.size), max(data_test5long$Population.size), length.out=100),
  Population.size = unique(data_test5long$Population.size),
  Generation = seq(min(data_test5long$Generation), max(data_test5long$Generation), length.out=100))
#Generation = unique(data_test5long$Generation))
newdata8$pred <- predict(model5_SY, newdata = newdata8, type = "response")

newdata8$Region <- factor(newdata8$Region, labels=c("X Sensitivity", "X Precision"), levels=c("Sensitivity.X", "Precision.X"))
newdata8$Population.size2 <- factor(newdata8$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(2, 3, 4, 5, 6))

test5_plot_popzise_gen_pairsX <- ggplot() +
  geom_line(data=newdata8, aes(x=Generation, y=pred, group=Pairs), color="black", size=2.8) +
  geom_line(data=newdata8, aes(x=Generation, y=pred, group=Pairs, colour = Pairs), size=2) +
  #  scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Population.size2, nrow=2, labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  labs(x = expression(log[10]*"(Generations)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))



newdata9 <- expand.grid(
  Method = unique(data_test6long$Method),
  Region = c("Sensitivity.Y2", "Precision.Y2"),
  #Population.size = seq(min(data_test5long$Population.size), max(data_test5long$Population.size), length.out=100),
  Population.size = unique(data_test5long$Population.size),
  Generation = seq(min(data_test5long$Generation), max(data_test5long$Generation), length.out=100))
#Generation = unique(data_test5long$Generation))
newdata9$pred <- predict(model6_SY, newdata = newdata9, type = "response")
newdata9$Region <- factor(newdata9$Region, labels=c("Y Sensitivity", "Y Precision"), levels=c("Sensitivity.Y2", "Precision.Y2"))

newdata9$Population.size2 <- factor(newdata9$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(2, 3, 4, 5, 6))

test6Y_plot_methodY <- ggplot() +
  geom_line(data=newdata9, aes(x=Generation, y=pred, group=Method), color="black", size=2.8) +
  geom_line(data=newdata9, aes(x=Generation, y=pred, group=Method, colour = Method), size=2) +
  # scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Population.size2, nrow=2, labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  labs(x = expression(log[10]*"(Generations)"), y = ("Predicted accuracy")) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))

png("Figures/Supplement/test6Y.png", width=5500, height=3000, res=300)
test6Y_plot_methodY
dev.off()

newdata10 <- expand.grid(
  Method = unique(data_test6long$Method),
  Region = c("Sensitivity.X", "Precision.X"),
  #Population.size = seq(min(data_test5long$Population.size), max(data_test5long$Population.size), length.out=100),
  Population.size = unique(data_test5long$Population.size),
  Generation = seq(min(data_test5long$Generation), max(data_test5long$Generation), length.out=100))
#Generation = unique(data_test5long$Generation))
newdata10$pred <- predict(model6_SY, newdata = newdata10, type = "response")
newdata10$Region <- factor(newdata10$Region, labels=c("X Sensitivity", "X Precision"), levels=c("Sensitivity.X", "Precision.X"))

newdata10$Population.size2 <- factor(newdata10$Population.size, order=T, labels=c("100", "1 000", "10 000", "100 000", "1 000 000"), levels=c(2, 3, 4, 5, 6))

test6X_plot_methodX <- ggplot() +
  geom_line(data=newdata10, aes(x=Generation, y=pred, group=Method), color="black", size=2.8) +
  geom_line(data=newdata10, aes(x=Generation, y=pred, group=Method, colour = Method), size=2) +
  # scale_color_manual(values = c("#d7191c", "#2c7bb6"), limits=c("OFF", "ON")) +
  scale_x_continuous(expand=c(0.1,0.1)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(Region~Population.size2, nrow=2, labeller = labeller(Population.size2 = as_labeller(ne_labels, label_parsed))) +
  labs(x = expression(log[10]*"(Generations)"), y = "Predicted accuracy") +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=25, face="bold"), #change legend title font size
        legend.text = element_text(size=20), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=25, color="black"),
        axis.text.x = element_text(size=25, color="black"))

png("Figures/Supplement/test6X.png", width=5500, height=3000, res=300)
test6X_plot_methodX
dev.off()


Simulations_plot <-  (heterozygosity_plot / test3_plot_whats_popzise_genY / test3_plot_whats_gen_windY / test3_plot_whats_popzise_windY) +
  plot_layout(heights = c(0.5, 1, 1, 1)) +
  plot_annotation(tag_levels = 'A') & theme(plot.title = element_text(size = 30, face="bold"), plot.tag = element_text(size = 25, face="bold"), plot.tag.position = c(0.05, 1.00), legend.position = "right", legend.box.spacing = unit(1, "cm"), plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )
png("Figures/Figure_simulations.png", width=5500, height=7000, res=300)
Simulations_plot
dev.off()




### Calculate average predictive contrasts ----

## Main terms
contrast_main <- function(model, data, region_name, var, level1, level2) {
  d <- data[data$Region == region_name, ]
  if (is.numeric(d[[var]])) {
    low  <- min(d[[var]], na.rm = TRUE)
    high <- max(d[[var]], na.rm = TRUE)
    d_low <- d
    d_high <- d
    d_low[[var]]  <- low
    d_high[[var]] <- high
    p_low <- suppressWarnings(predict(model, type = "response", newdata = d_low))
    p_high <- suppressWarnings(predict(model, type = "response", newdata = d_high))
    return(mean(p_high - p_low, na.rm = TRUE))
  }
  if(is.factor(d[[var]])) {
    levs <- levels(d[[var]])[c(level1,level2)]
    effects <- numeric(length(levs))
    names(effects) <- levs
    for(i in seq_along(levs)){
      tmp <- d
      tmp[[var]] <- factor(levs[i], levels=levs)
      p <- suppressWarnings(predict(model, type="response", newdata=tmp))
      effects[i] <- mean(p, na.rm=TRUE)
    }
    return(effects[2]-effects[1])
  }
  stop("Unsupported variable type")
}

## Interactions
contrast_interaction_numeric <- function(model, data, region_name, var1, var2) {
  d <- data[data$Region == region_name, ]
  low1  <- min(d[[var1]], na.rm = TRUE)
  high1 <- max(d[[var1]], na.rm = TRUE)
  low2  <- min(d[[var2]], na.rm = TRUE)
  high2 <- max(d[[var2]], na.rm = TRUE)
  make <- function(v1, v2) {
    tmp <- d
    tmp[[var1]] <- v1
    tmp[[var2]] <- v2
    predict(model, type = "response", newdata = tmp)
  }
  LL <- make(low1, low2)
  LH <- make(low1, high2)
  HL <- make(high1, low2)
  HH <- make(high1, high2)
  diffdiff <- (HH - HL) - (LH - LL)
  mean(diffdiff, na.rm = TRUE)
}

contrast_interaction_1fac <- function(model, data, region_name, factor_var, numeric_var, level1, level2) {
  d <- data[data$Region == region_name,]
  levs <- levels(d[[factor_var]])[c(level1,level2)]
  low  <- min(d[[numeric_var]], na.rm = TRUE)
  high <- max(d[[numeric_var]], na.rm = TRUE)
  slopes <- numeric(length(levs))
  names(slopes) <- levs
  for(i in seq_along(levs)) {
    tmp_low  <- d
    tmp_high <- d
    tmp_low[[factor_var]]  <- factor(levs[i], levels = levs)
    tmp_high[[factor_var]] <- factor(levs[i], levels = levs)
    tmp_low[[numeric_var]]  <- low
    tmp_high[[numeric_var]] <- high
    p_low  <- predict(model, type = "response", newdata = tmp_low)
    p_high <- predict(model, type = "response", newdata = tmp_high)
    slopes[i] <- mean(p_high - p_low, na.rm = TRUE)
  }
  pairwise <- c(slopes[2] - slopes[1])
  return(pairwise)
}

# Model 3
data_model <- data_test3long

predictors <- c(
  "Population.size",
  "Generation",
  "WhatsHap",
  "MAC",
  "Window.size"
)

interaction_pairs <- list(
  c("Population.size", "Generation"),
  c("Population.size", "WhatsHap"),
  c("Population.size", "MAC"),
  c("Population.size", "Window.size"),
  c("Generation", "WhatsHap"),
  c("Generation", "MAC"),
  c("Generation", "Window.size"),
  c("WhatsHap", "MAC"),
  c("WhatsHap", "Window.size"),
  c("MAC", "Window.size")
)

regions <- c(
  "Sensitivity.A",
  "Precision.A",
  "Sensitivity.X",
  "Precision.X",
  "Sensitivity.Y2",
  "Precision.Y2"
)

models <- list(
  "Sensitivity.A"  = model3_SA,
  "Precision.A"    = model3_PA,
  "Sensitivity.X"  = model3_SX,
  "Precision.X"    = model3_PX,
  "Sensitivity.Y2" = model3_SY,
  "Precision.Y2"   = model3_PY
)

effect_results <- data.frame()
for (reg in regions) {
  cat("Processing:", reg, "\n")
  # MAIN EFFECTS
  for (var in predictors) {
    if(is.factor(data_model[[var]])) {
      lv <- levels(data_model[[var]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_main(models[[reg]], data_model, reg, var, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var, ": ", lv[lv1], " -> ", lv[lv2], sep=""), Effect = val))
          }
        }
      }
    } else {
     val <- contrast_main(models[[reg]], data_model, reg, var, NA, NA)
     effect_results <- rbind(effect_results, data.frame(Region = reg, Term = var, Effect = val))
    }
  }
  # INTERACTIONS
  for (pair in interaction_pairs) {
    var1 <- pair[1]
    var2 <- pair[2]
    if (is.numeric(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      val <- contrast_interaction_numeric(models[[reg]], data_model, reg, var1, var2)
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    } else if(is.factor(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      lv <- levels(data_model[[var1]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var1, var2, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, ": ", lv[lv1], " -> ", lv[lv2], " × ", var2, sep=""), Effect = val))
          }
        }
      }
    } else if(is.factor(data_model[[var2]]) && is.numeric(data_model[[var1]])) {
      lv <- levels(data_model[[var2]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var2, var1, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var2, ": ", lv[lv1], " -> ", lv[lv2], " × ", var1, sep=""), Effect = val))
          }
        }
      }
    } else {
      val <- NA
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    }
  }
}
effect_results_3 <- effect_results


# Model 4
predictors <- c(
  "Population.size",
  "Generation",
  "Window.size",
  "Model"
)

interaction_pairs <- list(
  c("Population.size", "Generation"),
  c("Population.size", "Window.size"),
  c("Population.size", "Model"),
  c("Generation", "Window.size"),
  c("Generation", "Model"),
  c("Window.size", "Model")
)

regions <- c(
  "Sensitivity.A",
  "Precision.A",
  "Sensitivity.X",
  "Precision.X",
  "Sensitivity.Y2",
  "Precision.Y2"
)

models <- list(
  "Sensitivity.A"  = model4_SA,
  "Precision.A"    = model4_PA,
  "Sensitivity.X"  = model4_SX,
  "Precision.X"    = model4_PX,
  "Sensitivity.Y2" = model4_SY,
  "Precision.Y2"   = model4_PY
)

data_model <- data_test4long
data_model$Model <- factor(data_model$Model)
effect_results <- data.frame()
for (reg in regions) {
  cat("Processing:", reg, "\n")
  # MAIN EFFECTS
  for (var in predictors) {
    if(is.factor(data_model[[var]])) {
      lv <- levels(data_model[[var]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_main(models[[reg]], data_model, reg, var, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var, ": ", lv[lv1], " -> ", lv[lv2], sep=""), Effect = val))
          }
        }
      }
    } else {
      val <- contrast_main(models[[reg]], data_model, reg, var, NA, NA)
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = var, Effect = val))
    }
  }
  # INTERACTIONS
  for (pair in interaction_pairs) {
    var1 <- pair[1]
    var2 <- pair[2]
    if (is.numeric(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      val <- contrast_interaction_numeric(models[[reg]], data_model, reg, var1, var2)
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    } else if(is.factor(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      lv <- levels(data_model[[var1]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var1, var2, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, ": ", lv[lv1], " -> ", lv[lv2], " × ", var2, sep=""), Effect = val))
          }
        }
      }
    } else if(is.factor(data_model[[var2]]) && is.numeric(data_model[[var1]])) {
      lv <- levels(data_model[[var2]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var2, var1, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var2, ": ", lv[lv1], " -> ", lv[lv2], " × ", var1, sep=""), Effect = val))
          }
        }
      }
    } else {
      val <- NA
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    }
  }
}
effect_results_4 <- effect_results


# Model 5

predictors <- c(
  "Population.size",
  "Generation",
  "Pairs",
  "Window.size"
)

interaction_pairs <- list(
  c("Population.size", "Generation"),
  c("Population.size", "Pairs"),
  c("Population.size", "Window.size"),
  c("Generation", "Pairs"),
  c("Generation", "Window.size"),
  c("Pairs", "Window.size")
)

regions <- c(
  "Sensitivity.A",
  "Precision.A",
  "Sensitivity.X",
  "Precision.X",
  "Sensitivity.Y2",
  "Precision.Y2"
)

models <- list(
  "Sensitivity.A"  = model5_SA,
  "Precision.A"    = model5_PA,
  "Sensitivity.X"  = model5_SX,
  "Precision.X"    = model5_PX,
  "Sensitivity.Y2" = model5_SY,
  "Precision.Y2"   = model5_PY
)

data_model <- data_test5long

effect_results <- data.frame()
for (reg in regions) {
  cat("Processing:", reg, "\n")
  # MAIN EFFECTS
  for (var in predictors) {
    if(is.factor(data_model[[var]])) {
      lv <- levels(data_model[[var]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_main(models[[reg]], data_model, reg, var, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var, ": ", lv[lv1], " -> ", lv[lv2], sep=""), Effect = val))
          }
        }
      }
    } else {
      val <- contrast_main(models[[reg]], data_model, reg, var, NA, NA)
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = var, Effect = val))
    }
  }
  # INTERACTIONS
  for (pair in interaction_pairs) {
    var1 <- pair[1]
    var2 <- pair[2]
    if (is.numeric(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      val <- contrast_interaction_numeric(models[[reg]], data_model, reg, var1, var2)
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    } else if(is.factor(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      lv <- levels(data_model[[var1]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var1, var2, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, ": ", lv[lv1], " -> ", lv[lv2], " × ", var2, sep=""), Effect = val))
          }
        }
      }
    } else if(is.factor(data_model[[var2]]) && is.numeric(data_model[[var1]])) {
      lv <- levels(data_model[[var2]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var2, var1, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var2, ": ", lv[lv1], " -> ", lv[lv2], " × ", var1, sep=""), Effect = val))
          }
        }
      }
    } else {
      val <- NA
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    }
  }
}
effect_results_5 <- effect_results

# Model 6

predictors <- c(
  "Population.size",
  "Generation",
  "Method"
)

interaction_pairs <- list(
  c("Population.size", "Generation"),
  c("Population.size", "Method"),
  c("Generation", "Method")
)

regions <- c(
  "Sensitivity.X",
  "Precision.X",
  "Sensitivity.Y2",
  "Precision.Y2"
)

models <- list(
  "Sensitivity.X"  = model6_SX,
  "Precision.X"    = model6_PX,
  "Sensitivity.Y2" = model6_SY,
  "Precision.Y2"   = model6_PY
)

data_model <- data_test6long

effect_results <- data.frame()
for (reg in regions) {
  cat("Processing:", reg, "\n")
  # MAIN EFFECTS
  for (var in predictors) {
    if(is.factor(data_model[[var]])) {
      lv <- levels(data_model[[var]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_main(models[[reg]], data_model, reg, var, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var, ": ", lv[lv1], " -> ", lv[lv2], sep=""), Effect = val))
          }
        }
      }
    } else {
      val <- contrast_main(models[[reg]], data_model, reg, var, NA, NA)
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = var, Effect = val))
    }
  }
  # INTERACTIONS
  for (pair in interaction_pairs) {
    var1 <- pair[1]
    var2 <- pair[2]
    if (is.numeric(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      val <- contrast_interaction_numeric(models[[reg]], data_model, reg, var1, var2)
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    } else if(is.factor(data_model[[var1]]) && is.numeric(data_model[[var2]])) {
      lv <- levels(data_model[[var1]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var1, var2, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, ": ", lv[lv1], " -> ", lv[lv2], " × ", var2, sep=""), Effect = val))
          }
        }
      }
    } else if(is.factor(data_model[[var2]]) && is.numeric(data_model[[var1]])) {
      lv <- levels(data_model[[var2]])
      for (lv1 in seq_along(lv)) {
        for (lv2 in seq_along(lv)) {
          if (lv2 > lv1) {
            val <- contrast_interaction_1fac(models[[reg]], data_model, reg, var2, var1, lv1, lv2)
            effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var2, ": ", lv[lv1], " -> ", lv[lv2], " × ", var1, sep=""), Effect = val))
          }
        }
      }
    } else {
      val <- NA
      effect_results <- rbind(effect_results, data.frame(Region = reg, Term = paste(var1, var2, sep = " × "), Effect = val))
    }
  }
}
effect_results_6 <- effect_results




# Plot heat map
effect_results <- rbind(effect_results_3, effect_results_4[grep("Model", effect_results_4$Term),], effect_results_5[grep("Pairs", effect_results_5$Term),])
effect_results <- rbind(effect_results, effect_results_6[grep("Method", effect_results_6$Term),])
effect_results <- effect_results[which(effect_results$Region != "Sensitivity.A" & effect_results$Region != "Precision.A"),]

effect_results$Term <- factor(effect_results$Term, order=T, levels=rev(c("Population.size", "Generation", "WhatsHap: OFF -> ON", "MAC",
                                                          "Window.size", "Model: Hamming -> Inverse_MAF", "Pairs",
                                                          "Population.size × Generation","WhatsHap: OFF -> ON × Population.size", "Population.size × MAC", 
                                                          "Population.size × Window.size", "Model: Hamming -> Inverse_MAF × Population.size", "Population.size × Pairs", 
                                                          "WhatsHap: OFF -> ON × Generation", "Generation × MAC", 
                                                          "Generation × Window.size", "Model: Hamming -> Inverse_MAF × Generation", "Generation × Pairs",
                                                          "WhatsHap: OFF -> ON × MAC", "WhatsHap: OFF -> ON × Window.size", "MAC × Window.size",
                                                          "Model: Hamming -> Inverse_MAF × Window.size", "Pairs × Window.size",
                                                          "Method: Python script -> PhaseWY 1 female and 2 males", "Method: Python script -> PhaseWY 2 females and 1 male", 
                                                          "Method: PhaseWY 1 female and 2 males -> PhaseWY 2 females and 1 male", "Method: Python script -> PhaseWY 1 female and 2 males × Population.size",
                                                          "Method: Python script -> PhaseWY 2 females and 1 male × Population.size", "Method: PhaseWY 1 female and 2 males -> PhaseWY 2 females and 1 male × Population.size",
                                                          "Method: Python script -> PhaseWY 1 female and 2 males × Generation", "Method: Python script -> PhaseWY 2 females and 1 male × Generation",
                                                          "Method: PhaseWY 1 female and 2 males -> PhaseWY 2 females and 1 male × Generation")),


                                                        labels=rev(c("log10(Population size)", "log10(Generations)", "WhatsHap (OFF -> ON)", "MAC filter",
                                                                     "log10(Window size)", "Model (Hamming -> Inverse_MAF)", "Female-Male pairs",
                                                                     "log10(Population size) × log10(Generations)","log10(Population size) × WhatsHap (OFF -> ON)", "log10(Population size) × MAC filter",
                                                                     "log10(Population size) × log10(Window size)", "log10(Population size) × Model (Hamming -> Inverse_MAF)", "log10(Population size) × Female-Male pairs",
                                                                     "log10(Generations) × WhatsHap (OFF -> ON)", "log10(Generations) × MAC filter", 
                                                                     "log10(Generations) × log10(Window size)", "log10(Generations) × Model (Hamming -> Inverse_MAF)", "log10(Generations) × Female-Male pairs",
                                                                     "WhatsHap (OFF -> ON) × MAC filter", "WhatsHap (OFF -> ON) × log10(Window size)", "MAC filter × log10(Window size)",
                                                                     "log10(Window size) × Model (Hamming -> Inverse_MAF)", "log10(Window size) × Female-Male pairs",
                                                                     "Method (Python script -> PhaseWY 1 female and 2 males)", "Method (Python script -> PhaseWY 2 females and 1 male)", 
                                                                     "Method (PhaseWY 1 female and 2 males -> PhaseWY 2 females and 1 male)", "Population.size × Method (Python script -> PhaseWY 1 female and 2 males)",
                                                                     "Population.size × Method (Python script -> PhaseWY 2 females and 1 male)", "Population.size × Method (PhaseWY 1 female and 2 males -> PhaseWY 2 females and 1 male)",
                                                                     "Generation × Method (Python script -> PhaseWY 1 female and 2 males)", "Generation × Method (Python script -> PhaseWY 2 females and 1 male)",
                                                                     "Generation × Method (PhaseWY 1 female and 2 males -> PhaseWY 2 females and 1 male)")))
effect_results$Region <- factor(effect_results$Region, levels = c("Sensitivity.Y2", "Precision.Y2", "Sensitivity.X", "Precision.X", "Sensitivity.A", "Precision.A"),
                                                       labels = c("Y Sensitivity", "Y Precision", "X Sensitivity", "X Precision", "A Sensitivity", "A Precision"))

#effect_results <- effect_results[-which((grepl("MAC", effect_results$Term) & grepl("×", effect_results$Term)) | (grepl("Model", effect_results$Term) & grepl("×", effect_results$Term)) | (grepl("pairs", effect_results$Term) & grepl("×", effect_results$Term))),]
#effect_results <- effect_results[-which(grepl("MAC", effect_results$Term) | grepl("Model", effect_results$Term)),]



effect_results_model <- effect_results[grep("Method", effect_results$Term),]
effect_results_param <- effect_results[-grep("Method", effect_results$Term),]
  
mx <- max(abs(effect_results_param$Effect), na.rm = TRUE)

plot_APC_param <- ggplot() +
  geom_tile(data=effect_results_param, aes(x = Region, y = Term, fill = Effect), color = "white", linewidth = 0.3) +
  #scale_fill_viridis_c(option = "viridis", limits = c(-mx, mx), name = "Δ Accuracy") +
  scale_fill_gradient2(low = "#1F77B4", mid = "#F8F8F8", high = "#FFB000", midpoint = 0, limits = c(-mx, mx), name = "Δ Accuracy") +
  scale_x_discrete(position = "top") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=12), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        title = element_text(size=20),
        axis.line = element_line(colour = "black"),
        axis.text.x.top  = element_text(size=15, color="black", face="bold"),
        axis.ticks.x.top = element_line(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.text.y = element_text(size=12, color="black"))

plot_APC_param

png("Figures/APC_parameters.png", width=3750, height=3000, res=300)
plot_APC_param
dev.off()

mx2 <- max(abs(effect_results_model$Effect), na.rm = TRUE)

plot_APC_model <- ggplot() +
  geom_tile(data=effect_results_model, aes(x = Region, y = Term, fill = Effect), color = "white", linewidth = 0.3) +
  #scale_fill_viridis_c(option = "viridis", limits = c(-mx, mx), name = "Δ Accuracy") +
  scale_fill_gradient2(low = "#1F77B4", mid = "#F8F8F8", high = "#FFB000", midpoint = 0, limits = c(-mx2, mx2), name = "Δ Accuracy") +
  scale_x_discrete(position = "top") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw()+
  theme(legend.key.size = unit(1, 'line'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15, face="bold"), #change legend title font size
        legend.text = element_text(size=12), #change legend text font size
        legend.box.spacing = unit(1, "cm"),
        strip.background = element_rect(color="black", fill="white", linewidth=1),
        strip.text = element_text(size=20, color="black", face="bold"),
        panel.spacing = unit(0, "lines"),
        title = element_text(size=20),
        axis.line = element_line(colour = "black"),
        axis.text.x.top  = element_text(size=12, color="black", face="bold"),
        axis.ticks.x.top = element_line(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.text.y = element_text(size=12, color="black"))

plot_APC_model
png("Figures/Supplement/APC_model.png", width=4250, height=3000, res=300)
plot_APC_model
dev.off()

#effect_results_6$formatted <- formatC(effect_results_6$Effect, digits = 3, format = "g")


