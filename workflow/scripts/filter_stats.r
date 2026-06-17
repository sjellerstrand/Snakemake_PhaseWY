#!/usr/bin/Rscript

### Plot statistics
## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library (ggplot2)
library(dplyr)
library(gridExtra)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
assign(args[[i]][1], args[[i]][2])
}

## Variant based statistics

# Variant missingness
var_miss <- read_delim(paste(VCF_OUT, ".lmiss", sep=""), delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a3 <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("var miss") + theme_light()
b3 <- summary(var_miss$fmiss)

# Minor allele frequency
var_freq <- read_delim(paste(VCF_OUT, ".frq", sep=""), delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
a4 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("var maf") + theme_light()
b4 <- summary(var_freq$maf)

## Individual based statistics

# Heterozygosity and inbreeding coefficient per individual
if(file.exists(paste(VCF_OUT, ".het", sep=""))) {
     ind_het <- read_delim(paste(VCF_OUT, ".het", sep=""), delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
     a5 <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("ind fstat") + theme_light()
     b5 <-summary(ind_het$f)
} else {
    a5 <- ggplot() + theme_void()
    b5 <- rep(NA, 6)
}

# Proportion of missing data per individual
ind_miss <- read_delim(paste(VCF_OUT, ".imiss", sep=""), delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a7 <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("ind miss") + theme_light()
b7<- summary(ind_miss$fmiss)

## Write plots as .png
sumstat <-rbind(b3, b4, b5, b7)
sumstat <- cbind(c("var miss", "var maf", "ind fstat", "ind miss"), sumstat)
write.table(sumstat, file=paste(VCF_OUT,"_stats.txt", sep=""), quote=F, sep='\t', row.names = F, col.names = T)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(sumstat, rows=NULL, theme=tt)
png(paste(VCF_OUT,"_stats.png", sep=""), width=1300, height=1200, res=300)
print(grid.arrange(a3, a7, a4, a5, tbl, nrow=3, ncol=2))
dev.off()
