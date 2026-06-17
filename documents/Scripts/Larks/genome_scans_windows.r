#!/usr/bin/Rscript

## Export variables and load libraries
rm(list=ls())
options(scipen=999)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

WINDOW <- as.numeric(WINDOW)
STEP <- as.numeric(STEP)
MIN_SCAFFOLD_SIZE <- as.numeric(MIN_SCAFFOLD_SIZE)

### Import data
regions <- read.delim(MASK, sep="\t", head=F)
colnames(regions) <- c("Scaffold", "Start", "End")

ref_reg <- read.delim(paste(REF, ".fai", sep=""), sep="\t", head=F)
colnames(ref_reg) <- c("Scaffold", "Size")

# Set up a data frame
data <- data.frame(matrix(NA, 1, 6))
colnames(data) <- c("scaffold", "start", "end", "N_callable_sites", "N_tot_sites", "Abs_pos")
abs_pos_scaff <- 0

### Loop over scaffold
for(scaff in ref_reg$Scaffold) {

  ### Calculate total size of genomic region that were callable
  scaff_reg <- regions[which(regions$Scaffold == scaff),]
  if(nrow(scaff_reg) > 0) {
    size <- 0
    for(i in 1:nrow(scaff_reg)) {
      size <- size + (scaff_reg$End[i] - scaff_reg$Start[i])
    }

    if(ref_reg$Size[which(ref_reg$Scaffold == scaff)]  >= MIN_SCAFFOLD_SIZE) {

      ### Translate genomic coordinates into a simplified and continous sequence
      trans_cord <- rep(NA, size)
      start_cord <- 1
      for(i in 1:nrow(scaff_reg)) {
        gen_cord <- seq(scaff_reg$Start[i]+1, scaff_reg$End[i])
        trans_cord[start_cord:(length(gen_cord)+start_cord-1)] <- gen_cord
        start_cord <- start_cord + length(gen_cord)
      }

      ### Calculate sliding windows based on translated coordinates
      if(size > WINDOW) {
        startpos <- seq(1, size - WINDOW, STEP)
        endpos <- startpos + WINDOW - 1
        if(endpos[length(endpos)] < size) {
          startpos <- c(startpos, size - WINDOW + 1)
          endpos <- c(endpos, size)
        }

        ### Translate coordinates for sliding window into genomic coordinates
        for(i in 1:length(startpos)) {
          startpos[i] <- trans_cord[startpos[i]]
          endpos[i] <- trans_cord[endpos[i]]
        }
        N_callable_sites_end <- WINDOW
      } else {
        ### If Genomic region is smaller than window, make only one window
        startpos <- trans_cord[1]
        endpos <- trans_cord[size]
        N_callable_sites_end <- size
      }

      ### Set up data frame for windows
      scaff_wind <- data.frame(matrix(NA, length(startpos), 6))
      colnames(scaff_wind) <- c("scaffold", "start", "end", "N_callable_sites", "N_tot_sites", "Abs_pos")
      scaff_wind$scaffold <- scaff
      scaff_wind$start <- startpos
      scaff_wind$end <- endpos
      scaff_wind$N_callable_sites <- WINDOW
      scaff_wind$N_callable_sites[nrow(scaff_wind)] <- N_callable_sites_end
      scaff_wind$N_tot_sites <- endpos - startpos + 1
      scaff_wind$Abs_pos <- abs_pos_scaff + startpos + round((endpos - startpos + 1)/2)
      data <- rbind(data, scaff_wind)
      abs_pos_scaff <- abs_pos_scaff + ref_reg$Size[which(ref_reg$Scaffold == scaff)]
    }
  }
}
data <- data[-1,]

### Write file
write.table(data, file=paste(OUTDIR1, "/windows/", DATA, "_windows_", WINDOW, "_steps_", STEP, "_exon_dist_", EXON_DIST, ".txt", sep=""), quote=FALSE, sep="\t", row.names = F, col.names = T)

quit()
