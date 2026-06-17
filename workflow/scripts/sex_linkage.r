#!/usr/bin/Rscript

# Version 2025-10-21
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

## Export variables and load libraries
rm(list=ls())
library(vcfR)
library(data.table)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

options(scipen=999)
set.seed(123456)


SEX_DEPTH_THRESH <- as.numeric(SEX_DEPTH_THRESH)
CONTIG_LENGTH <- as.numeric(CONTIG_LENGTH)
if(!exists("MODEL")) {
  MODEL <- "Hamming"
}

### Import data
vcf1 <- read.vcfR(paste(INDIR1, "/", CONTIG, "_phased_all_variants_mac.vcf.gz", sep=""), verbose = F)
if(nrow(extract.gt(vcf1)) == 1) {
  vcf1 <- cbind(t(as.matrix(getFIX(vcf1))), as.matrix(extract.gt(vcf1, element="GT")))
} else {
  vcf1 <- cbind(getFIX(vcf1), extract.gt(vcf1, element="GT"))
}
sample_info <- read.table(INDS, sep='\t', head=T)
sex_depth <- read.table(paste(INDIR2, "/", CONTIG, "/", CONTIG, "_sex_depth_windows.bed", sep=""), sep='\t', head=F)

# Find heterogametic individuals
heterogametes <- sample_info[which(sample_info$sex == "HETGAM"),1]
homogametes <- sample_info[which(sample_info$sex == "HOMGAM"),1]

# Seperate haplotypes
haplotypes <- as.data.frame(t(cbind(matrix(vcf1[,1]), matrix(vcf1[,2]), matrix(NA, nrow(vcf1), (ncol(vcf1)-7)*2))))
rownames(haplotypes)[c(1,2)] <- c("Contig", "Pos")
for(i in 1:(ncol(vcf1)-7)) {
  rownames(haplotypes)[c(2+(i*2-1),2+(i*2))] <- c(paste(colnames(vcf1)[7+i], "_left", sep="") ,paste(colnames(vcf1)[7+i], "_right", sep=""))
  haplotypes[2+(i*2-1),] <- unlist(strsplit(vcf1[,7+i], "|", fixed=T))[c(T,F)]
  haplotypes[2+(i*2),] <- unlist(strsplit(vcf1[,7+i], "|", fixed=T))[c(F,T)]
}
rm(vcf1)
positions <- as.numeric(haplotypes[2,])
haplotype_names <- rownames(haplotypes)[3:nrow(haplotypes)]

### Set up data frame
startpos <- sex_depth[,2]+1
endpos <- sex_depth[,3]
WIND_OVERLAP <- ((endpos[1]-startpos[1]+1)/(startpos[2]-startpos[1]))

sexlink <- matrix(NA, nrow(sex_depth), (21+length(heterogametes)+1))
colnames(sexlink) <- c("Start pos [bp]", "End pos [bp]", "Sex-linked status", "Number of variants in window", "Ratio of total heterogametic individuals in smallest cluster [#inds/#inds]", "Ratio of haplotypes in smallest cluster [#haps/#haps]", "Smallest cluster [#inds]", "Largest cluster [#inds]", "TotSS", "Smallest cluster WithinSS", "Largest cluster WithinSS", "BetweenSS", "Window length [bp]", "Variants/bp", "Window mid position [bp]", "Phase information available", "Border change", "Phase switch", "Unknown", "Depth difference", "Heterozygosity difference", heterogametes, "sexdepth")
sexlink[,1] <- sex_depth[,2]+1
sexlink[,2] <- sex_depth[,3]
sexlink[,20] <- sex_depth[,4]
rm(sex_depth)
phswitch <- matrix(nrow=0, ncol=3)

### Evaluate windows for sex linkage in sliding windows
i <- startpos[1]
stop_search <- NA
phase_switch2 <- NA
border2 <- NA
set.seed(123456) # As k-means clustering algorithm starts with k randomly selected centroids, set a seed for R’s random number generator to make the results reproducible
while(is.na(stop_search)) {
  
  # Reset some parameters
  i0 <- i
  status <- NA
  phase <- matrix("Unknown", 1, length(heterogametes))
  colnames(phase) <- colnames(sexlink)[22:(ncol(sexlink)-1)]
  
  # Evaluate K-means
  snps <- as.data.frame(haplotypes[3:nrow(haplotypes), which(between(positions, startpos[i], endpos[i]))], row.names=haplotype_names)
  if(ncol(snps) == 0) {
    status <- "Autosomal: no variation in window"
  } else {
    snps_in <- as.data.frame(lapply(snps, as.numeric), row.names = rownames(snps))
    
    if(MODEL== "Inverse_MAF") { #Inverse maf weighting, giving rare alleles higher weight
      maf <- apply(snps_in, 2, function(x) {
        freq <- mean(x, na.rm = TRUE)
        pmin(freq, 1 - freq)
      })
      weights <- 1 / (maf +  1e-6) # small number to avoid division by zero
      snps_in <- sweep(snps_in, 2, sqrt(weights), `*`)
    }
    
    km <- kmeans(snps_in, centers=2, nstart=10)
    # Find smallest cluster
    if(km$size[1] == km$size[2]) {
      min <- 1
      max <- 2
    } else {
      min <- which(km$size == min(km$size))
      max <- which(km$size == max(km$size))
    }
    min_cluster_temp <- strsplit(names(km$cluster)[which(km$cluster == min)], "_", fixed=T)
    min_cluster <- rbind(sapply(lapply(min_cluster_temp, head, -1), paste, collapse="_"), sapply(min_cluster_temp, tail, 1))
    
    # Are any homogametes clustered in the smallest cluster?
    for(j in 1:length(homogametes)) {
      if(length(which(min_cluster[1,] == homogametes[j])) > 0) {
        status <- "Autosomal: homogamete in cluster"
        break
      }
    }
    
    # Count and evaluate heterogamete haplotypes in the smallest cluster
    heterogametes_count <- 0
    for(j in 1:length(heterogametes)) {
      if(length(which(min_cluster[1,] == heterogametes[j])) > 0) {
        
        # Do both haploypes of any heterogametes occur in cluster?
        if(length(which(min_cluster[1,] == heterogametes[j])) > 1 && is.na(status)) {
          status <- "Autosomal: both haplotypes of heterogamete in cluster"
        }
        heterogametes_count <- heterogametes_count + 1
      }
    }
    
    # Does the fraction of heterogametes in the smallest cluster meet the set threshold for it to be classed as sex-linked?
    if(heterogametes_count/length(heterogametes) == 1 && is.na(status)) {
      status <- "Sex-linked"
    } else if(is.na(status)) {
      status <- "Autosomal: too few heterogametes in cluster"
    }
  }
  # If the current window is sex-linked, do the following:
  if(status == "Sex-linked") {
    # Find phase
    for(j in 1:ncol(min_cluster)) {
      phase[1,which(colnames(phase) == min_cluster[1,j])] <- min_cluster[2,j]
    }
  }

  # Calculate sex differences in heterozygosity
  hetgam_hets <- rep(NA, length(heterogametes))
  for(j in 1:length(heterogametes)) {
    genotypes <- t(snps[which(sapply(strsplit(rownames(snps), "_", fixed = TRUE), function(x) paste(head(x, -1), collapse = "_")) == heterogametes[j]),])
    hetgam_hets[j] <- nrow(genotypes[apply(genotypes, 1, function(x) x[1] != x[2]), , drop = FALSE])
  }

  homgam_hets <- rep(NA, length(homogametes))
  for(j in 1:length(heterogametes)) {
    genotypes <- t(snps[which(sapply(strsplit(rownames(snps), "_", fixed = TRUE), function(x) paste(head(x, -1), collapse = "_")) == homogametes[j]),])
    homgam_hets[j] <- nrow(genotypes[apply(genotypes, 1, function(x) x[1] != x[2]), , drop = FALSE])
  }
  
  # Register data for window
  sexlink[i,1] <- startpos[i]
  sexlink[i,2] <- endpos[i]
  sexlink[i,3] <- status
  sexlink[i,4] <- ncol(snps)
  if(status == "Autosomal: no variation in window") {
    sexlink[i,5:12] <- NA
  } else {
    sexlink[i,5] <- heterogametes_count/length(heterogametes)
    sexlink[i,6] <- km$size[min]/(km$size[min] + km$size[max])
    sexlink[i,7] <- km$size[min]
    sexlink[i,8] <- km$size[max]
    sexlink[i,9] <- km$totss
    sexlink[i,10] <- km$withinss[min]
    sexlink[i,11] <- km$withinss[max]
    sexlink[i,12] <- km$betweenss
  }
  sexlink[i,13] <-  endpos[i] - startpos[i] + 1
  sexlink[i,14] <- as.numeric(sexlink[i,4]) / as.numeric(sexlink[i,13])
  sexlink[i,15] <- as.numeric(startpos[i] + (as.numeric(sexlink[i,13])-1)/2)
  if(length(which(is.na(phase))) == 0 && length(which(phase == "Unknown")) == 0) {
    sexlink[i,16] <- "Yes"
  } else {
    sexlink[i,16] <- "No"
  }
  
  if(i > 1 && ((length(which(sexlink[i-1, 22:(ncol(sexlink)-1)] != phase[1,])) > 0) || (sexlink[i-1,3] == "Sex-linked" && sexlink[i,3] != "Sex-linked"))) {
    
    if(i-WIND_OVERLAP >= 0) {
      startoverlap <- floor(i-WIND_OVERLAP)
    } else {
      startoverlap <- 1
    }
    if(i+WIND_OVERLAP-1 > nrow(sexlink)) {
      endoverlap <- nrow(sexlink)
    } else {
      endoverlap <- ceiling(i+WIND_OVERLAP-1)
    }
    
    # If border change
    
    # If transition from autosomal region into sex-linked region
    if(sexlink[i,16] == "Yes" && sexlink[i-1,16] != "Yes") {
      
      # Is autosomal region actually sex linked due to sex-depth difference?
      if(length(which(!is.na(sexlink[(startoverlap+1):(i-1),20]))) > 0 && mean(as.numeric(sexlink[(startoverlap+1):(i-1),20]), na.rm=T) >= SEX_DEPTH_THRESH) {
        border1 <- startoverlap
        border2 <- endoverlap
      } else {
        sexlink[i,(ncol(sexlink))] <- sexlink[i,3]
        sexlink[i,3] <- "Sex-linked"
        sexlink[i,16] <- "Yes"
      }
      
      # If transition from sex-linked region into autosomal region
    } else if(sexlink[i,16] != "Yes" && sexlink[i-1,16] == "Yes") {
      
      # Is autosomal region actually sex linked due to sex-depth difference?
      if(length(which(!is.na(sexlink[i:(endoverlap-1),20]))) > 0 && mean(as.numeric(sexlink[i:(endoverlap-1),20]), na.rm=T) >= SEX_DEPTH_THRESH) {
        border1 <- startoverlap
        border2 <- endoverlap
      } else {
        sexlink[i,(ncol(sexlink))] <- sexlink[i,3]
        sexlink[i,3] <- "Sex-linked"
        sexlink[i,16] <- "Yes"
      }
      
      # If phase switch
    } else if(sexlink[i,16] == "Yes" && sexlink[i-1,16] == "Yes" && length(which(sexlink[i-1, 22:(ncol(sexlink)-1)] == "Unknown")) != length(heterogametes) && length(which(phase[1,] == "Unknown")) != length(heterogametes)) {
      phase_switch1 <- startoverlap
      phase_switch2 <- endoverlap
      phase_switch_mat <- sexlink[i-1, 22:(ncol(sexlink)-1)]
    } else if(sexlink[i,16] == "Yes" && sexlink[i-1,16] == "Yes" && length(which(sexlink[i-1, 22:(ncol(sexlink)-1)] == "Unknown")) == length(heterogametes) && length(which(phase[1,] == "Unknown")) != length(heterogametes)) {
      phase_mat <- sexlink[(startoverlap+1):(i-1), 22:(ncol(sexlink)-1)] 
      phase_mat <- phase_mat[rowSums(phase_mat != phase[1,]) > 0 & rowSums(phase_mat == "Unknown") < length(heterogametes),, drop=F]
      if(nrow(phase_mat) > 0) {
        phase_switch1 <- startoverlap - floor(WIND_OVERLAP-nrow(phase_mat))
        if(phase_switch1 < 1) {
          phase_switch1 <- 1
        }
        phase_switch2 <- endoverlap
        phase_switch_mat <- phase_mat[nrow(phase_mat),]
      }
    }
  }
  if(!is.na(border2) && i <= border2) {
    sexlink[border1:i, 17] <- "Yes"
    if(i == border2) {
      border2 <- NA
    }
  }
  if(!is.na(phase_switch2) && i <= phase_switch2) {
    sexlink[phase_switch1:i, 18] <- "Yes"
    phswitch <- rbind(phswitch, cbind(rep(phase_switch1, length(which(phase_switch_mat != phase[1,]))),rep(phase_switch2, length(which(phase_switch_mat != phase[1,]))), unname(which(phase_switch_mat != phase[1,])) + 21))
    if(i == phase_switch2) {
      phase_switch2 <- NA
    }
  }
  sexlink[i,21] <- mean(hetgam_hets)/mean(homgam_hets)
  sexlink[i,22:(ncol(sexlink)-1)] <- phase[1,]
  
  # Continue loop?
  if(endpos[i] == endpos[length(endpos)]) {
    stop_search <- "Stop"
  } else {
    i <- i0+1
  }
}

# Set phase switches to unknown for respecive individual
if(nrow(phswitch) > 0) {
  for(i in 1:nrow(phswitch)) {
    startp <- phswitch[i,1]
    endp <- phswitch[i,2]
    sexlink[startp:endp, phswitch[i,3]] <- "Unknown"
  }
}

sexlink[is.na(sexlink[,17]), 17] <- "No"
sexlink[which(!is.na(sexlink[,18]) & sexlink[,16] == "No"), 18] <- "No"
sexlink[is.na(sexlink[,18]), 18] <- "No"
sexlink[which(sexlink[,17] == "Yes" | sexlink[,18] == "Yes"), 19] <- "Yes"
sexlink[is.na(sexlink[,19]), 19] <- "No"

# Summarize coherent regions, including unknown at borders
regions <- matrix(NA, 1, 4)
colnames(regions) <- c("chrom", "chromStart", "chromEnd", "Status")

# If there is any sex-linkage
if(length(which(sexlink[,3] == "Sex-linked")) > 0) {
  stop_search <- NA
  i <- 1
  
  while(is.na(stop_search)) {
    i0 <- i
    beginreg2 <- NA
    endreg2 <- NA
    status3 <- NA
    
    while(i != nrow(sexlink) && sexlink[i0,17] == sexlink[i,17]) {
      i <- i + 1
    }
    
    # Define sex-linkage status of region
    if(sexlink[i0, 17] == "Yes") {
      status3 <- "Unknown"
    } else if(sexlink[i0,3] == "Sex-linked") {
      status3 <- "Sex-linked"
    } else {
      status3 <- "Autosomal"
    }
    
    # Define start of region
    if(i0 == 1) {
      beginreg2 <- 0
    } else {
      beginreg2 <- as.numeric(regions[nrow(regions)-1,3]) + 1 - 1
    }
    
    # Define end of region
    if(i == nrow(sexlink) && sexlink[i0, 17] == sexlink[i, 17]) {
      regions[nrow(regions),] <- c(CONTIG, beginreg2, CONTIG_LENGTH, status3)
    } else if(i == nrow(sexlink) && sexlink[i0, 17] != sexlink[i, 17]) {
      endreg2 <- as.numeric(sexlink[i-1,2])
      if(status3 == "Unknown") {
        endreg2 <- as.numeric(sexlink[i-WIND_OVERLAP,2])
        regions[nrow(regions),] <- c(CONTIG, beginreg2, endreg2, status3)
        regions <- rbind(regions, c(CONTIG, endreg2, CONTIG_LENGTH, status3))
      } else {
        endreg2 <- as.numeric(sexlink[i-1,2])
        regions[nrow(regions),] <- c(CONTIG, beginreg2, endreg2, status3)
        regions <- rbind(regions, c(CONTIG, endreg2, CONTIG_LENGTH, status3))
      }
    } else if(status3 == "Unknown") {
      endreg2 <- as.numeric(sexlink[i-WIND_OVERLAP,2])
      regions[nrow(regions),] <- c(CONTIG, beginreg2, endreg2, status3)
    } else {
      endreg2 <- as.numeric(sexlink[i-1,2])
      regions[nrow(regions),] <- c(CONTIG, beginreg2, endreg2, status3)
    }
    
    # Continue loop?
    if(i == nrow(sexlink)) {
      stop_search <- "Stop"
    } else {
      regions <- rbind(regions, matrix(NA, 1, 4))
    }
  }
  # If there is no evidence of sex-linkage, write whole contig as autosomal
} else {
  regions[1,] <- c(CONTIG, 0, CONTIG_LENGTH, "Autosomal")
}
if(is.matrix(regions)) {
  regions <- as.data.frame(regions)
} else {
  regions <- t(as.data.frame(regions))
}

colnames(regions)[1] <- paste("#", colnames(regions)[1], sep="")

# Write phase and sex-linkage information from search as bed files
sexlink[which(!is.na(sexlink[,ncol(sexlink)])),3] <- sexlink[which(!is.na(sexlink[,ncol(sexlink)])),ncol(sexlink)]
sexlink <- sexlink[,-ncol(sexlink)]
sexlink <- cbind("#chrom"=CONTIG, sexlink)
sexlink[,2] <- as.numeric(sexlink[,2])-1
write.table(sexlink, file=paste(OUTDIR, "/", CONTIG, "_phase_windows.bed", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)
write.table(regions, file=paste(OUTDIR, "/", CONTIG, "_phase_info.bed", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)
if(length(which(sexlink[,4] == "Sex-linked")) > 0) {
  for(IND in colnames(sexlink)[22:length(colnames(sexlink))]) {
    het_left <- sexlink[which(sexlink[,which(colnames(sexlink)==IND)] == "left"), c("#chrom", "Start pos [bp]", "End pos [bp]"), drop = FALSE]
    het_right <- sexlink[which(sexlink[, which(colnames(sexlink) == IND)] == "right"), c("#chrom", "Start pos [bp]", "End pos [bp]"), drop = FALSE]
    write.table(het_left, file=paste(OUTDIR, "/", CONTIG, "_", IND, "_het_left.bed", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)
    write.table(het_right, file=paste(OUTDIR, "/", CONTIG, "_", IND, "_het_right.bed", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)
    
  }
}
quit()