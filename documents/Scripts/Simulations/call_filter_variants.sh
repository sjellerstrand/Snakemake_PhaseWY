q#!/bin/bash -l

#SBATCH -A naiss2025-5-344
#SBATCH -p main
#SBATCH -t 05:00:00
#SBATCH -J call_filter_variants
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=simulations_2025;
PROJECT2=simulations_2025;
WORKDIR=$MAINDIR/data/$PROJECT;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/$PROJECT;
METADATA=$WORKDIR/metadata;
REF=$WORKDIR/simulate_reads/alignments/reference.fasta;
FUNCTIONS=$MAINDIR/scripts/$PROJECT;

# Load modules
ml PDCOLD/23.12 R/4.4.0;

# Load vt
export PATH=$PATH:/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/softwares/vt;

### Define functions
filter_stats=$FUNCTIONS/quality_control/filter_stats.sh;
pca=$FUNCTIONS/quality_control/pca.sh;

# Activate conda environment
conda activate call_variants_bcftools;

### Create folders
mkdir $OUTDIR/call_filter_variants;
OUTDIR=$OUTDIR/call_filter_variants;

### Setup file info
INDS=$(cat $METADATA/sample_info.txt | cut -f1);
for IND in ${INDS[@]}; do
  echo $WORKDIR/simulate_reads/alignments/$IND\_merged.bam;
done > $OUTDIR/BAM_list.txt;
BAMS=$OUTDIR/BAM_list.txt;

## Set maximum coverage allowed
echo "A maximum of 1000 reads are processed per site and sample";

### Call variants with bcftools
bcftools mpileup -f $REF -b $BAMS -Q 30 -q 20 -d 1000 -B \
-a FORMAT/DP,INFO/AD,INFO/ADF,INFO/ADR -Ou | \
bcftools call -m -M -v --threads 128 -Ou | \
bcftools filter -e 'INFO/DP < 3' -Oz \
> $OUTDIR/$PROJECT\_raw.vcf.gz;
tabix $OUTDIR/$PROJECT\_raw.vcf.gz;

conda deactivate;

# Filter variants
conda activate filter_variants;

## Excess heterozygosity
vcftools --gzvcf $OUTDIR/$PROJECT\_raw.vcf.gz --hardy --stdout |tail -n+2 | cut -f1,2,3 | \
awk -F'\t|/' 'BEGIN {print "#Header"} {if($3 == 0 && $5 == 0) print $1"\t"$2-1"\t"$2}' \
> $OUTDIR/excess_het_with_header.bed;

## Filter
bcftools view $OUTDIR/$PROJECT\_raw.vcf.gz | \
vcftools --vcf - \
--exclude-bed $OUTDIR/excess_het_with_header.bed \
--max-missing 0.95 \
--min-meanDP 3 \
--max-meanDP 30 \
--minQ 30 \
--minDP 3 \
--recode --recode-INFO-all --stdout | \
vcfallelicprimitives --keep-info --keep-geno | \
vt decompose_blocksub - -o + | \
vt normalize + -m -r $REF -o + | \
bcftools norm --rm-dup all -Ov | \
awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("\\|","/",$i)} print}' | \
vcffixup - | \
vcfclassify - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $OUTDIR/$PROJECT\_filtered.vcf.gz;
tabix $OUTDIR/$PROJECT\_filtered.vcf.gz;
VCF_OUT=$OUTDIR/$PROJECT\_filtered;
source $filter_stats;
source $pca;

# Subset to smaller datasets
bcftools view $VCF_OUT.vcf.gz -s ^Ind_9_Female,Ind_10_Female,Ind_19_Male,Ind_20_Male | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $VCF_OUT\8.vcf.gz;
tabix $VCF_OUT\8.vcf.gz;


bcftools view $VCF_OUT\8.vcf.gz -s ^Ind_7_Female,Ind_8_Female,Ind_17_Male,Ind_18_Male | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $VCF_OUT\6.vcf.gz;
tabix $VCF_OUT\6.vcf.gz;

bcftools view $VCF_OUT\6.vcf.gz -s ^Ind_5_Female,Ind_6_Female,Ind_15_Male,Ind_16_Male | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $VCF_OUT\4.vcf.gz;
tabix $VCF_OUT\4.vcf.gz;

bcftools view $VCF_OUT\4.vcf.gz -s ^Ind_3_Female,Ind_4_Female,Ind_13_Male,Ind_14_Male | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $VCF_OUT\2.vcf.gz;
tabix $VCF_OUT\2.vcf.gz;

bcftools view $VCF_OUT\2.vcf.gz -s ^Ind_12_Male | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $VCF_OUT\1_het_2_hom.vcf.gz;
tabix $VCF_OUT\1_het_2_hom.vcf.gz;

bcftools view $VCF_OUT\2.vcf.gz -s ^Ind_2_Female | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $VCF_OUT\2_het_1_hom.vcf.gz;
tabix $VCF_OUT\2_het_1_hom.vcf.gz;

bcftools view $VCF_OUT\2.vcf.gz -s Ind_11_Male,Ind_1_Female | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $VCF_OUT\1.vcf.gz;
tabix $VCF_OUT\1.vcf.gz;
