#!/bin/bash -l

#SBATCH -A naiss2025-5-344
#SBATCH -p shared
#SBATCH -n 128
#SBATCH -t 100:00:00
#SBATCH -J coinheritance
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

# Set parameters
MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=larks_2026;
WORKDIR=$MAINDIR/data/$PROJECT;
VCFS1=$WORKDIR/Rasolark_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/vcfs;
VCFS2=$WORKDIR/Skylark_Europe_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/vcfs;
VCFS3=$WORKDIR/A7_filter_variants_mito/filter0;
OUTDIR=$MAINDIR/working/$PROJECT;
REF=$MAINDIR/data/reference/Alauda_arvensis\
/Alauda_arvensis_M_hifiasm-purged-default_hap0.purged_no_mito.yahs_r2.scf.FINAL_mito.fasta;
METADATA=$WORKDIR/metadata;
FUNCTIONS=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/bin_general;

# Define functions
VCF2PHYLIP=$FUNCTIONS/vcf2phylip-accessed-2024-02-22/vcf2phylip.py;

### Load modules
conda activate filter_variants;

# Load vt
export PATH=$PATH:/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/softwares/vt;


### Create folders
mkdir $OUTDIR/coinheritance;
OUTDIR=$OUTDIR/coinheritance;

# Females
FEMALES1=$(cat $METADATA/Rasolark_sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}' | tr '\n' ',' | head -c-1);
FEMALES2=$(cat $METADATA/Skylark_sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}' | tr '\n' ',' | head -c-1);

# Filter data

# Z
bcftools view $VCFS1/homogametic_filtered.vcf.gz -r scaffold_1:15917500-213977500 -s $FEMALES1 | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/Rasolark_Z_filtered_snps.vcf.gz;
tabix $OUTDIR/Rasolark_Z_filtered_snps.vcf.gz;


bcftools view $VCFS2/homogametic_filtered.vcf.gz -r scaffold_1:15917500-213977500 -s $FEMALES2 | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/Skylark_Z_filtered_snps.vcf.gz;
tabix $OUTDIR/Skylark_Z_filtered_snps.vcf.gz;


bcftools merge $VCFS1/homogametic_filtered.vcf.gz $VCFS2/homogametic_filtered.vcf.gz -r scaffold_1:15917500-213977500 | \
bcftools view -s $FEMALES1,$FEMALES2 | \
vcffixup - | \
vcfclassify - | \
vcffilter \
-s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/Z_filtered_snps.vcf.gz;
tabix $OUTDIR/Z_filtered_snps.vcf.gz;


# W
bcftools view $VCFS1/heterogametic_filtered.vcf.gz -r scaffold_1:15917500-213977500 -s $FEMALES1 | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/Rasolark_W_filtered_snps.vcf.gz;
tabix $OUTDIR/Rasolark_W_filtered_snps.vcf.gz;


bcftools view $VCFS2/heterogametic_filtered.vcf.gz -r scaffold_1:15917500-213977500 -s $FEMALES2 | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/Skylark_W_filtered_snps.vcf.gz;
tabix $OUTDIR/Skylark_W_filtered_snps.vcf.gz;


bcftools merge $VCFS1/heterogametic_filtered.vcf.gz $VCFS2/heterogametic_filtered.vcf.gz -r scaffold_1:15917500-213977500 | \
bcftools view -s $FEMALES1,$FEMALES2 | \
vcffixup - | \
vcfclassify - | \
vcffilter \
-s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/W_filtered_snps.vcf.gz;
tabix $OUTDIR/W_filtered_snps.vcf.gz;


# Mitochondrion
bcftools view $VCFS3/Introgression_2025_0.vcf.gz -r ptg000978l_1:1-18000 -s $FEMALES1,$FEMALES2 | \
vcftools --gzvcf - \
--minQ 30 \
--minDP 3 \
--recode --recode-INFO-all --stdout | \
vcffilter \
-f "AC > 0" | \
vcfallelicprimitives --keep-info --keep-geno | \
vt decompose_blocksub - -o + | \
vt normalize + -m -r $REF -o + | \
bcftools norm --rm-dup all -Ov | \
awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("\\|","/",$i)} print}' | \
vcffixup - | \
vcfclassify - | \
vcffilter \
-f "AC > 0" | \
bgzip -c > $OUTDIR/Mitochondrion_filtered.vcf.gz;
tabix $OUTDIR/Mitochondrion_filtered.vcf.gz;

bcftools view $OUTDIR/Mitochondrion_filtered.vcf.gz -s $FEMALES1 | \
vcffixup - | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/Rasolark_M_filtered_snps.vcf.gz;
tabix $OUTDIR/Rasolark_M_filtered_snps.vcf.gz;

bcftools view $OUTDIR/Mitochondrion_filtered.vcf.gz -s $FEMALES2 | \
vcffixup - | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/Skylark_M_filtered_snps.vcf.gz;
tabix $OUTDIR/Skylark_M_filtered_snps.vcf.gz;

bcftools view $OUTDIR/Mitochondrion_filtered.vcf.gz | \
vcffixup - | \
vcfclassify - | \
vcffilter \
-s -f "!( INS | DEL | MNP )" \
-f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/M_filtered_snps.vcf.gz;
tabix $OUTDIR/M_filtered_snps.vcf.gz;

# Make tree
cd $OUTDIR;
python3 $VCF2PHYLIP -i $OUTDIR/Rasolark_Z_filtered_snps.vcf.gz;
python3 $VCF2PHYLIP -i $OUTDIR/Skylark_Z_filtered_snps.vcf.gz;
python3 $VCF2PHYLIP -i $OUTDIR/Z_filtered_snps.vcf.gz;
python3 $VCF2PHYLIP -i $OUTDIR/Rasolark_W_filtered_snps.vcf.gz;
python3 $VCF2PHYLIP -i $OUTDIR/Skylark_W_filtered_snps.vcf.gz;
python3 $VCF2PHYLIP -i $OUTDIR/M_filtered_snps.vcf.gz;
python3 $VCF2PHYLIP -i $OUTDIR/Rasolark_M_filtered_snps.vcf.gz;
python3 $VCF2PHYLIP -i $OUTDIR/Skylark_M_filtered_snps.vcf.gz;
conda deactivate;
conda activate tree_reconstruction;

# Run IQtree
iqtree2 -s $OUTDIR/Rasolark_Z_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/Rasolark_Z_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/Rasolark_Z_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/Skylark_Z_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/Skylark_Z_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/Skylark_Z_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/Z_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/Z_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/Z_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/Rasolark_W_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/Rasolark_W_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/Rasolark_W_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/Skylark_W_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/Skylark_W_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/Skylark_W_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/W_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/W_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/W_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/Rasolark_M_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/Rasolark_M_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/Rasolark_M_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/Skylark_M_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/Skylark_M_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/Skylark_M_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
iqtree2 -s $OUTDIR/M_filtered_snps.min4.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
if [ $(ls $OUTDIR/M_filtered_snps.min4.phy.varsites.phy | wc -l) -eq 1 ] ; then
  iqtree2 -s $OUTDIR/M_filtered_snps.min4.phy.varsites.phy -m TEST+ASC -nt AUTO -B 1000 -wbt;
fi;
cd -;
