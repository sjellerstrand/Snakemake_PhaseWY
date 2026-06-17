#!/bin/bash -l

#SBATCH -A naiss2025-5-344
#SBATCH -p shared
#SBATCH -n 128
#SBATCH -t 100:00:00
#SBATCH -J divergence
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

# Set parameters
MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=larks_2026;
WORKDIR=$MAINDIR/data/$PROJECT;
VCFS1=$WORKDIR/Rasolark_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/vcfs;
VCFS2=$WORKDIR/Skylark_Europe_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/vcfs;
OUTDIR0=$MAINDIR/working/$PROJECT;
CALLABLE=$WORKDIR/Skylark_Europe_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/beds/heterogametic.bed;
REF=$MAINDIR/data/reference/Alauda_arvensis\
/Alauda_arvensis_M_hifiasm-purged-default_hap0.purged_no_mito.yahs_r2.scf.FINAL_mito.fasta;
CDS=$WORKDIR/B3_annotation_lift_over/Aarv_cds.bed;
GTF=$WORKDIR/B3_annotation_lift_over/Aarv_Tgut_liftover_polished_rmdup.gtf;
METADATA=$WORKDIR/metadata;
FUNCTIONS=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/bin_general;

### Load modules
conda activate filter_variants;

### Create folders
mkdir $OUTDIR0/divergence_CDS;
OUTDIR0=$OUTDIR0/divergence_CDS;

# Females
FEMALES2=$(cat $METADATA/Skylark_sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}' | tr '\n' ',' | head -c-1);

# Find ancestral strata
> $OUTDIR0/anc_strata_match.bed;

# Match genes to ancestral strata
for i in $(seq 2 1 $(cat $METADATA/Ancestral_strata_genes.tsv | tail -n+2 | wc -l)); do
  GENE=$(cat $METADATA/Ancestral_strata_genes.tsv | head -n$i | tail -n1 | cut -f1);
  STRATA=$(cat $METADATA/Ancestral_strata_genes.tsv | head -n$i | tail -n1 | cut -f2);
  cat $GTF | grep $GENE | awk -F'\t' -v GENE=$GENE -v STRATA=$STRATA '{if($9 ~ "^gene_id \"" GENE "\";" && $3 == "CDS") print $1"\t"$4-1"\t"$5"\t"GENE"\t"STRATA}'
done >> $OUTDIR0/anc_strata_match.bed;

bedtools intersect -a $CALLABLE -b $OUTDIR0/anc_strata_match.bed -wb | cut -f1,2,3,7,8 > $OUTDIR0/anc_strata_match_callable.bed;
cat $OUTDIR0/anc_strata_match_callable.bed | awk -F'\t' '{if($5 == "S3") print}' | bedtools sort > $OUTDIR0/S3.bed;
cat $OUTDIR0/anc_strata_match_callable.bed | awk -F'\t' '{if($5 == "S2") print}' | bedtools sort > $OUTDIR0/S2.bed;
cat $OUTDIR0/anc_strata_match_callable.bed | awk -F'\t' '{if($5 == "S1") print}' | bedtools sort > $OUTDIR0/S1.bed;
cat $OUTDIR0/anc_strata_match_callable.bed | awk -F'\t' '{if($5 == "S0") print}' | bedtools sort > $OUTDIR0/S0.bed;


# Set up mask files
cat $REF | awk 'BEGIN {notprint="F"} {if($0 ~ /^>/ && $0 != ">scaffold_1") {notprint="T"}; if(notprint=="F") print}' \
> $OUTDIR0/scaffold_1.fa;
cat $REF.fai | head -n1 > $OUTDIR0/scaffold_1.fa.fai;
REF=$OUTDIR0/scaffold_1.fa;
cat $CDS | awk -F '\t' '{if($1=="scaffold_1") print}' | bedtools intersect -a - -b $CALLABLE | bedtools sort > $OUTDIR0/cds.bed;
CDS=$OUTDIR0/cds.bed;
bedtools complement -i $CDS -g $REF.fai > $OUTDIR0/missing.bed;
MISSING=$OUTDIR0/missing.bed;
bedtools complement -i $OUTDIR0/S3.bed -g $REF.fai > $OUTDIR0/S3_mask.bed;
bedtools complement -i $OUTDIR0/S2.bed -g $REF.fai > $OUTDIR0/S2_mask.bed;
bedtools complement -i $OUTDIR0/S1.bed -g $REF.fai > $OUTDIR0/S1_mask.bed;
bedtools complement -i $OUTDIR0/S0.bed -g $REF.fai > $OUTDIR0/S0_mask.bed;


# Filter data

# Skylark
bcftools view $VCFS2/homogametic_filtered.vcf.gz -R $CDS | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" | \
vcftools --vcf - \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR0/Skylark_snps_Z.vcf.gz;
tabix $OUTDIR0/Skylark_snps_Z.vcf.gz;

bcftools view $VCFS2/heterogametic_filtered.vcf.gz -h | \
tail -n1 | cut -f 10- | tr '\t' '\n' | awk '{print $1" "$1"_W"}' > $OUTDIR0/W_head.txt;

bcftools reheader $VCFS2/heterogametic_filtered.vcf.gz -s $OUTDIR0/W_head.txt | \
bcftools view -T $CDS | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" | \
vcftools --vcf - \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR0/Skylark_snps_W.vcf.gz;
tabix $OUTDIR0/Skylark_snps_W.vcf.gz;

bcftools merge $OUTDIR0/Skylark_snps_Z.vcf.gz $OUTDIR0/Skylark_snps_W.vcf.gz -R $CDS | \
vcffixup - | \
vcffilter -f "( AC > 0 & AF < 1 )" | \
vcftools --vcf - \
--max-missing 0.95 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR0/Skylark_filtered_snps.vcf.gz;
tabix $OUTDIR0/Skylark_filtered_snps.vcf.gz;


conda deactivate;
conda activate gene_extraction;


### Stratum 5
mkdir $OUTDIR0/5;
OUTDIR=$OUTDIR0/5;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -r scaffold_1:15900000-52100000 | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_5.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_5.vcf.gz;
echo scaffold_1:15900000-52100000 > $OUTDIR/5_reg.txt;
> $OUTDIR/5.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/5.fasta;
  samtools faidx $REF -r $OUTDIR/5_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_5.vcf.gz -m $MISSING --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?" >> $OUTDIR/5.fasta;
done;


### Stratum S3
mkdir $OUTDIR0/S3;
OUTDIR=$OUTDIR0/S3;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -R $OUTDIR0/S3.bed | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_S3.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_S3.vcf.gz;
echo scaffold_1:$(cat $OUTDIR0/S3.bed | cut -f2 | sort | head -n1)-$(cat $OUTDIR0/S3.bed | cut -f3 | sort | tail -n1) > $OUTDIR/S3_reg.txt;
> $OUTDIR/S3.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/S3.fasta;
  samtools faidx $REF -r $OUTDIR/S3_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_S3.vcf.gz -m $OUTDIR0/S3_mask.bed --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?" >> $OUTDIR/S3.fasta;
done;


### Stratum S2
mkdir $OUTDIR0/S2;
OUTDIR=$OUTDIR0/S2;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -R $OUTDIR0/S2.bed | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_S2.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_S2.vcf.gz;
echo scaffold_1:$(cat $OUTDIR0/S2.bed | cut -f2 | sort | head -n1)-$(cat $OUTDIR0/S2.bed | cut -f3 | sort | tail -n1) > $OUTDIR/S2_reg.txt;
> $OUTDIR/S2.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/S2.fasta;
  samtools faidx $REF -r $OUTDIR/S2_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_S2.vcf.gz -m $OUTDIR0/S2_mask.bed --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?" >> $OUTDIR/S2.fasta;
done;

### Stratum S1
mkdir $OUTDIR0/S1;
OUTDIR=$OUTDIR0/S1;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -R $OUTDIR0/S1.bed | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_S1.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_S1.vcf.gz;
echo scaffold_1:$(cat $OUTDIR0/S1.bed | cut -f2 | sort | head -n1)-$(cat $OUTDIR0/S1.bed | cut -f3 | sort | tail -n1) > $OUTDIR/S1_reg.txt;
> $OUTDIR/S1.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/S1.fasta;
  samtools faidx $REF -r $OUTDIR/S1_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_S1.vcf.gz -m $OUTDIR0/S1_mask.bed --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?" >> $OUTDIR/S1.fasta;
done;


### Stratum S0
mkdir $OUTDIR0/S0;
OUTDIR=$OUTDIR0/S0;
cd $OUTDIR;


# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -R $OUTDIR0/S0.bed | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_S0.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_S0.vcf.gz;
echo scaffold_1:$(cat $OUTDIR0/S0.bed | cut -f2 | sort | head -n1)-$(cat $OUTDIR0/S0.bed | cut -f3 | sort | tail -n1) > $OUTDIR/S0_reg.txt;
> $OUTDIR/S0.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/S0.fasta;
  samtools faidx $REF -r $OUTDIR/S0_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_S0.vcf.gz -m $OUTDIR0/S0_mask.bed --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?" >> $OUTDIR/S0.fasta;
done;

### Stratum 4A
mkdir $OUTDIR0/4A;
OUTDIR=$OUTDIR0/4A;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -r scaffold_1:128450000-137750000 | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_4A.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_4A.vcf.gz;
echo scaffold_1:128450000-137750000 > $OUTDIR/4A_reg.txt;
> $OUTDIR/4A.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/4A.fasta;
  samtools faidx $REF -r $OUTDIR/4A_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_4A.vcf.gz -m $MISSING --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?" >> $OUTDIR/4A.fasta;
done;



### Stratum 3a
mkdir $OUTDIR0/3a;
OUTDIR=$OUTDIR0/3a;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -r scaffold_1:137750000-146200000 | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_3a.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_3a.vcf.gz;
echo scaffold_1:137750000-146200000 > $OUTDIR/3a_reg.txt;
> $OUTDIR/3a.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/3a.fasta;
  samtools faidx $REF -r $OUTDIR/3a_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_3a.vcf.gz -m $MISSING --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?" >> $OUTDIR/3a.fasta;
done;

### Stratum 3b
mkdir $OUTDIR0/3b;
OUTDIR=$OUTDIR0/3b;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -r scaffold_1:146200000-149750000 | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_3b.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_3b.vcf.gz;
echo scaffold_1:146200000-149750000 > $OUTDIR/3b_reg.txt;
> $OUTDIR/3b.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/3b.fasta;
  samtools faidx $REF -r $OUTDIR/3b_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_3b.vcf.gz -m $MISSING --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?"  >> $OUTDIR/3b.fasta;
done;


### Stratum 3c
mkdir $OUTDIR0/3c;
OUTDIR=$OUTDIR0/3c;
cd $OUTDIR;

# Skylark
bcftools view $OUTDIR0/Skylark_filtered_snps.vcf.gz -r scaffold_1:149750000-213950000 | \
bgzip -c > $OUTDIR/Skylark_filtered_snps_3c.vcf.gz;
tabix $OUTDIR/Skylark_filtered_snps_3c.vcf.gz;
echo scaffold_1:149750000-213950000 > $OUTDIR/3c_reg.txt;
> $OUTDIR/3c.fasta;

for sample in $(bcftools query -l $OUTDIR0/Skylark_filtered_snps.vcf.gz); do
  echo ">$sample" >> $OUTDIR/3c.fasta;
  samtools faidx $REF -r $OUTDIR/3c_reg.txt | \
  bcftools consensus -s $sample -f $REF -H I $OUTDIR/Skylark_filtered_snps_3c.vcf.gz -m $MISSING --mask-with "?" | \
  grep -v "^>" | tr -d "\n" | awk '{print  $0"\n"}' | tr -d "?"  >> $OUTDIR/3c.fasta;
done;

### Make trees
conda deactivate;
conda activate tree_reconstruction;

# 5
OUTDIR=$OUTDIR0/5;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/5.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# S3
OUTDIR=$OUTDIR0/S3;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/S3.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# S2
OUTDIR=$OUTDIR0/S2;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/S2.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# S1
OUTDIR=$OUTDIR0/S1;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/S1.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# S0
#OUTDIR=$OUTDIR0/S0;
#cd $OUTDIR;

# Skylark
#iqtree2 -s $OUTDIR/S0.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# 4A
OUTDIR=$OUTDIR0/4A;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/4A.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# 3a
OUTDIR=$OUTDIR0/3a;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/3a.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# 3b
OUTDIR=$OUTDIR0/3b;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/3b.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;


# 3c
OUTDIR=$OUTDIR0/3c;
cd $OUTDIR;

# Skylark
iqtree2 -s $OUTDIR/3c.fasta -m GTR+F+G -nt AUTO -B 1000 -wbtl;
