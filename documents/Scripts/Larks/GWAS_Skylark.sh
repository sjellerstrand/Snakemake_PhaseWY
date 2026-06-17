#!/bin/bash -l

#SBATCH -A naiss2025-22-1246
#SBATCH -p shared
#SBATCH -n 40
#SBATCH -t 100:00:00
#SBATCH -J GWAS
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL


# Set parameters
MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=larks_2026;
WORKDIR=$MAINDIR/data/$PROJECT;
VCFS=$WORKDIR/Skylark_Europe_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/vcfs;
OUTDIR0=$MAINDIR/working/$PROJECT;
REF=$MAINDIR/data/reference/Alauda_arvensis\
/Alauda_arvensis_M_hifiasm-purged-default_hap0.purged_no_mito.yahs_r2.scf.FINAL_mito.fasta;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT/analyses;

### Load modules
conda activate GWAS;

### Create folders
mkdir $OUTDIR0/GWAS_Skylark;
OUTDIR0=$OUTDIR0/GWAS_Skylark;

# Filter data
bcftools view $VCFS/homogametic_filtered.vcf.gz | \
vcfclassify - | \
vcffilter -s -f "!( INS | DEL | MNP )" -f "AC > 0 & AF < 1" | \
vcftools --vcf - \
--max-missing 0.95 \
--max-alleles 2 \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR0/homogametic_filtered_snps.vcf.gz;
tabix $OUTDIR0/homogametic_filtered_snps.vcf.gz;

## Create bed files
plink --vcf $OUTDIR0/homogametic_filtered_snps.vcf.gz \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--vcf-half-call missing \
--make-bed --out $OUTDIR0/GWAS;

# Add sex to fam file
mv $OUTDIR0/GWAS.fam $OUTDIR0/GWAS_tmp.fam;
for i in $(seq 1 1 $(cat $OUTDIR0/GWAS_tmp.fam | wc -l)); do
    IND=$(cat $OUTDIR0/GWAS_tmp.fam | head -n $i | tail -n1 | cut -d' ' -f1);
    SEX=$(cat $METADATA/Skylark_sample_info.txt | awk -F'\t' -v IND=$IND '{if($1==IND) print $3}');
    echo $(cat $OUTDIR0/GWAS_tmp.fam | head -n $i | tail -n1 | awk -v SEX=$SEX '{if(SEX=="Female") {sex=1} else {sex=2}; print $1,$2,$3,$4,$5,sex}');
done > $OUTDIR0/GWAS.fam;

# Run GWAS
plink --bfile $OUTDIR0/GWAS --logistic --out $OUTDIR0/GWAS_sex --maf 0.05 --allow-extra-chr --allow-no-sex;
cat $OUTDIR0/GWAS_sex.assoc.logistic | tr -s ' ' | tr ' ' '\t' > $OUTDIR0/GWAS_Sex_assoc.logistic_tabs;

# Create kinship matrix
cd $OUTDIR0;
gemma -bfile $OUTDIR0/GWAS -gk 1 -o kinship;

# Run GWAS
gemma -bfile $OUTDIR0/GWAS -k $OUTDIR0/output/kinship.cXX.txt -maf 0.05 -lmm 4 -o GWAS_out;