#!/bin/bash -l

#SBATCH -A naiss2025-5-344
#SBATCH -p shared
#SBATCH -n 5
#SBATCH -t 15:00:00
#SBATCH -J evaluate2
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=simulations_2025;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR0=$MAINDIR/working/$PROJECT;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT;
REF=$WORKDIR/simulate_reads/alignments/reference.fasta;

# Activate conda environment
conda activate filter_variants

# Phase pairs
phasing_female_male_pair=$MAINDIR/scripts/simulations_2025/process_data/phasing_female_male_pair.py;

### Create folders
mkdir $OUTDIR0/evaluate2;
OUTDIR0=$OUTDIR0/evaluate2;
mkdir $OUTDIR0/runs \
$OUTDIR0/metadata;

# Phase single pair with python script
python3 $phasing_female_male_pair -i $WORKDIR/call_filter_variants/simulations_2025_filtered1.vcf.gz \
-o $OUTDIR0/1_pair_phase.vcf;
bgzip $OUTDIR0/1_pair_phase.vcf;
tabix $OUTDIR0/1_pair_phase.vcf.gz;

# Extract Y
bcftools view $OUTDIR0/1_pair_phase.vcf.gz -s Ind_11_Male_W | \
awk -F'\t' -v OFS='\t' '{if($0 ~ /^#CHROM/) {print $1, $2, $3, $4, $5, $6, $7, $8, $9, "Ind_11_Male"} else print}' | \
bcftools annotate -x INFO | \
bgzip -c > $OUTDIR0/Y_1_pair_phase.vcf.gz;
tabix $OUTDIR0/Y_1_pair_phase.vcf.gz;

# Extract X
bcftools view $OUTDIR0/1_pair_phase.vcf.gz -s Ind_11_Male_Z | \
awk -F'\t' -v OFS='\t' '{if($0 ~ /^#CHROM/) {print $1, $2, $3, $4, $5, $6, $7, $8, $9, "Ind_11_Male"} else print}' | \
bcftools annotate -x INFO | \
bgzip -c > $OUTDIR0/X_1_pair_phase.vcf.gz;
tabix $OUTDIR0/X_1_pair_phase.vcf.gz;

# Set up replicates info
echo -e "Seed\tPopulation size\tGenerationt\tWhatsHap\tMAC\tWindow size\tStep size\tModel\tPairs" \
> $OUTDIR0/metadata/replicates.tsv;

# Replicates
cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON",10000,2500,"Hamming"}' \
>> $OUTDIR0/metadata/replicates.tsv;

# Set up output file
echo -e "Simulation\tSeed\tPopulation size\tGeneration\tMethod\tScaffold\tSensitivity X\tSensitivity Y1\tSensitivity Y2\tPrecision X\tPrecision Y1\tPrecision Y2" \
> $OUTDIR0/simulation_info_method.tsv;

for i in $(seq 1 1 $(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | wc -l)); do

  ### Settings
  SEED=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f1);
  POPSIZE=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f2);
  GEN=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f3);
  WHATSHAP=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f4);
  WINDOW=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f5);
  STEP=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f6);
  MODEL=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f7);
  
  ### Info
  SIMULATION=$(echo seed_$SEED\_popsize_$POPSIZE);
  SCAFFOLD=$(echo seed_$SEED\_popsize_$POPSIZE\_gen_$GEN);
  
  ### Input
  VCF_TRUTH_IN=$WORKDIR/slim_sim/runs/seed_$SEED\_popsize_$POPSIZE/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.vcf.gz;
  HETGAM_DROPOUT1=$WORKDIR/Phase_WY_results/sims_pairs_1_het_2_hom_sex_depth_0.75_whatshap_$WHATSHAP/mac_0_window_$WINDOW\_step_$STEP\_model_$MODEL/beds/hetgam_dropout.bed;
  HETGAM_DROPOUT2=$WORKDIR/Phase_WY_results/sims_pairs_2_het_1_hom_sex_depth_0.75_whatshap_$WHATSHAP/mac_1_window_$WINDOW\_step_$STEP\_model_$MODEL/beds/hetgam_dropout.bed;
  HETGAM_DROPOUT3=$WORKDIR/Phase_WY_results/sims_pairs_1_sex_depth_0.75_whatshap_$WHATSHAP/mac_0_window_$WINDOW\_step_$STEP\_model_$MODEL/beds/hetgam_dropout.bed;
  PHASE_INFO1=$WORKDIR/Phase_WY_results/sims_pairs_1_het_2_hom_sex_depth_0.75_whatshap_$WHATSHAP/mac_0_window_$WINDOW\_step_$STEP\_model_$MODEL/beds/phase_windows.bed;
  PHASE_INFO2=$WORKDIR/Phase_WY_results/sims_pairs_2_het_1_hom_sex_depth_0.75_whatshap_$WHATSHAP/mac_1_window_$WINDOW\_step_$STEP\_model_$MODEL/beds/phase_windows.bed;

  # Make directory
  mkdir $OUTDIR0/runs/$SIMULATION\_gen_$GEN;
  OUTDIR=$OUTDIR0/runs/$SIMULATION\_gen_$GEN;

  ### Extract male X from truth set and update reference allele
  bcftools view $VCF_TRUTH_IN \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  bcftools norm -f $REF --check-ref s | \
  awk -F'\t' -v OFS='\t' '{if($0 ~ /^#/) {print; next} else \
  {for (i=10; i<=NF; i++) {split($i, alleles, "|"); $i = alleles[1]}; print}}' | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/truth_set_X.vcf.gz;
  VCF_TRUTH_X=$OUTDIR/truth_set_X.vcf.gz;
  tabix $VCF_TRUTH_X;
  
  ### Extract male Y from truth set and update reference allele
  bcftools view $VCF_TRUTH_IN \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  bcftools norm -f $REF --check-ref s | \
  awk -F'\t' -v OFS='\t' '{if($0 ~ /^#/) {print; next} else \
  {for (i=10; i<=NF; i++) {split($i, alleles, "|"); $i = alleles[2]}; print}}' | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/truth_set_Y.vcf.gz;
  VCF_TRUTH_Y=$OUTDIR/truth_set_Y.vcf.gz;
  tabix $VCF_TRUTH_Y;
  
  ### PhaseWY 2 females and 1 male

  ### Extract male X from analysed data PhaseWY 2 females and 1 male
  bcftools view $WORKDIR/Phase_WY_results/sims_pairs_1_het_2_hom\_sex_depth_0.75_whatshap_$WHATSHAP/mac_0_window_$WINDOW\_step_$STEP\_model_$MODEL/vcfs/homogametic.vcf.gz \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_X_PhaseWY1.vcf.gz;
  VCF_IN_X=$OUTDIR/analysed_set_X_PhaseWY1.vcf.gz;
  tabix $VCF_IN_X;
  
  ### Extract male Y from analysed data PhaseWY 2 females and 1 male
  bcftools view $WORKDIR/Phase_WY_results/sims_pairs_1_het_2_hom_sex_depth_0.75_whatshap_$WHATSHAP/mac_0_window_$WINDOW\_step_$STEP\_model_$MODEL/vcfs/heterogametic.vcf.gz \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_Y_PhaseWY1.vcf.gz;
  VCF_IN_Y=$OUTDIR/analysed_set_Y_PhaseWY1.vcf.gz;
  tabix $VCF_IN_Y;

  # Check X Sensitivity PhaseWY 2 female and 1 males
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN_X | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_X -H | wc -l) -eq 0 ]; then
    SENSITIVITY_X_PhaseWY1=NA;
    PRECISION_X_PhaseWY1=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH_X $VCF_IN_X | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_X_PhaseWY1=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_X_PhaseWY1=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_X_PhaseWY1=0;
    PRECISION_X_PhaseWY1=NA;
  fi;

  # Check Y sensitivity PhaseWY 2 females and 1 male
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN_Y | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_Y -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y1_PhaseWY1=NA;
    PRECISION_Y1_PhaseWY1=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH_Y $VCF_IN_Y | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y1_PhaseWY1=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y1_PhaseWY1=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y1_PhaseWY1=0;
    PRECISION_Y1_PhaseWY1=NA;
  fi;

  # Check Y sensitivity without heterogametic dropouts PhaseWY 2 females and 1 male
    INFERRED=$(bcftools view $VCF_IN_Y -T ^$HETGAM_DROPOUT1 | bcftools query -f '[%GT\n]' | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_Y -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y2_PhaseWY1=NA;
    PRECISION_Y2_PhaseWY1=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -T ^$HETGAM_DROPOUT1 -s - $VCF_TRUTH_Y $VCF_IN_Y | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y2_PhaseWY1=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y2_PhaseWY1=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y2_PhaseWY1=0;
    PRECISION_Y2_PhaseWY1=NA;
  fi;
  
  
  ### PhaseWY 1 females and 2 male

  ### Extract male X from analysed data PhaseWY 1 female and 2 males
  bcftools view $WORKDIR/Phase_WY_results/sims_pairs_2_het_1_hom\_sex_depth_0.75_whatshap_$WHATSHAP/mac_1_window_$WINDOW\_step_$STEP\_model_$MODEL/vcfs/homogametic.vcf.gz \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_X_PhaseWY2.vcf.gz;
  VCF_IN_X=$OUTDIR/analysed_set_X_PhaseWY2.vcf.gz;
  tabix $VCF_IN_X;
  
  ### Extract male Y from analysed data PhaseWY 1 female and 2 males
  bcftools view $WORKDIR/Phase_WY_results/sims_pairs_2_het_1_hom_sex_depth_0.75_whatshap_$WHATSHAP/mac_1_window_$WINDOW\_step_$STEP\_model_$MODEL/vcfs/heterogametic.vcf.gz \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_Y_PhaseWY2.vcf.gz;
  VCF_IN_Y=$OUTDIR/analysed_set_Y_PhaseWY2.vcf.gz;
  tabix $VCF_IN_Y;
  
  # Check X sensitivity PhaseWY 1 females and 2 male
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN_X | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_X -H | wc -l) -eq 0 ]; then
    SENSITIVITY_X_PhaseWY2=NA;
    PRECISION_X_PhaseWY2=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH_X $VCF_IN_X | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_X_PhaseWY2=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_X_PhaseWY2=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_X_PhaseWY2=0;
    PRECISION_X_PhaseWY2=NA;
  fi;

  # Check Y sensitivity PhaseWY 1 females and 2 male
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN_Y | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_Y -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y1_PhaseWY2=NA;
    PRECISION_Y1_PhaseWY2=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH_Y $VCF_IN_Y | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y1_PhaseWY2=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y1_PhaseWY2=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y1_PhaseWY2=0;
    PRECISION_Y1_PhaseWY2=NA;
  fi;

  # Check Y sensitivity without heterogametic dropouts PhaseWY 1 females and 2 male
    INFERRED=$(bcftools view $VCF_IN_Y -T ^$HETGAM_DROPOUT2 | bcftools query -f '[%GT\n]' | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_Y -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y2_PhaseWY2=NA;
    PRECISION_Y2_PhaseWY2=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -T ^$HETGAM_DROPOUT2 -s - $VCF_TRUTH_Y $VCF_IN_Y | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y2_PhaseWY2=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y2_PhaseWY2=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y2_PhaseWY2=0;
    PRECISION_Y2_PhaseWY2=NA;
  fi;
  
  ### python script

  ### Extract male X from analysed data Python script
  bcftools view $OUTDIR0/X_1_pair_phase.vcf.gz \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_X_script.vcf.gz;
  VCF_IN_X=$OUTDIR/analysed_set_X_script.vcf.gz;
  tabix $VCF_IN_X;

  ### Extract male Y from analysed data Python script
  bcftools view $OUTDIR0/Y_1_pair_phase.vcf.gz \
  -s Ind_11_Male -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_Y_script.vcf.gz;
  VCF_IN_Y=$OUTDIR/analysed_set_Y_script.vcf.gz;
  tabix $VCF_IN_Y;
 
  # Check X sensitivity Python script
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN_X | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_X -H | wc -l) -eq 0 ]; then
    SENSITIVITY_X_script=NA;
    PRECISION_X_script=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH_X $VCF_IN_X | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_X_script=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_X_script=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_X_script=0;
    PRECISION_X_script=NA;
  fi;

  # Check Y sensitivity Python script
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN_Y | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_Y -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y1_script=NA;
    PRECISION_Y1_script=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH_Y $VCF_IN_Y | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y1_script=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y1_script=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y1_script=0;
    PRECISION_Y1_script=NA;
  fi;

  # Check Y sensitivity without heterogametic dropouts Python script
  INFERRED=$(bcftools view $VCF_IN_Y -T ^$HETGAM_DROPOUT3 | bcftools query -f '[%GT\n]' | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH_Y -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y2_script=NA;
    PRECISION_Y2_script=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -T ^$HETGAM_DROPOUT3 -s - $VCF_TRUTH_Y $VCF_IN_Y | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y2_script=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y2_script=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y2_script=0;
    PRECISION_Y2_script2=NA;
  fi;


  # Print data
  echo -e "${SIMULATION}\t${SEED}\t${POPSIZE}\t${GEN}\tPhaseWY 2 females and 1 male\t${SCAFFOLD}\t${SENSITIVITY_X_PhaseWY1}\t${SENSITIVITY_Y1_PhaseWY1}\t${SENSITIVITY_Y2_PhaseWY1}\t${PRECISION_X_PhaseWY1}\t${PRECISION_Y1_PhaseWY1}\t${PRECISION_Y2_PhaseWY1}";
  echo -e "${SIMULATION}\t${SEED}\t${POPSIZE}\t${GEN}\tPhaseWY 1 female and 2 males\t${SCAFFOLD}\t${SENSITIVITY_X_PhaseWY2}\t${SENSITIVITY_Y1_PhaseWY2}\t${SENSITIVITY_Y2_PhaseWY2}\t${PRECISION_X_PhaseWY2}\t${PRECISION_Y1_PhaseWY2}\t${PRECISION_Y2_PhaseWY2}";
  echo -e "${SIMULATION}\t${SEED}\t${POPSIZE}\t${GEN}\tPython script\t${SCAFFOLD}\t${SENSITIVITY_X_script}\t${SENSITIVITY_Y1_script}\t${SENSITIVITY_Y2_script}\t${PRECISION_X_script}\t${PRECISION_Y1_script}\t${PRECISION_Y2_script}";

  rm -r $OUTDIR;
  
done >> $OUTDIR0/simulation_info_method.tsv;
