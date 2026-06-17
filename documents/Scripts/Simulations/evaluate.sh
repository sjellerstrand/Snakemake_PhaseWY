#!/bin/bash -l

#SBATCH -A naiss2025-5-344
#SBATCH -p shared
#SBATCH -n 5
#SBATCH -t 15:00:00
#SBATCH -J evaluate
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

### Create folders
mkdir $OUTDIR0/evaluate;
OUTDIR0=$OUTDIR0/evaluate;
mkdir $OUTDIR0/runs \
$OUTDIR0/metadata;

# Set up replicates info
echo -e "Seed\tPopulation size\tGenerationt\tWhatsHap\tMAC\tWindow size\tStep size\tModel\tPairs" \
> $OUTDIR0/metadata/replicates.tsv;

# WhatsHap ON
for WINDOW in $(echo 100 1000 10000); do
  # MAC 1 & Pairs 10
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","1",WINDOW,WINDOW/4,"Hamming",10}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 5 & Pairs 10
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","5",WINDOW,WINDOW/4,"Hamming",10}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 9 & Pairs 10
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","9",WINDOW,WINDOW/4,"Hamming",10}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 1 & Pairs 8
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","1",WINDOW,WINDOW/4,"Hamming",8}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 1 & Pairs 6
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","1",WINDOW,WINDOW/4,"Hamming",6}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 1 & Pairs 4
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","1",WINDOW,WINDOW/4,"Hamming",4}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 1 & Pairs 2
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","1",WINDOW,WINDOW/4,"Hamming",2}' \
  >> $OUTDIR0/metadata/replicates.tsv;
done;
# WhatsHap OFF
for WINDOW in $(echo 100 1000 10000); do
  # MAC 1
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"OFF","1",WINDOW,WINDOW/4,"Hamming", 10}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 5
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"OFF","5",WINDOW,WINDOW/4,"Hamming",10}' \
  >> $OUTDIR0/metadata/replicates.tsv;
  # MAC 9
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"OFF","9",WINDOW,WINDOW/4,"Hamming",10}' \
  >> $OUTDIR0/metadata/replicates.tsv;
done;
# Inverse MAF model
for WINDOW in $(echo 100 1000 10000); do
  cat $WORKDIR/simulate_reads/metadata/replicates.txt | awk -F'\t' -v OFS='\t' -v WINDOW=$WINDOW '{print $1,$2,$3,"ON","1",WINDOW,WINDOW/4,"Inverse_MAF",10}' \
  >> $OUTDIR0/metadata/replicates.tsv;
done;



# Set up output file
echo -e "Simulation\tSeed\tPopulation size\tGeneration\tWhatsHap\tMAC\tWindow size\tStep size\tModel\tPairs\tScaffold\tDropout A\tDropout sex\tSensitivity A\tSensitivity X\tSensitivity Y1\tSensitivity Y2\tPrecision A\tPrecision X\tPrecision Y1\tPrecision Y2\tPi PAR\tPi Sex\tHeterozygosity A\tHeterozygosity sex" \
> $OUTDIR0/simulation_info.tsv;

for i in $(seq 1 1 $(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | wc -l)); do

  ### Settings
  SEED=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f1);
  POPSIZE=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f2);
  GEN=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f3);
  WHATSHAP=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f4);
  MAC=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f5);
  WINDOW=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f6);
  STEP=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f7);
  MODEL=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f8);
  PAIRS=$(cat $OUTDIR0/metadata/replicates.tsv | tail -n+2 | head -n$i | tail -n1 | cut -f9);
  
  ### Info
  SIMULATION=$(echo seed_$SEED\_popsize_$POPSIZE);
  SCAFFOLD=$(echo seed_$SEED\_popsize_$POPSIZE\_gen_$GEN);
  PIPAR=$(cat $WORKDIR/slim_sim/runs/seed_$SEED\_popsize_$POPSIZE/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_pi.txt | head -n2 | tail -n1 |cut -f3);
  PISEX=$(cat $WORKDIR/slim_sim/runs/seed_$SEED\_popsize_$POPSIZE/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_pi.txt | tail -n1 | cut -f3);
  
  ### Input
  VCF_TRUTH_IN=$WORKDIR/slim_sim/runs/seed_$SEED\_popsize_$POPSIZE/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.vcf.gz;
  HETGAM_DROPOUT=$WORKDIR/Phase_WY_results/sims_pairs_$PAIRS\_sex_depth_0.75_whatshap_$WHATSHAP/mac_$MAC\_window_$WINDOW\_step_$STEP\_model_$MODEL/beds/hetgam_dropout.bed;
  PHASE_INFO=$WORKDIR/Phase_WY_results/sims_pairs_$PAIRS\_sex_depth_0.75_whatshap_$WHATSHAP/mac_$MAC\_window_$WINDOW\_step_$STEP\_model_$MODEL/beds/phase_windows.bed;
  HETGAM=$(cat $METADATA/sample_info.txt | grep "Male" | head -n $PAIRS | tr '\n' ',' | head -c-1);

  # Make directory
  mkdir $OUTDIR0/runs/$SIMULATION\_gen_$GEN\_pairs_$PAIRS\_sex_depth_0.75_whatshap_$WHATSHAP\_mac_$MAC\_window_$WINDOW\_step_$STEP\_model_$MODEL;
  OUTDIR=$OUTDIR0/runs/$SIMULATION\_gen_$GEN\_pairs_$PAIRS\_sex_depth_0.75_whatshap_$WHATSHAP\_mac_$MAC\_window_$WINDOW\_step_$STEP\_model_$MODEL;
  
  ### Extract A from truth set and update reference allele
  bcftools view $VCF_TRUTH_IN \
  -s $HETGAM -r $SCAFFOLD:1-100000 | \
  bcftools norm -f $REF --check-ref s | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/truth_set_A.vcf.gz;
  VCF_TRUTH=$OUTDIR/truth_set_A.vcf.gz;
  tabix $VCF_TRUTH;

  ### Extract male A from analysed data
  bcftools view $WORKDIR/Phase_WY_results/sims_pairs_$PAIRS\_sex_depth_0.75_whatshap_$WHATSHAP/mac_$MAC\_window_$WINDOW\_step_$STEP\_model_$MODEL/vcfs/autosomal_filtered.vcf.gz \
  -s $HETGAM -r $SCAFFOLD:1-100000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_A.vcf.gz;
  VCF_IN=$OUTDIR/analysed_set_A.vcf.gz;
  tabix $VCF_IN;

  # Check A accuracy
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH -H | wc -l) -eq 0 ]; then
    SENSITIVITY_A=1;
    PRECISION_A=1;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH $VCF_IN | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$4+$5+$6}} END {print count" "total}');
    SENSITIVITY_A=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_A=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_A=0;
    PRECISION_A=NA;
  fi;

  ### Extract male X from truth set and update reference allele
  bcftools view $VCF_TRUTH_IN \
  -s $HETGAM -r $SCAFFOLD:100001-200000 | \
  bcftools norm -f $REF --check-ref s | \
  awk -F'\t' -v OFS='\t' '{if($0 ~ /^#/) {print; next} else \
  {for (i=10; i<=NF; i++) {split($i, alleles, "|"); $i = alleles[1]}; print}}' | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/truth_set_X.vcf.gz;
  VCF_TRUTH=$OUTDIR/truth_set_X.vcf.gz;
  tabix $VCF_TRUTH;

  ### Extract male X from analysed data
  bcftools view $WORKDIR/Phase_WY_results/sims_pairs_$PAIRS\_sex_depth_0.75_whatshap_$WHATSHAP/mac_$MAC\_window_$WINDOW\_step_$STEP\_model_$MODEL/vcfs/homogametic_filtered.vcf.gz \
  -s $HETGAM -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_X.vcf.gz;
  VCF_IN=$OUTDIR/analysed_set_X.vcf.gz;
  tabix $VCF_IN;

  # Check X accuracy
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH -H | wc -l) -eq 0 ]; then
    SENSITIVITY_X=NA;
    PRECISION_X=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH $VCF_IN | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_X=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_X=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_X=0;
    PRECISION_X=NA;
  fi;
  
  ### Extract male Y from truth set and update reference allele
  bcftools view $VCF_TRUTH_IN \
  -s $HETGAM -r $SCAFFOLD:100001-200000 | \
  bcftools norm -f $REF --check-ref s | \
  awk -F'\t' -v OFS='\t' '{if($0 ~ /^#/) {print; next} else \
  {for (i=10; i<=NF; i++) {split($i, alleles, "|"); $i = alleles[2]}; print}}' | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/truth_set_Y.vcf.gz;
  VCF_TRUTH=$OUTDIR/truth_set_Y.vcf.gz;
  tabix $VCF_TRUTH;

  ### Extract male Y from analysed data
  bcftools view $WORKDIR/Phase_WY_results/sims_pairs_$PAIRS\_sex_depth_0.75_whatshap_$WHATSHAP/mac_$MAC\_window_$WINDOW\_step_$STEP\_model_$MODEL/vcfs/heterogametic_filtered.vcf.gz \
  -s $HETGAM -r $SCAFFOLD:100001-200000 | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/analysed_set_Y.vcf.gz;
  VCF_IN=$OUTDIR/analysed_set_Y.vcf.gz;
  tabix $VCF_IN;

  # Check Y accuracy
  INFERRED=$(bcftools query -f '[%GT\n]' $VCF_IN | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view $VCF_TRUTH -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y1=NA;
    PRECISION_Y1=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -s - $VCF_TRUTH $VCF_IN | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y1=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y1=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y1=0;
    PRECISION_Y1=NA;
  fi;
  
  # Check Y accuracy without heterogametic dropouts
  INFERRED=$(bcftools view -T ^$HETGAM_DROPOUT $VCF_IN | bcftools query -f '[%GT\n]' | grep -vE '^(\./\.|\.\|\.)$' | wc -l);
  if [ $(bcftools view -T ^$HETGAM_DROPOUT $VCF_TRUTH -H | wc -l) -eq 0 ]; then
    SENSITIVITY_Y2=NA;
    PRECISION_Y2=NA;
  elif [ $INFERRED -gt 0 ]; then
    RESULTS=$(bcftools stats -T ^$HETGAM_DROPOUT -s - $VCF_TRUTH $VCF_IN | awk 'BEGIN {count=0; total=0} {if($1=="GCTs" || $1=="GCTi") {count+=$3+$9+$15}; if($1=="PSC" && $2!=1) {total+=$12+$13}} END {print count" "total}');
    SENSITIVITY_Y2=$(echo $RESULTS | awk '{print $1/$2}');
    PRECISION_Y2=$(echo $INFERRED $RESULTS | awk '{print $2/$1}');
  else
    SENSITIVITY_Y2=0;
    PRECISION_Y2=NA;
  fi;

  # Get proportion hetgam dropout
  cat $HETGAM_DROPOUT  | awk -F'\t' -v SCAFFOLD=$SCAFFOLD '{if($1==SCAFFOLD) {print}}' \
  > $OUTDIR/dropout.bed;
  
  DROPOUT_A=$(cat $OUTDIR/dropout.bed | awk -F'\t' 'BEGIN{SUM=0} {if($2 < 100000 && $3 > 100000) {SUM+=100000-$2} else if($3 <= 100000) {SUM+=$3-$2}} END {print SUM/100000}');
  DROPOUT_SEX=$(cat $OUTDIR/dropout.bed | awk -F'\t' 'BEGIN{SUM=0} {if($2 < 100000 && $3 > 100000) {SUM+=$3-100000} else if($2 >= 100000) {SUM+=$3-$2}} END {print SUM/100000}');

  # Get sex differences in heterozygosity
  cat $PHASE_INFO | awk -F'\t' -v SCAFFOLD=$SCAFFOLD '{if($1==SCAFFOLD) {print}}' \
  > $OUTDIR/phase_info.bed;
  
  if [ $(cat $OUTDIR/phase_info.bed | wc -l) -eq 0 ]; then
    HET_A=1;
    HET_SEX=1;
  else
    HET_A=$(cat $OUTDIR/phase_info.bed | awk -F'\t' 'BEGIN{SUM=0; N=0} {if($3 <= 100000 && $22 != NaN && $22 != "Inf") {SUM+=$22; N+=1}} END {if(SUM == 0 || N == 0) {print 1} else {print SUM/N}}');
    HET_SEX=$(cat $OUTDIR/phase_info.bed | awk -F'\t' 'BEGIN{SUM=0; N=0} {if($2 >= 100000 && $22 != NaN && $22 != "Inf") {SUM+=$22; N+=1}} END {if(SUM == 0 || N == 0) {print 1} else {print SUM/N}}');
  fi;

  # Print data
  echo -e "${SIMULATION}\t${SEED}\t${POPSIZE}\t${GEN}\t${WHATSHAP}\t${MAC}\t${WINDOW}\t${STEP}\t${MODEL}\t${PAIRS}\t${SCAFFOLD}\t${DROPOUT_A}\t${DROPOUT_SEX}\t${SENSITIVITY_A}\t${SENSITIVITY_X}\t${SENSITIVITY_Y1}\t${SENSITIVITY_Y2}\t${PRECISION_A}\t${PRECISION_X}\t${PRECISION_Y1}\t${PRECISION_Y2}\t${PIPAR}\t${PISEX}\t${HET_A}\t${HET_SEX}";
  
  rm -r $OUTDIR;
  
done >> $OUTDIR0/simulation_info.tsv;
