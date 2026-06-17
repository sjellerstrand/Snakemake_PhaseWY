#!/bin/bash -l

#SBATCH -A naiss2025-22-1246
#SBATCH -p shared
#SBATCH -n 40
#SBATCH -t 100:00:00
#SBATCH -J Genome_scans_windows
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Pairwise statistics are performed between populations, not sex

# Set parameters
SEX_SYSTEM=ZW
NONDIPLOIDS=No

## Fixed number of callable sites for windows
WINDOW=10000
STEP=$WINDOW

## Distance from exons
EXON_DIST=20000;

## Minimum callable sites per genome type (autosomal, homogametic, heterogametic) for scaffold to be analysed
MIN_SCAFFOLD_SIZE=0

MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=larks_2026;
WORKDIR=$MAINDIR/data/$PROJECT;
VCFS=$WORKDIR/Rasolark_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/vcfs;
BEDS=$WORKDIR/Rasolark_sex_depth_0.75_whatshap_ON/mac_1_window_10000_step_2500_model_Hamming/beds;
EXONS=$WORKDIR/B3_annotation_lift_over/*_exonic.bed;
OUTDIR0=$MAINDIR/working/$PROJECT;
REF=$MAINDIR/data/reference/Alauda_arvensis\
/Alauda_arvensis_M_hifiasm-purged-default_hap0.purged_no_mito.yahs_r2.scf.FINAL_mito.fasta;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT/analyses;
genomics_general=/cfs/klemming/projects/supr/snic2020-2-25/user_data/simon/bin_general/genomics_general-accessed-2025-11-10;

conda activate filter_variants;

# Load modules
ml PDC R;

# Define functions
genome_scans_windows=$FUNCTIONS/genome_scans_windows.r;
merge_scan_data=$FUNCTIONS/merge_scan_data.r;

# Create folders
mkdir $OUTDIR0/E2a_genome_scans_Rasolark_W;
OUTDIR0=$OUTDIR0/E2a_genome_scans_Rasolark_W;
mkdir $OUTDIR0/Reuseable_files \
$OUTDIR0/windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST \
$OUTDIR0/windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST/windows;
OUTDIR1=$OUTDIR0/windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST;

# Perform window based scans across the genome
DATASETS=W;

## Run through datasets
for DATA in $DATASETS; do

  ### Create subfolders
  mkdir $OUTDIR1/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST;
  OUTDIR2=$OUTDIR1/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST;

  ### Prepare BED files
  if [[ $EXON_DIST =~ ^[0-9]+$ ]]; then
    bedtools slop -i $EXONS -g $REF.fai -b $EXON_DIST | bedtools merge \
    > $OUTDIR0/Reuseable_files/exons_with_flank_$EXON_DIST.bed;
    bedtools subtract -a $BEDS/heterogametic.bed \
    -b $OUTDIR0/Reuseable_files/exons_with_flank_$EXON_DIST.bed | \
    awk -F'\t' '{if($1 == "scaffold_1") print}' \
    > $OUTDIR0/Reuseable_files/$DATA\_callable_$EXON_DIST.bed;
  else
    cat $BEDS/heterogametic.bed | awk -F'\t' '{if($1 == "scaffold_1") print}' \
    > $OUTDIR0/Reuseable_files/$DATA\_callable_$EXON_DIST.bed;
  fi;
  MASK=$OUTDIR0/Reuseable_files/$DATA\_callable_$EXON_DIST.bed;

  if [ $(find $OUTDIR1 | grep "$OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST\_infile.txt" | wc -l) == 0 ]; then

    ### Calculate windows with consideration to missing data
    Rscript $genome_scans_windows --args WINDOW=$WINDOW STEP=$STEP MIN_SCAFFOLD_SIZE=$MIN_SCAFFOLD_SIZE \
    EXON_DIST=$EXON_DIST PROJECT=$PROJECT MASK=$MASK REF=$REF OUTDIR1=$OUTDIR1 OUTDIR2=$OUTDIR2 DATA=$DATA;
    cat $OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST.txt | tail -n+2 | cut -f1,2,3 \
    > $OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST\_infile.txt;

  fi;

  ### Calculate populations statistics

  #### Get ploidy info
  if [ $DATA == "autosomal" ]; then

    if [ $NONDIPLOIDS == "Yes" ]; then
      # Get ploidy file from metadata folder
      ploidy=$METADATA/ploidy.txt;
    else
      # Set all individuals as diploid
      cat $METADATA/Rasolark_sample_info.txt | tail -n +2 | awk -F'\t' '{print $1"\t2"}' \
      > $OUTDIR2/ploidy.txt;
      ploidy=$OUTDIR2/ploidy.txt;
    fi;
  else
    if [ $SEX_SYSTEM != NA ]; then

      if [ $SEX_SYSTEM == ZW ]; then
        HOMGAM=Male;
        HETGAM=Female;
      elif [ $SEX_SYSTEM == XY ]; then
        HOMGAM=Female;
        HETGAM=Male;
      fi;
      if [ $DATA == "Z" ]; then
        ### Set all homogametes as diploids and all heterogametes as haploid
        cat $METADATA/Rasolark_sample_info.txt | awk -F'\t' -v HETGAM=$HETGAM '$3==HETGAM {print $1"\t1"}' \
        > $OUTDIR2/ploidy1.txt;
        cat $METADATA/Rasolark_sample_info.txt | awk -F'\t' -v HOMGAM=$HOMGAM '$3==HOMGAM {print $1"\t2"}' \
        > $OUTDIR2/ploidy2.txt;
        cat $OUTDIR2/ploidy1.txt $OUTDIR2/ploidy2.txt > $OUTDIR2/ploidy.txt;
        rm $OUTDIR2/ploidy1.txt $OUTDIR2/ploidy2.txt;
        ploidy=$OUTDIR2/ploidy.txt;

      elif [ $DATA == "W" ]; then
        ### Set all heterogametes as haploid
        cat $METADATA/Rasolark_sample_info.txt | awk -F'\t' -v HETGAM=$HETGAM '$3==HETGAM {print $1"\t1"}' \
        > $OUTDIR2/ploidy.txt;
        ploidy=$OUTDIR2/ploidy.txt;
      fi;
    fi;
  fi;

  #### Convert input vcf to correct format
  if [ $(find $OUTDIR0 | grep "$OUTDIR0/Reuseable_files/Rasolark_$DATA\_$EXON_DIST.geno.gz" | wc -l) == 0 ]; then
    bcftools view $VCFS/heterogametic_filtered.vcf.gz -R $MASK | \
    vcfclassify - | \
    vcffilter -s -f "!( INS | DEL | MNP )" -f "AC > 0 & AF < 1" | \
    vcftools --vcf - \
    --max-missing 0.95 \
    --max-alleles 2 \
    --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIR0/Reuseable_files/Rasolark_$DATA\_snps.vcf.gz;
    python3 $genomics_general/VCF_processing/parseVCF.py -i $OUTDIR0/Reuseable_files/Rasolark_$DATA\_snps.vcf.gz \
    --ploidyFile $ploidy | bgzip -c > $OUTDIR0/Reuseable_files/Rasolark\_$DATA\_$EXON_DIST.geno.gz;
  fi;

  #### Set window parameters
  PARAMETERS=$(echo --windType predefined --windCoords $OUTDIR1/windows/$DATA\_windows_$WINDOW\_steps_$STEP\_exon_dist_$EXON_DIST\_infile.txt);

  #### Calculate Pi and Tajima's D
  python3 $genomics_general/popgenWindows.py -g $OUTDIR0/Reuseable_files/Rasolark_$DATA\_$EXON_DIST.geno.gz \
 -f phased --writeFailedWindows --ploidyFile $ploidy -o $OUTDIR2/Rasolark_$DATA.pi_tajimas_D \
 --analysis popDist popFreq $PARAMETERS;

done;
