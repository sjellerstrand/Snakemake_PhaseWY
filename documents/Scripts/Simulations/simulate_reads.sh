#!/bin/bash -l

#SBATCH -A naiss2025-5-344
#SBATCH -p shared
#SBATCH -n 64
#SBATCH -t 02:00:00
#SBATCH -J simulate_reads
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

MINGEN=2;
MAXGEN=7;
GENERATIONS=$(for i in $(seq $MINGEN 1 $MAXGEN); do echo "10 ^ $i" | bc; done;);

### Set parameters
MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=simulations_2025;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR0=$MAINDIR/working/$PROJECT;
METADATA=$WORKDIR/metadata;

# Activate conda environment
conda activate simulate_reads;

### Create folders
mkdir $OUTDIR0/simulate_reads;
OUTDIR0=$OUTDIR0/simulate_reads;
mkdir $OUTDIR0/runs \
$OUTDIR0/metadata \
$OUTDIR0/alignments;
cd $OUTDIR0;

### Set up data
for GEN in $GENERATIONS; do
  cat $WORKDIR/slim_sim/metadata/replicates.txt | awk -F'\t' -v GEN=$GEN '{print $1"\t"$2"\t"GEN}'
done | sort -k2,2 -k3,3 -nr | \
awk -F'\t' '{if($2<1000000) {print} else if($3 < 10000000) {print}}' \
> $OUTDIR0/metadata/replicates.txt;

REF_IND=Ind_REF_Female;
cat $WORKDIR/slim_sim/metadata/sample_names.txt | grep -v $REF_IND \
> $OUTDIR0/metadata/sample_names.txt;


export WORKDIR OUTDIR0 METADATA REF_IND

## Run simulations in parallel
simulate_reads() {

  SEED=$1;
  POPSIZE=$2;
  GEN=$3;

  mkdir $OUTDIR0/runs/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN;
  OUTDIR=$OUTDIR0/runs/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN;
  mkdir $OUTDIR/reads \
  $OUTDIR/alignments;

  # index reference reads
  echo ">seed_${SEED}_popsize_${POPSIZE}_gen_${GEN}" \
  > $OUTDIR/alignments/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_reference.fa;
  cat $WORKDIR/slim_sim/runs/seed_$SEED\_popsize_$POPSIZE/gen_$GEN/Sequences/$REF_IND\_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN:0.fa | tail -n+2 \
  >> $OUTDIR/alignments/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_reference.fa;  
  bwa index $OUTDIR/alignments/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_reference.fa;

  # Simulate illumina reads
  for IND in $(cat $OUTDIR0/metadata/sample_names.txt); do
    ngsngs -i $WORKDIR/slim_sim/runs/seed_$SEED\_popsize_$POPSIZE/gen_$GEN/Sequences/$IND\_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN:0.fa \
    -c 10 -ld Norm,350,20 -seq PE -f fastq.gz -q1 $METADATA/AccFreqL150R1.txt -q2 $METADATA/AccFreqL150R2.txt \
    -s "$(echo $SEED$(echo $IND |cut -d'_' -f2)0)" -o $OUTDIR/reads/$IND\_0;
    ngsngs -i $WORKDIR/slim_sim/runs/seed_$SEED\_popsize_$POPSIZE/gen_$GEN/Sequences/$IND\_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN:1.fa \
    -c 10 -ld Norm,350,20 -seq PE -f fastq.gz -q1 $METADATA/AccFreqL150R1.txt -q2 $METADATA/AccFreqL150R2.txt\
    -s "$(echo $SEED$(echo $IND |cut -d'_' -f2)0)" -o $OUTDIR/reads/$IND\_1;
    
    # Concatenate haplotypes
    zcat $OUTDIR/reads/$IND\_0_R1.fq.gz $OUTDIR/reads/$IND\_1_R1.fq.gz | \
    gzip > $OUTDIR/reads/$IND\_R1.fq.gz;
    zcat $OUTDIR/reads/$IND\_0_R2.fq.gz $OUTDIR/reads/$IND\_1_R2.fq.gz | \
    gzip > $OUTDIR/reads/$IND\_R2.fq.gz;
    
    ## Remove temporary files;
    rm $OUTDIR/reads/$IND\_0_* \
    $OUTDIR/reads/$IND\_1_*;

    ## Align reads
    bwa mem -t 1 -M $OUTDIR/alignments/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_reference.fa \
    $OUTDIR/reads/$IND\_R1.fq.gz $OUTDIR/reads/$IND\_R2.fq.gz \
    -R "@RG\tID:$(echo $IND\_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN)\tSM:$IND\tLB:$(echo seed_$SEED\_popsize_$POPSIZE\_gen_$GEN)\tPU:$SEED\tPL:NGSNGS" | \
    samtools view -f2 -F260 -q20 -b -@ 1 | \
    samtools sort -@ 1 \
    > $OUTDIR/alignments/$IND.bam;

  done;

  cd $OUTDIR0;

};


## Excecute function in parallell
export -f simulate_reads;
parallel --colsep '\t' 'simulate_reads {}' :::: $OUTDIR0/metadata/replicates.txt;

## Merge reference sequences
cat $(find $OUTDIR0/runs/ -wholename "*_reference.fa" | awk -F'[_/]' '{for (i = 1; i <= NF; i++) { \
    if ($i == "seed") seed = $(i + 1); if ($i == "popsize") popsize = $(i + 1); if ($i == "gen") gen = $(i + 1);} \
    print seed, popsize, gen, $0;}' | sort -k1n -k2n -k3n | cut -d' ' -f4-) \
> $OUTDIR0/alignments/reference.fasta;
samtools faidx $OUTDIR0/alignments/reference.fasta;

## Merge bams per individual
for IND in $(cat $OUTDIR0/metadata/sample_names.txt); do
  BAMS=$(find $OUTDIR0/runs/ -wholename "*alignments/$IND.bam" | awk -F'[_/]' '{for (i = 1; i <= NF; i++) { \
    if ($i == "seed") seed = $(i + 1); if ($i == "popsize") popsize = $(i + 1); if ($i == "gen") gen = $(i + 1);} \
    print seed, popsize, gen, $0;}' | sort -k1n -k2n -k3n | cut -d' ' -f4-);
  samtools merge $OUTDIR0/alignments/$IND\_merged.bam -@ 120 $BAMS;
  samtools index $OUTDIR0/alignments/$IND\_merged.bam -@ 120;
  rm $BAMS;
done;

## Remove temporary files
rm -r $OUTDIR0/runs;
