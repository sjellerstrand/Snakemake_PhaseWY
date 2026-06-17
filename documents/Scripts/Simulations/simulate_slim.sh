#!/bin/bash -l

#SBATCH -A naiss2025-5-344
#SBATCH -p memory
#SBATCH --mem=880GB
#SBATCH -t 3-00:00:00
#SBATCH -J simulate_slim
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL


### --mem=1760GB

NFEMALES=10;
NMALES=10;
RECRATE=1e-8;
MUTRATE=1e-8;
REPLICATES=10;
MINGEN=2;
MAXGEN=7;
GENERATIONS=$(for i in $(seq $MINGEN 1 $MAXGEN); do echo "10 ^ $i" | bc; done;);
MINPOPSIZE=2;
MAXPOPSIZE=6;
POPULATIONSIZES=$(for i in $(seq $MAXPOPSIZE -1 $MINPOPSIZE); do echo "10 ^ $i" | bc; done;);
CHROMSIZE=200000;

### Set parameters
MAINDIR=/cfs/klemming/projects/snic/snic2020-2-25/user_data/simon/PhaseWY;
PROJECT=simulations_2025;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR0=$MAINDIR/working/$PROJECT;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT;

# Define functions
slim_sim=$FUNCTIONS/slim/PhaseWY_simulation.slim;
slim_sim_resume=$FUNCTIONS/slim/PhaseWY_simulation_resume.slim;
sim_mut=$FUNCTIONS/slim/simulate_mutations.py;

# Activate conda environment
conda activate SLiM;

### Create folders
mkdir $OUTDIR0/slim_sim;
OUTDIR0=$OUTDIR0/slim_sim;
mkdir $OUTDIR0/runs \
$OUTDIR0/metadata;
cd $OUTDIR0;

if [ ! -f $OUTDIR0/metadata/replicates.txt ]; then

  ### Set up metadata
  for SIZE in $POPULATIONSIZES; do
  	for i in $(seq 1 1 $REPLICATES); do
      echo -e "$((RANDOM % 100000))\t$SIZE";
  	done;
  done > $OUTDIR0/metadata/replicates.txt;
  
  
  echo Ind_REF_Female \
  > $OUTDIR0/metadata/sample_names.txt;
  REF_IND=Ind_REF_Female;
  for IND in $(seq 1 1 $NFEMALES); do
      echo Ind_$IND\_Female;
  done >> $OUTDIR0/metadata/sample_names.txt;
  for IND in $(seq $(echo $NFEMALES + 1 | bc) 1 $(echo $NFEMALES + $NMALES | bc)); do
      echo Ind_$IND\_Male;
  done >> $OUTDIR0/metadata/sample_names.txt;

fi;
SLIMGENS=$(echo \"$(echo $GENERATIONS | tr ' ' ',')\");

export OUTDIR0 NFEMALES NMALES REF_IND RECRATE MUTRATE MINGEN GENERATIONS SLIMGENS CHROMSIZE slim_sim slim_sim_resume sim_mut;


## Run simulations in parallel
slim_sim() {

  SEED=$1;
  POPSIZE=$2;
 
  OUTDIR=$OUTDIR0/runs/seed_$SEED\_popsize_$POPSIZE;

  FIRST_GEN=$(echo $GENERATIONS | tr ' ' '\n' | head -n1);
  LAST_GEN=$(echo $GENERATIONS | tr ' ' '\n' | tail -n1);
  
    
  if [ $POPSIZE -eq 1000000 ]; then
    GENERATIONS=$(for i in $(seq $MINGEN 1 6); do echo "10 ^ $i" | bc; done;);
    SLIMGENS=$(echo \"$(echo $GENERATIONS | tr ' ' ',')\");
  fi;
 
  # Check if simulation is not finished
  if [ ! -f $OUTDIR/gen_$LAST_GEN/Sequences/$REF_IND\_seed_$SEED\_popsize\_$POPSIZE\_gen_$LAST_GEN:0.fa ]; then
  
    # If simualtion has been initiated, continue on previous run.
    if [ -f $OUTDIR/gen_$FIRST_GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$FIRST_GEN.trees ]; then
    
      OUTDIR=$OUTDIR0/runs/seed_$SEED\_popsize_$POPSIZE;
      cd $OUTDIR;
      
      LAST_CHECK=$(find $OUTDIR/ -name "*.trees" | rev | cut -d"." -f2 | cut -d"_" -f1 | rev | sort -n | tail -n1);
      CHECK_POP=$(find $OUTDIR/ -name "*$LAST_CHECK.trees");

      slim -s $SEED -d seed=$SEED -d recRate=$RECRATE -d popSize=$POPSIZE -d chromSize=$CHROMSIZE -d generations=$SLIMGENS -d checkpoint=$LAST_CHECK -d check_pop="'$CHECK_POP'" $slim_sim_resume;
  
    # Start new simulation
    else

      mkdir $OUTDIR0/runs/seed_$SEED\_popsize_$POPSIZE;
      OUTDIR=$OUTDIR0/runs/seed_$SEED\_popsize_$POPSIZE;
      cd $OUTDIR;
    
    	for GEN in $GENERATIONS; do
          mkdir $OUTDIR/gen_$GEN \
          $OUTDIR/gen_$GEN/Sequences;
    	done;
    
      # Generate random ancestral nuclotide sequence
      echo ">seed_${SEED}_popsize_${POPSIZE}" \
      > anc_sequence_seed_$SEED\_popsize_$POPSIZE.fasta;
      awk -v SEED=$SEED -v CHROMSIZE=$CHROMSIZE 'BEGIN {srand(SEED); bases = "ACGT"; \
      for(i=0; i<CHROMSIZE; i++) {printf "%s", substr(bases, int(rand()*4)+1, 1);}print "";}' \
      >> anc_sequence_seed_$SEED\_popsize_$POPSIZE.fasta;

      slim -s $SEED -d seed=$SEED -d recRate=$RECRATE -d popSize=$POPSIZE -d chromSize=$CHROMSIZE -d generations=$SLIMGENS $slim_sim;

    fi;
  
    for GEN in $GENERATIONS; do
   
      # Simulate neutral mutations
      python3 $sim_mut -i $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.trees \
      -s $SEED -r $RECRATE -m $MUTRATE -N $POPSIZE -g $GEN \
      -o $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_overlaid.trees \
      -x $OUTDIR/anc_sequence_seed_$SEED\_popsize_$POPSIZE.fasta \
      -z $OUTDIR/gen_$GEN/Sequences/anc_sequence_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.fasta;
      samtools faidx $OUTDIR/gen_$GEN/Sequences/anc_sequence_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.fasta;
  
      # Sample random males and females
      FEMALES=$(python3 -m tskit individuals $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_overlaid.trees | \
      cut -f1,5 | awk -F"[\t,:]" '{if($13==0) print "tsk_"$1, $13}'  | shuf | head -n $(echo "$NFEMALES + 1" | bc) | cut -d' ' -f1);
      MALES=$(python3 -m tskit individuals $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_overlaid.trees | \
      cut -f1,5 | awk -F"[\t,:]" '{if($13==1) print "tsk_"$1, $13}'  | shuf | head -n $NMALES | cut -d' ' -f1);
      echo $FEMALES $MALES | tr ' ' '\n'\
      > $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_individuals.txt;
      INDIVIDUALS=$OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_individuals.txt;
  
      # Sample individuals and convert to vcf
      python3 -m tskit vcf --contig-id seed_$SEED\_popsize_$POPSIZE\_gen_$GEN \
      $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN\_overlaid.trees | \
      bcftools view -S $INDIVIDUALS | \
      bcftools reheader -s $OUTDIR0/metadata/sample_names.txt | \
      vcffixup - | \
      vcffilter -f "AC > 0" | \
      awk 'BEGIN{OFS="\t"} /^#/ {print $0; next} {$2 = $2 + 1; print}' | \
      bgzip -c > $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.vcf.gz;
      tabix $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.vcf.gz;
  
      # Extract sequences
      if [ $(bcftools view $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.vcf.gz -H | wc -l) -gt 0 ]; then
        cd $OUTDIR/gen_$GEN/Sequences;
        vcf2fasta -f $OUTDIR/gen_$GEN/Sequences/anc_sequence_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.fasta \
        $OUTDIR/gen_$GEN/seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.vcf.gz;
        cd $OUTDIR;
      else
        for IND in $(cat $OUTDIR0/metadata/sample_names.txt); do
         cp $OUTDIR/gen_$GEN/Sequences/anc_sequence_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.fasta \
         $OUTDIR/gen_$GEN/Sequences/$IND\_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN:0.fa;
         cp $OUTDIR/gen_$GEN/Sequences/anc_sequence_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN.fasta \
         $OUTDIR/gen_$GEN/Sequences/$IND\_seed_$SEED\_popsize_$POPSIZE\_gen_$GEN:1.fa;
        done;
      fi;
  
    done;
  
    cd $OUTDIR0;
    
  fi;

};


## Excecute function in parallell
export -f slim_sim;
parallel --colsep '\t' 'slim_sim {}' :::: $OUTDIR0/metadata/replicates2.txt;
