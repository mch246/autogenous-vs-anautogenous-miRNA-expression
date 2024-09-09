#!/bin/bash
#SBATCH --job-name=miRDeep2 --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mch284@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs miRDeep2 #
#-----------------------------------------------------------------------------#

#- Go into Conda environment where miRDeep2.pl is  --------------------------------#

source activate mirDeep2

cd /home/mch284/albopictus_auto_miRNA/miRDeep_dir2/

#- Set variables ----------------------------------------------------------------#

genome=/home/mch284/genomes/AalbF3/aedes_albopictus_AalbF3.fa
reads=/home/mch284/albopictus_auto_miRNA/miRDeep_dir2/processed_reads.fa
arf=/home/mch284/albopictus_auto_miRNA/miRDeep_dir2/reads_vs_genome.arf
known_miRNA=/home/mch284/albopictus_auto_miRNA/miRNAs/mature_miRNAs_known_nospace.fa


#- RUN miRDeep2.pl ----------------------------------------------------------------#

miRDeep2.pl ${reads} ${genome} ${arf} ${known_miRNA} none none 2> miRDeep2_report.log


#- FIN -----------------------------------------------------------------------#
