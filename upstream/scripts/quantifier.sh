#!/bin/bash
#SBATCH --job-name=quantifier_17Dec --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mch284@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs the quantifier module of miRDeep2 #
#-----------------------------------------------------------------------------#

#- Go into Conda environment where quantifier.pl is  --------------------------------#

source activate mirDeep2

cd /home/mch284/albopictus_auto_miRNA/quantified_miRNA_dir_17Dec/

#- Set variables ----------------------------------------------------------------#

reads=/home/mch284/albopictus_auto_miRNA/miRDeep_dir2/processed_reads.fa
mature_miRNA=/home/mch284/albopictus_auto_miRNA/miRNAs/FULL_mature_star_auto_vienna_downgrade_nospace.fa
precursor_miRNA=/home/mch284/albopictus_auto_miRNA/miRNAs/FULL_Precursor_auto_vienna_downgrade_nospace.fa
config_file=/home/mch284/albopictus_auto_miRNA/miRDeep_dir2/config_file.txt


#- RUN quantifier.pl ----------------------------------------------------------------#

quantifier.pl -p ${precursor_miRNA} -m ${mature_miRNA} -r ${reads} -c ${config_file} -d -k

#- FIN -----------------------------------------------------------------------#
