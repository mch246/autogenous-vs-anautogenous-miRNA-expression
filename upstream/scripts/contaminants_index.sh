#!/bin/bash
#SBATCH --job-name=contaminants_index --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mch284@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs makes a genome index with bowtie #
#-----------------------------------------------------------------------------#

module load bowtie2/2.4.4

#- Set variables ----------------------------------------------------------------#

ref_genome=/home/mch284/genomes/AalbF3/albo_contaminants.fasta
index_path_prefix=/home/mch284/genomes/AalbF3/contaminant_index/contaminant


#- RUN bowtie ----------------------------------------------------------------#

bowtie2-build ${ref_genome} ${index_path_prefix}


#- Unload module----------------#
module unload bowtie2/2.4.4

#- FIN -----------------------------------------------------------------------#
