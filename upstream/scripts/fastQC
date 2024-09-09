#!/bin/bash
#SBATCH --job-name=fastqc --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mch284@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=1G

#-----------------------------------------------------------------------------#
# This script gives quality control of fastq files #
#-----------------------------------------------------------------------------#


#- RUN fastqc ----------------------------------------------------------------#

/home/mch284/FastQC/fastqc -o /home/mch284/albopictus_auto_miRNA/fastqc_dir /home/mch284/albopictus_auto_miRNA/trim_dir/*cln.fastq.gz

#- FIN -----------------------------------------------------------------------#
