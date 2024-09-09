#!/bin/bash
#SBATCH --job-name=trim --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mch284@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs trimmomatic to clean up raw reads #
#-----------------------------------------------------------------------------#


#- Set variables --------------------------------------------------------------$

raw_dir=/home/mch284/albopictus_auto_miRNA/raw_data

trim_dir=/home/mch284/albopictus_auto_miRNA/trim_dir

trim=/home/mch284/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar

adapter=/home/mch284/trimmomatic/Trimmomatic-0.39/adapters/smRNA_NexFlex_adapters.fa:2:30:10

#- RUN Trimmomatic-------------------------------------------------------------$

files=(${raw_dir}/*_R1_001.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L002_R1_001.fastq.gz`
java -Xmx2G -jar ${trim} SE \
          ${raw_dir}/${base}_L002_R1_001.fastq.gz \
          ${trim_dir}/${base}_cln.fastq.gz \
          ILLUMINACLIP:${adapter} \
          HEADCROP:4 \
          TRAILING:6 \
          SLIDINGWINDOW:4:15 \
          MINLEN:17
done
