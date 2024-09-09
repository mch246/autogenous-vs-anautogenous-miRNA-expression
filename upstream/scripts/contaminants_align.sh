#!/bin/bash
#SBATCH --job-name=contaminants_align --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=mch284@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This runs bowtie2 for contaminant mapping#
#---------------------------------------------------------------------#
 
#- Load Module--------------------------------------------#
module load bowtie2/2.4.4
 
#- DEFINE FILE LOCATIONS--------------------------------------------#
 
index_path=/home/mch284/genomes/AalbF3/contaminant_index/contaminant
input_dir=/home/mch284/albopictus_auto_miRNA/trim_dir
filtered_dir=/home/mch284/albopictus_auto_miRNA/filtered_dir
sam_dir=/home/mch284/albopictus_auto_miRNA/sam_dir
 
#- RUN command ----------------#
files=(${input_dir}/*_cln.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _cln.fastq.gz`

bowtie2 -x ${index_path} -U ${input_dir}/${base}_cln.fastq.gz -S ${sam_dir}/${base}_flt.cln.sam --un-gz ${filtered_dir}/${base}_flt.fastq.gz

done 

#- Unload module----------------#
module unload bowtie2/2.4.4

#- FIN -----------------------#
