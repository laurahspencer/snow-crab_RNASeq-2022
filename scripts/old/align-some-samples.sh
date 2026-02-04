#!/bin/bash

#SBATCH --job-name=genotype-snow
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/genotype-20231111.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH -p himem
#SBATCH -t 30-0:0:0
#SBATCH --mem=1200GB

module load aligners/bowtie2/2.4.2
module load bio/gatk/4.2.0.0
module load bio/samtools/1.15.1
source /home/lspencer/venv/bin/activate
module load bio/vcftools/0.1.16

REF=/home/lspencer/references/tanner #tanner crab genome directory
INPUT=/scratch/lspencer/snowcrab-OA-2022/trimmed #trimmed reads
OUTPUT=/home/lspencer/snowcrab-OA-2022/genotype # destination for all output files

# Run this alignment script on few remaining files that didn't complete before due to space issues 
FILES=$(ls /scratch/lspencer/snowcrab-OA-2022/trimmed/*.trimmed.R1.v1.fastq.gz | tail -n 8)

## Align trimmed reads to concatenated genome
echo "Aligning reads"

# Run Bowtie2 over each RNASeq paired sample
for file in ${FILES}
do
sample="$(basename -a ${file} | cut -d "." -f 1)"
file_R1="${sample}.trimmed.R1.v1.fastq.gz"
file_R2="${sample}.trimmed.R2.v1.fastq.gz"
map_file="${sample}.bowtie.concat.sam"

# run Bowtie2 on each file
bowtie2 \
-x tanner.asm.hic.p_ctg_concat.fa \
--sensitive \
--threads 24 \
--no-unal \
-1 ${INPUT}/${file_R1} \
-2 ${INPUT}/${file_R2} \
-S ${OUTPUT}/${map_file}; \
done >> ${OUTPUT}/02-bowtieout.txt 2>&1
