#!/bin/bash
#SBATCH --job-name=snowcrab_trim
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/snow_trim_%a.out
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -t 3-00:00:00
#SBATCH --array=1-63%7

module load bio/fastqc
source /home/lspencer/venv/bin/activate

IN=/share/afsc/snowcrab-OA-2022/raw-data/concat
OUT=/scratch/lspencer/snowcrab-OA-2022/trimmed
FASTQC=/scratch/lspencer/snowcrab-OA-2022/fastqc/trimmed

mkdir -p $OUT $FASTQC

# Get sample name from array index
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" /home/lspencer/snowcrab-OA-2022/scripts/sample_list.txt)

# Trim adapters and low-quality bases
cutadapt \
    -o ${OUT}/${SAMPLE}.trimmed.R1.fastq.gz \
    -p ${OUT}/${SAMPLE}.trimmed.R2.fastq.gz \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20,15 \
    -m 50 \
    --trim-n \
    --cores=10 \
    ${IN}/${SAMPLE}_R1.fastq.gz \
    ${IN}/${SAMPLE}_R2.fastq.gz \
    &> ${OUT}/cutadapt.${SAMPLE}.log

# Run FastQC
fastqc \
    --threads 10 \
    -o ${FASTQC} \
    ${OUT}/${SAMPLE}.trimmed.R1.fastq.gz \
    ${OUT}/${SAMPLE}.trimmed.R2.fastq.gz \
    &> ${FASTQC}/fastqc.${SAMPLE}.log
