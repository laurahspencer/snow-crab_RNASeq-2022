#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/align_snowgenome_star%a.out
#SBATCH --array=1-63%15
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00

module load bio/samtools/1.19
module load aligners/star/2.7.10a

GENOME=/home/lspencer/references/snowcrab/GCA_016584305.1_ASM1658430v1_genomic.fna
INDEX_DIR=/home/lspencer/references/snowcrab/star_index
READ_DIR=/scratch/lspencer/snowcrab-OA-2022/trimmed
OUT_DIR=/scratch/lspencer/snowcrab-OA-2022/aligned_snow_star

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" /home/lspencer/snowcrab-OA-2022/scripts/sample_list.txt)

STAR --runThreadN 8 \
     --genomeDir $INDEX_DIR \
     --readFilesIn $READ_DIR/${SAMPLE}.trimmed.R1.fastq.gz $READ_DIR/${SAMPLE}.trimmed.R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 40000000000 \
     --outSAMstrandField intronMotif \
     --outFileNamePrefix $OUT_DIR/${SAMPLE}_
