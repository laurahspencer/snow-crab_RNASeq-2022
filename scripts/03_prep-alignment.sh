#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/salmon_index.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

# create sample list 
READ_DIR=/scratch/lspencer/snowcrab-OA-2022/trimmed
ls ${READ_DIR}/*.trimmed.R1.fastq.gz \
  | sed 's|.*/||; s/\.trimmed\.R1\.fastq\.gz$//' \
  > /home/lspencer/snowcrab-OA-2022/scripts/sample_list.txt


# prep transcriptome for alignment

module load aligners/salmon/1.10.0

TRANSCRIPTS=/home/lspencer/snowcrab-OA-2022/trinity/Trinity-GG.fasta
INDEX=/home/lspencer/snowcrab-OA-2022/trinity/salmon_index_GG

mkdir -p $(dirname "$INDEX")
salmon index -t "$TRANSCRIPTS" -i "$INDEX" -p 8

