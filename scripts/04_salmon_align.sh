#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/salmon_quant_%a.out
#SBATCH --array=1-63%5
#SBATCH --cpus-per-task=8
#SBATCH --time=5-00:00:00

module load aligners/salmon/1.10.0

READ_DIR=/scratch/lspencer/snowcrab-OA-2022/trimmed
INDEX=/home/lspencer/snowcrab-OA-2022/trinity/salmon_index_GG
OUT_DIR=/home/lspencer/snowcrab-OA-2022/salmon
LIST=/home/lspencer/snowcrab-OA-2022/scripts/sample_list.txt
LIBTYPE=A   # set to A if unstranded, ISF if FR

mkdir -p "$OUT_DIR"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

salmon quant \
  -i "$INDEX" \
  -l "$LIBTYPE" \
  -1 "$READ_DIR/${SAMPLE}.trimmed.R1.fastq.gz" \
  -2 "$READ_DIR/${SAMPLE}.trimmed.R2.fastq.gz" \
  -p 8 \
  --validateMappings \
  --gcBias \
  --seqBias \
  -o "$OUT_DIR/${SAMPLE}"

