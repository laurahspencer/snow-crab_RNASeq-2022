#!/bin/bash
#SBATCH --job-name=multiqc_trim
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/multiqc_trim.out
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH -c 4

source /home/lspencer/venv/bin/activate

multiqc /scratch/lspencer/snowcrab-OA-2022/fastqc/trimmed \
  --outdir /scratch/lspencer/snowcrab-OA-2022/fastqc/trimmed

