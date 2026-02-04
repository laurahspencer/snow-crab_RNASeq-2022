#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/multiqc_aligned_snowcrab_genome.out
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH -c 4

source /home/lspencer/venv/bin/activate

multiqc /scratch/lspencer/snowcrab-OA-2022/aligned_snow_star/ \
  --outdir /scratch/lspencer/snowcrab-OA-2022/aligned_snow_star/

