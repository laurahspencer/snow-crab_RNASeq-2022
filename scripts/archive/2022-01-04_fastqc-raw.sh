#!/bin/bash

#SBATCH --job-name=snow_fastqc-raw
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/snow_fastqc-raw.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10

# This script is for running FastQC and MultiQC on raw (but concatenated)
# red king crab RNA-Seq data

# Load modules and virtual environments
module load bio/fastqc
source /home/lspencer/venv/bin/activate

# run fastqc on each raw read file for data from lane "5010"
cd /home/lspencer/snowcrab-OA-2022/raw-data/5010/

fastqc \
--threads 8 \
*.fastq.gz \
--outdir /scratch/lspencer/snowcrab-OA-2022/fastqc/raw/

# run fastqc on each raw read file for data from lane "5329"
cd /home/lspencer/snowcrab-OA-2022/raw-data/5329/

fastqc \
--threads 10 \
*.fastq.gz \
--outdir /scratch/lspencer/snowcrab-OA-2022/fastqc/raw/


# Run multiqc to summarize fastqc reports
multiqc \
/scratch/lspencer/snowcrab-OA-2022/fastqc/raw/ \
--outdir /scratch/lspencer/snowcrab-OA-2022/fastqc/raw/
