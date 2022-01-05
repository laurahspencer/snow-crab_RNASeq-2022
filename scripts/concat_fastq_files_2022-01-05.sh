#!/bin/bash

#SBATCH --job-name=snowcrab_concat
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/snowcrab_concat.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 4
#SBATCH -t 1000

# This script is for concatenating two data files from the same sample
# but different lanes into one .fastq file.
# It has been slightly adapted from code written by Giles Goetz

BASE=/home/lspencer/snowcrab-OA-2022/raw-data
IN=${BASE}
OUT=${BASE}/concat

FOLDERS=$(ls -1d ${IN}/5* | awk -F "/" '{ print $NF }')

for folder in ${FOLDERS}
do
    echo ${folder}
    SUB_FOLDER=${IN}/${folder}
    SAMPLES=$(ls ${SUB_FOLDER}/5*_R1_*.fastq.gz | \
        awk -F "/" '{ print $NF }' | \
        awk -F "_" '{ print $2 }')

    for sample in ${SAMPLES}
    do
        echo ${sample}
        zcat ${SUB_FOLDER}/${folder}_${sample}_*_R1_*.fastq.gz \
            >> ${OUT}/sample_${sample}_R1.fastq
        zcat ${SUB_FOLDER}/${folder}_${sample}_*_R2_*.fastq.gz \
            >> ${OUT}/sample_${sample}_R2.fastq
    done
done
