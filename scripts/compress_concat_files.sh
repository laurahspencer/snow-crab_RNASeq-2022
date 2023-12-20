#!/bin/bash

#SBATCH --job-name=compress_snow
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/compress.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10
#SBATCH -t 1000


IN=/share/afsc/snowcrab-OA-2022/raw-data/concat

FILES=$(ls ${IN}/*.fastq | awk -F "/" '{print $NF}')

for file in ${FILES}
do
    echo ${file}
    gzip -c ${IN}/${file} > ${IN}/${file}.gz
done
