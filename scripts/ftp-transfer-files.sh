#!/bin/bash

#SBATCH --job-name=ftp-transfer
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/20230722-ftp-transfer.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 1000

# This script is for copying ranner raw rnaseq files
# from Sedna to the Phase Genomics server via FTP

module load tools/lftp

lftp 34.208.80.201
login tanner-rna gJxwDxa*EmL$4#

echo "transferring raw data in 5010/"
mirror -R /share/afsc/snowcrab-OA-2022/raw-data/5010/

echo "transferring raw data in 5329/"
mirror -R /share/afsc/snowcrab-OA-2022/raw-data/5329/

echo "transferring wget-log"
mirror -R /share/afsc/snowcrab-OA-2022/raw-data/wget-log
