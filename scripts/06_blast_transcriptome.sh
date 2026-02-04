#!/bin/bash

#SBATCH --job-name=diam-blast
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/diamond-blast_transcriptome.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0


module load bio/diamond   # or use a conda env with diamond

# Make DB once (Swiss-Prot fasta)
diamond makedb --in /home/lspencer/references/blast/uniprot_sprot_20250812.fasta \
               -d /home/lspencer/references/blast/uniprot_sprot_20250812

# Search
diamond blastx \
  -q /home/lspencer/snowcrab-OA-2022/trinity/Trinity-GG.fasta \
  -d /home/lspencer/references/blast/uniprot_sprot_20250812.dmnd \
  -o /home/lspencer/snowcrab-OA-2022/trinity/transcriptome_diamond_blast.tab \
  -p 20 \
  --evalue 1e-10 \
  --max-target-seqs 1 \
  --sensitive \
  --outfmt 6


# NOTE: old code, I abandoned BLASTX because it was taking WAY too long!  
# This script is for blasting the snow crab juvenile transcriptome sequences against
# the Uniqprot/Swissprot database

#module load bio/blast/2.15.0+

## Blast genes against uniprot/swissprot
#blastx \
#-query /home/lspencer/snowcrab-OA-2022/trinity/Trinity-GG.fasta \
#-db /home/lspencer/references/blast/uniprot_sprot_20250812_protein \
#-out /home/lspencer/snowcrab-OA-2022/trinity/transcriptome_blastx.tab \
#-evalue 1E-5 \
#-num_threads 20 \
#-outfmt 6 #\
#-max_target_seqs 1 \
#-max_hsps 1
