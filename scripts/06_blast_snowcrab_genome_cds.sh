#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/diamond-blast_snow_genome.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0


module load bio/diamond

# Make DB once (Swiss-Prot fasta)
#diamond makedb --in /home/lspencer/references/blast/uniprot_sprot_20250812.fasta \
#               -d /home/lspencer/references/blast/uniprot_sprot_20250812

# Search
diamond blastx \
  -q /share/afsc/assemblies/snowcrab/GCA_016584305.1_ASM1658430v1_cds_from_genomic.fna \
  -d /share/afsc/assemblies/blast/uniprot_sprot_20250812.dmnd \
  -o /share/afsc/assemblies/snowcrab/GCA_016584305.1_ASM1658430v1_cds_diamond_blast.tab \
  -p 20 \
  --evalue 1e-10 \
  --max-target-seqs 1 \
  --sensitive \
  --outfmt 6


# This script is for blasting the snow crab juvenile transcriptome sequences against
# the Uniqprot/Swissprot database

module load bio/blast/2.15.0+

# Blast genes against uniprot/swissprot
blastx \
-query /share/afsc/assemblies/snowcrab/GCA_016584305.1_ASM1658430v1_cds_from_genomic.fna \
-db /share/afsc/assemblies/blast/uniprot_sprot_20250812_protein \
-out /share/afsc/assemblies/snowcrab/GCA_016584305.1_ASM1658430v1_cds_blastx.tab \
-evalue 1E-10 \
-num_threads 20 \
-outfmt 6 #\
-max_target_seqs 1 \
-max_hsps 1
