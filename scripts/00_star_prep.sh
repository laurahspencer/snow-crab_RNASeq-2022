#!/bin/bash
#SBATCH --job-name=genome_prep
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=16
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/snow_genome_star_prep.out

module load aligners/star
export SINGULARITY_BINDPATH="/share:/share"

#GENOME=/home/lspencer/references/tanner/tanner.asm.hic.p_ctg.fa
#INDEX_DIR=/home/lspencer/references/tanner/star_index
GENOME=/share/afsc/assemblies/snowcrab/GCA_016584305.1_ASM1658430v1_genomic.fna
GTF=/share/afsc/assemblies/snowcrab/GCA_016584305.1_ASM1658430v1_genomic.gtf
INDEX_DIR=/share/afsc/assemblies/snowcrab/star_index
mkdir -p $INDEX_DIR

## Convert GFF to GTF 
#/opt/bioinformatics/build/gffread-0.12.6.Linux_x86_64/gffread \
#/home/lspencer/references/tanner/final_annotation.gff -T -o /home/lspencer/references/tanner/final_annotation.gtf

#GTF=/home/lspencer/references/tanner/final_annotation.gtf

## Estimate read length from trimmed reads — assume 150 bp here
## Set sjdbOverhang to read length - 1
sjdbOverhang=149

# Generate STAR genome index with annotation
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir $INDEX_DIR \
     --genomeFastaFiles $GENOME \
     --sjdbGTFfile $GTF \
     --sjdbOverhang $sjdbOverhang
