#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --account=srlab
#SBATCH --partition=srlab
#SBATCH --nodes=1
#SBATCH --time=14-00:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lhs3@uw.edu
#SBATCH --output=/gscratch/srlab/lhs3/snowcrab-OA-2022/blast_tanner.txt

# This script is for blasting the tanner annotation against the
# the Uniqprot/Swissprot database

module load contrib/blast+/2.7.1
module load contrib/bedtools/2.28.0

BASE=/gscratch/scrubbed/lhs3

## Make blast database (only need to do this once!)
#makeblastdb \
#-in /gscratch/scrubbed/lhs3/blast/uniprot_sprot_20231025.fasta \
#-dbtype prot \
#-out uniprot_sprot_database_20231025

## Before blasting, extract sequences for just genes in fasta format
#bedtools getfasta \
#-fi ${BASE}/tanner-assembly/tanner.asm.hic.p_ctg.fa \
#-bed ${BASE}/tanner-assembly/tanner_genes.bed \
#-fo ${BASE}/tanner-assembly/tanner_genes.fasta

# Blast genes against uniprot/swissprot
blastx \
-query ${BASE}/tanner-assembly/tanner_genes.fasta \
-db ${BASE}/blast/uniprot_sprot_database_20231025 \
-out ${BASE}/tanner-assembly/tanner_genes_blastx.tab \
-evalue 1E-5 \
-num_threads 20 \
-outfmt 6
