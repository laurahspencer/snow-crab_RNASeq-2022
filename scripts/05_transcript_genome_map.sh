#!/bin/bash
#SBATCH --job-name=map_transcriptome
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/transcriptome_to_genome_map.out
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00

# Load environment (so gmap is in PATH)
source ~/.bashrc

# Set paths
GENOME_FA=/share/afsc/tanner-assembly/assembly/tanner.asm.hic.p_ctg.fa
GENOME_NAME="tanner"
TRANSCRIPTS_FA="/home/lspencer/snowcrab-OA-2022/trinity/Trinity-GG.fasta"                    # your genome-guided transcriptome fasta
OUT_DIR="/home/lspencer/snowcrab-OA-2022/trinity/gmap_alignments"
GMAP_INDEX_DIR="/home/lspencer/snowcrab-OA-2022/trinity/gmap_index"

# 1. Build GMAP genome index
mkdir -p "$GMAP_INDEX_DIR"
gmap_build -D "$GMAP_INDEX_DIR" -d "$GENOME_NAME" "$GENOME_FA"

# 2. Align Trinity transcripts to genome
mkdir -p "$OUT_DIR"
gmap -D "$GMAP_INDEX_DIR" -d "$GENOME_NAME" \
  -f gff3_gene --no-chimeras -n 0 -t 8 \
  "$TRANSCRIPTS_FA" > "$OUT_DIR/Trinity-GG.gff3"

# Done! Your output is now: gmap_alignments/Trinity-GG.gff3
