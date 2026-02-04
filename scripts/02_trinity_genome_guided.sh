#!/bin/bash
#SBATCH --job-name=trinity
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/trinity_genome_guided2.out
#SBATCH --cpus-per-task=24
#SBATCH --partition himem
#SBATCH --time=9-00:00:00

HOME=/home/lspencer/snowcrab-OA-2022
BAM_DIR=/scratch/lspencer/snowcrab-OA-2022/aligned_star
GENOME=/share/afsc/tanner-assembly/assembly/tanner.asm.hic.p_ctg.fa
OUTDIR=$HOME/trinity

module load bio/samtools/1.19

## Merge BAMs
#mkdir -p $OUTDIR/bam_merge
#samtools merge -@ 8 \
#  $OUTDIR/bam_merge/merged.sorted.bam \
#  $BAM_DIR/*Aligned.sortedByCoord.out.bam

# Index merged BAM
#samtools index $OUTDIR/bam_merge/merged.sorted.bam

source ~/.bashrc
conda activate trinity-2.14.0 
module load bio/jellyfish/2.3.0
module load aligners/bowtie2/2.4.2
module load aligners/salmon/1.4.0
module load bio/samtools/1.15.1

export TRINITY_HOME=/opt/bioinformatics/build/trinityrnaseq-v2.14.0
export PATH=/opt/bioinformatics/build/trinityrnaseq-v2.14.0:${PATH}

# Run Trinity genome-guided mode
Trinity \
  --genome_guided_bam $OUTDIR/bam_merge/merged.sorted.bam \
  --genome_guided_max_intron 10000 \
  --genome_guided_min_coverage 1 \
  --CPU 24 \
  --max_memory 1200G \
  --output $OUTDIR
