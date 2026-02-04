#!/bin/bash

#SBATCH --job-name=snowcrab-trim
#SBATCH --account=srlab
#SBATCH --partition=srlab
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lhs3@uw.edu
#SBATCH --output=/gscratch/srlab/lhs3/snowcrab-OA-2022/snow_trim.txt

# This script is for trimming raw (but concatenated) RNA-Seq data and
# filtering for quality and length.
# It has been slightly adapted from code written by Giles Goetz

module load /sw/modules-1.775/modulefiles/contrib/fastqc/0.11.5
conda activate cutadaptenv

IN=/gscratch/srlab/lhs3/snowcrab-OA-2022/raw-data/concat
OUT=/gscratch/scrubbed/lhs3/snow-crab-OA-2022/trimmed
FASTQC=/gscratch/scrubbed/lhs3/snow-crab-OA-2022/fastqc/trimmed
VER=1

SAMPLES=$(ls ${IN}/*_R1.fastq.gz | \
awk -F "/" '{print $NF}' | \
awk -F "." '{print $1}' | \
sed -e 's/_R1//')

for sample in ${SAMPLES}
do
    # Trimming the Illumina adapters
    # Quality-trim  5’ end with cutoff=20 & 3’ end with cutoff=15
    # Trimming out leftover N's
    # Filtering out sequences shorter then 50bp
    cutadapt \
        -o ${OUT}/${sample}.trimmed.R1.v${VER}.fastq.gz \
        -p ${OUT}/${sample}.trimmed.R2.v${VER}.fastq.gz \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -q 20,15 \
        -m 50 \
        --trim-n \
        --cores=20 \
        ${IN}/${sample}_R1.fastq.gz \
        ${IN}/${sample}_R2.fastq.gz \
        &> ${OUT}/cutadapt.${sample}.v${VER}.log

# Run fastqc on trimmed data files
    fastqc \
        --threads 10 \
        -o ${FASTQC} \
        ${OUT}/${sample}.trimmed.R1.v${VER}.fastq.gz \
        ${OUT}/${sample}.trimmed.R2.v${VER}.fastq.gz \
        &> ${FASTQC}/fastqc.${sample}.v${VER}.log
done

# Run multiqc to summarize fastqc reports
/gscratch/srlab/programs/anaconda3/bin/multiqc \
${FASTQC} \
--outdir ${FASTQC}
