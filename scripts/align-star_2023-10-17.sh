#!/bin/bash

#SBATCH --job-name=snowcrab-align
#SBATCH --account=srlab
#SBATCH --partition=srlab
#SBATCH --nodes=1
#SBATCH --time=14-00:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lhs3@uw.edu
#SBATCH --output=/gscratch/srlab/lhs3/snowcrab-OA-2022/snow_align.txt

REF=/gscratch/scrubbed/lhs3/tanner-assembly
IN=/gscratch/scrubbed/lhs3/snowcrab-OA-2022/trimmed
OUT=/gscratch/scrubbed/lhs3/snowcrab-OA-2022/aligned-star/
STAR=/gscratch/srlab/programs/STAR-2.7.6a/bin/Linux_x86_64/STAR

# ========= Build STAR genome index
# Use 20 threads, -runThreadN 20
# specify that I want to generate genome, --runMode genomeGenerate
# specify path to save STAR genome directory. Must already exist, --genomeDir
# specify path to genome, --genomeFastaFiles
# specify path to annotation file, --sjdbGTFfile
# specify length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database.
#    Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads.
#    My reads are 150bp, so I'll use 149 (Default is 100). --sjdbOverhang 149

${STAR} \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeDir ${REF}/STAR \
--genomeFastaFiles ${REF}/tanner.asm.hic.p_ctg.fa \
--sjdbGTFfile /gscratch/scrubbed/lhs3/snowcrab-OA-2022/final_annotation_LHS.gtf \
--sjdbOverhang 149 \
done

# ========= Run STAR alignment
# Use most of STAR's default settings, which include (but aren't limited to) ...
# Max number of multiple alignments is 10, if exceeded read is considered unmapped --outFilterMultimapNmax 10
# Min numer of bp overlap to assign read to a gene is 1nt
# Also ... count number of reads per gene while mapping, using --quantMode GeneCounts

# store sample names to variable
# Sample IDs are numbers 1:63
SAMPLES=$(seq 63)

# loop through sample names
for sample in ${SAMPLES}
do
echo "Started mapping ${sample}"
${STAR} \
--runThreadN 20 \
--genomeDir ${REF}/STAR \
--readFilesIn ${IN}/sample_${sample}.trimmed.R1.v1.fastq.gz ${IN}/sample_${sample}.trimmed.R2.v1.fastq.gz \
--readFilesCommand gunzip -c \
--outFilterMultimapNmax 50 \
--outFileNamePrefix ${OUT}/${sample}. \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outReadsUnmapped Fastx \
&> ${OUT}/star-gadmor.${sample}.log
done
