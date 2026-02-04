#!/bin/bash

#SBATCH --job-name=concat-genome
#SBATCH --output=/home/lspencer/snowcrab-OA-2022/sbatch_logs/concat-genome.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH -p himem
#SBATCH -t 30-0:0:0
#SBATCH --mem=1200GB

module load aligners/bowtie2/2.4.2
module load bio/gatk/4.2.0.0
module load bio/samtools/1.15.1
source /home/lspencer/venv/bin/activate
module load bio/vcftools/0.1.16

REF=/home/lspencer/references/tanner #tanner crab genome directory

# Move to the directory containing tanner crab genome.
# Bowtie2 looks for index in current directory - easier just to be in it from the get-go
cd ${REF}

########### -------------- CREATE CONCATENATED RED KING CRAB GENOME WITH ~100 SUPER-CONTIGS -------------------############################
# This step is to speed up this job. Otherwise it takes weeks/months/years!

# Check to see if there are instances of 1000 N in the genome already - this returns the number of contigs that contain exactly 1000 consecutive N's
echo "This is the number of contigs that contain exactly 1000 consecutive N's"
grep N tanner.asm.hic.p_ctg.fa | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | grep 1000$ | wc -l

# Concatenate contigs together, where i=the number of sequences to combine +1
echo "Starting genome concatenation"
awk -v i=247 -f /home/lspencer/snowcrab-OA-2022/scripts/merge_contigs.awk.txt tanner.asm.hic.p_ctg.fa > tanner.asm.hic.p_ctg_concat_temp.fa
sed 's/> ptg000001l/>ptg000001l/' tanner.asm.hic.p_ctg_concat_temp.fa > tanner.asm.hic.p_ctg_concat.fa  # remove the leading character before ptg000001l

## Edit fasta headers from concatenated IDs to the range (e.g. "Seq1 Seq2 Seq3" to "Seq1:Seq3")
# create variables
n_contigs_old=$(grep ">" tanner.asm.hic.p_ctg.fa | wc -l) #number of contigs in original fasta
d="${n_contigs_old//[^[:digit:]]/}" #Return ID of last/highest number contig
digits=$(echo ${#d}) #number of digits in the last contig
n_contigs_new=$(grep ">" tanner.asm.hic.p_ctg_concat.fa | wc -l) #number of contigs in concatenated genome
n=$(seq 1 1 $n_contigs_new) #vector of 1:# of contigs in concat. genome

# loop over each new contig header, paste first and last ID # together and replace concat header echo with new ID
for i in $n
do
old=$(grep ">" tanner.asm.hic.p_ctg_concat.fa | sed "${i}q;d")
start=$(echo $old | grep -o -E 'ptg[0-9]{1,6}l' | head -n 1) #If needed, change the 6 to match the value in $digits
start_old=$(echo ">"${start})
end=$(echo $old | grep -o -E 'ptg[0-9]{1,6}l' | tail -n 1) #If needed, change the 6 to match the value in $digits
new=$(echo ">"$start":"$end)
echo $new >> tmp.txt
sed -i "s|"$start_old".*|$new|" tanner.asm.hic.p_ctg_concat.fa # edit fasta in place.
done

# Index concatenated genome
echo "Indexing concatenated and re-headed fasta"
samtools faidx tanner.asm.hic.p_ctg_concat.fa

