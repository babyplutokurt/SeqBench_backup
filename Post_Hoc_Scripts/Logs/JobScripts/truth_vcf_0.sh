#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N truth_vcf
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/truth_vcf_0_output.log
#PBS -e /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/truth_vcf_0_error.log



# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Step 1: Index the reference fasta file (if not already indexed)
if [ ! -f /work/tus53997/RefSeq/hg38.fa.bwt ]; then
    bwa index /work/tus53997/RefSeq/hg38.fa
fi

# Step 2: Alignment
bwa mem -t 28 /work/tus53997/RefSeq/hg38.fa /work/tus53997/FASTQ/ERR016162_1.fastq /work/tus53997/FASTQ/ERR016162_2.fastq > /work/tus53997/SAM/ERR016162_1.fastq.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /work/tus53997/SAM/ERR016162_1.fastq.sam | samtools sort -@ 28 -o /work/tus53997/BAM/ERR016162_1.fastq_sorted.bam
samtools index --threads 28 /work/tus53997/BAM/ERR016162_1.fastq_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 28 -f /work/tus53997/RefSeq/hg38.fa /work/tus53997/BAM/ERR016162_1.fastq_sorted.bam | bcftools call -mv -Ov -o /work/tus53997/VCF/ERR016162_1.fastq.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 28 -c /work/tus53997/VCF/ERR016162_1.fastq.vcf > /work/tus53997/VCF/ERR016162_1.fastq.vcf.gz
tabix --threads 28 -p vcf /work/tus53997/VCF/ERR016162_1.fastq.vcf.gz

conda deactivate