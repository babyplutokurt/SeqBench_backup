#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -N truth_vcf
#PBS -l nodes=1:ppn=8
#PBS -M @
#PBS -o /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/truth_vcf_0_output.log
#PBS -e /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/truth_vcf_0_error.log


# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Step 1: Index the reference fasta file (if not already indexed)
if [ ! -f /home/tus53997/SeqBench/RefSeq/hg38.fa.bwt ]; then
    bwa index /home/tus53997/SeqBench/RefSeq/hg38.fa
fi

# Step 2: Alignment
bwa mem -t 8 /home/tus53997/SeqBench/RefSeq/hg38.fa /home/tus53997/SeqBench/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq /home/tus53997/SeqBench/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2.fastq > /home/tus53997/SeqBench/SAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /home/tus53997/SeqBench/SAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.sam | samtools sort -@ 8 -o /home/tus53997/SeqBench/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_sorted.bam
samtools index --threads 8 /home/tus53997/SeqBench/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 8 -f /home/tus53997/SeqBench/RefSeq/hg38.fa /home/tus53997/SeqBench/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_sorted.bam | bcftools call -mv -Ov -o /home/tus53997/SeqBench/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 8 -c /home/tus53997/SeqBench/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf > /home/tus53997/SeqBench/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf.gz
tabix --threads 8 -p vcf /home/tus53997/SeqBench/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf.gz

conda deactivate