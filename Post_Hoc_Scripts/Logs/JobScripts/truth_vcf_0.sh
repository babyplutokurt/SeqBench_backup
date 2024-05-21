#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N truth_vcf
#PBS -q normal
#PBS -l nodes=1:ppn=24
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
bwa mem -t 48 /work/tus53997/RefSeq/hg38.fa /work/tus53997/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq /work/tus53997/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2.fastq > /work/tus53997/SAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /work/tus53997/SAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.sam | samtools sort -@ 48 -o /work/tus53997/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_sorted.bam
samtools index --threads 48 /work/tus53997/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 48 -f /work/tus53997/RefSeq/hg38.fa /work/tus53997/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_sorted.bam | bcftools call -mv -Ov -o /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 48 -c /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf > /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf.gz
tabix --threads 48 -p vcf /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf.gz

conda deactivate