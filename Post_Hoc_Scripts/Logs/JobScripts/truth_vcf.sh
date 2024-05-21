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
if [ ! -f /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna.bwt ]; then
    bwa index /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna
fi

# Step 2: Alignment
<<<<<<< HEAD
bwa mem -t 24 /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/SeqBench/FASTQ/ERR103405_1.fastq /home/tus53997/SeqBench/FASTQ/ERR103405_2.fastq > /home/tus53997/SeqBench/SAM/ERR103405_1.fastq.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /home/tus53997/SeqBench/SAM/ERR103405_1.fastq.sam | samtools sort -@ 24 -o /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_sorted.bam
samtools index --threads 24 /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 24 -f /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_sorted.bam | bcftools call -mv -Ov -o /home/tus53997/SeqBench/VCF/ERR103405_1.fastq.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 24 -c /home/tus53997/SeqBench/VCF/ERR103405_1.fastq.vcf > /home/tus53997/SeqBench/VCF/ERR103405_1.fastq.vcf.gz
tabix --threads 24 -p vcf /home/tus53997/SeqBench/VCF/ERR103405_1.fastq.vcf.gz
=======
bwa mem -t 48 /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/work/work/FASTQ/ERR103405_1.fastq /home/tus53997/work/work/FASTQ/ERR103405_2.fastq > /home/tus53997/work/work/SAM/ERR103405_1.fastq.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /home/tus53997/work/work/SAM/ERR103405_1.fastq.sam | samtools sort -@ 48 -o /home/tus53997/work/work/BAM/ERR103405_1.fastq_sorted.bam
samtools index --threads 48 /home/tus53997/work/work/BAM/ERR103405_1.fastq_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 48 -f /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/work/work/BAM/ERR103405_1.fastq_sorted.bam | bcftools call -mv -Ov -o /home/tus53997/work/work/VCF/ERR103405_1.fastq.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 48 -c /home/tus53997/work/work/VCF/ERR103405_1.fastq.vcf > /home/tus53997/work/work/VCF/ERR103405_1.fastq.vcf.gz
tabix --threads 48 -p vcf /home/tus53997/work/work/VCF/ERR103405_1.fastq.vcf.gz
>>>>>>> 1f0bd7d1b255f712c9ab85d63df20a581d9db948

conda deactivate