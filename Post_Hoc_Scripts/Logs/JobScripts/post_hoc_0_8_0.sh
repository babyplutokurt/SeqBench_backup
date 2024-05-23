#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N post_hoc_0_8_0
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_8_0_output.log
#PBS -e /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_8_0_error.log

#PBS -W depend=afterok:1893784:1893783:1893815


# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Step 1: Index the reference fasta file (if not already indexed)
if [ ! -f /work/tus53997/RefSeq/hg38.fa.bwt ]; then
    bwa index /work/tus53997/RefSeq/hg38.fa
fi

# Step 2: Alignment
bwa mem -t 28 /work/tus53997/RefSeq/hg38.fa /work/tus53997/DecompressedOutput/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.spring.fastq /work/tus53997/DecompressedOutput/ERR016162_2.fastq_-c_-l_-q_qvz_0.4.spring.fastq > /work/tus53997/SAM/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /work/tus53997/SAM/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.sam | samtools sort -@ 28 -o /work/tus53997/BAM/ERR016162_1.fastq_-c_-l_-q_qvz_0.4_sorted.bam
samtools index --threads 28 /work/tus53997/BAM/ERR016162_1.fastq_-c_-l_-q_qvz_0.4_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 28 -f /work/tus53997/RefSeq/hg38.fa /work/tus53997/BAM/ERR016162_1.fastq_-c_-l_-q_qvz_0.4_sorted.bam | bcftools call -mv -Ov -o /work/tus53997/VCF/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 28 -c /work/tus53997/VCF/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.vcf > /work/tus53997/VCF/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.vcf.gz
tabix --threads 28 -p vcf /work/tus53997/VCF/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.vcf.gz

# Step 6: Compare the new VCF file with the original VCF file
bcftools isec --threads 28  -p /work/tus53997/VCF/comparison/0_8_0 -Oz /work/tus53997/VCF/ERR016162_1.fastq.vcf.gz /work/tus53997/VCF/ERR016162_1.fastq_-c_-l_-q_qvz_0.4.vcf.gz

# Step 7: Calculate precision, recall, and F1 score
TP=$(bcftools view /work/tus53997/VCF/comparison/0_8_0/0002.vcf.gz | grep -v '^#' | wc -l)
FP=$(bcftools view /work/tus53997/VCF/comparison/0_8_0/0001.vcf.gz | grep -v '^#' | wc -l)
FN=$(bcftools view /work/tus53997/VCF/comparison/0_8_0/0000.vcf.gz | grep -v '^#' | wc -l)

precision=$(echo "scale=6; $TP / ($TP + $FP)" | bc)
recall=$(echo "scale=6; $TP / ($TP + $FN)" | bc)
f1_score=$(echo "scale=6; 2 * ($precision * $recall) / ($precision + $recall)" | bc)

echo "Precision: $precision"
echo "Recall: $recall"
echo "F1 Score: $f1_score"


echo "post_hoc_0_8_0,spring_-c_-l_-q_qvz_0.4,$TP, $FP, $FN, $precision, $recall, $f1_score" >> "/home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/metrics/post_hoc_metrics_ERR016162_1.fastq.csv"

# Deactivate the conda environment
conda deactivate