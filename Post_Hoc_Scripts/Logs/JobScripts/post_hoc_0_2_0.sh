#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N post_hoc_0_2_0
#PBS -l nodes=1:ppn=24
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_2_0_output.log
#PBS -e /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_2_0_error.log
#PBS -W depend=afterok:72287:72288:72295


# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Step 1: Index the reference fasta file (if not already indexed)
if [ ! -f /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna.bwt ]; then
    bwa index /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna
fi

# Step 2: Alignment
bwa mem -t 24 /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/SeqBench/DecompressedOutput/ERR103405_1.fastq_-Q_5.fqz.fastq /home/tus53997/SeqBench/DecompressedOutput/ERR103405_2.fastq_-Q_5.fqz.fastq > /home/tus53997/SeqBench/SAM/ERR103405_1.fastq_-Q_5.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /home/tus53997/SeqBench/SAM/ERR103405_1.fastq_-Q_5.sam | samtools sort -@ 24 -o /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_-Q_5_sorted.bam
samtools index --threads 24 /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_-Q_5_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 24 -f /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_-Q_5_sorted.bam | bcftools call -mv -Ov -o /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-Q_5.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 24 -c /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-Q_5.vcf > /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-Q_5.vcf.gz
tabix --threads 24 -p vcf /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-Q_5.vcf.gz

# Step 6: Compare the new VCF file with the original VCF file
bcftools isec --threads 24  -p /home/tus53997/SeqBench/VCF/comparison/0_2_0 -Oz /home/tus53997/SeqBench/VCF/ERR103405_1.fastq.vcf.gz /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-Q_5.vcf.gz

# Step 7: Calculate precision, recall, and F1 score
TP=$(bcftools view /home/tus53997/SeqBench/VCF/comparison/0_2_0/0002.vcf.gz | grep -v '^#' | wc -l)
FP=$(bcftools view /home/tus53997/SeqBench/VCF/comparison/0_2_0/0001.vcf.gz | grep -v '^#' | wc -l)
FN=$(bcftools view /home/tus53997/SeqBench/VCF/comparison/0_2_0/0000.vcf.gz | grep -v '^#' | wc -l)

precision=$(echo "scale=6; $TP / ($TP + $FP)" | bc)
recall=$(echo "scale=6; $TP / ($TP + $FN)" | bc)
f1_score=$(echo "scale=6; 2 * ($precision * $recall) / ($precision + $recall)" | bc)

echo "Precision: $precision"
echo "Recall: $recall"
echo "F1 Score: $f1_score"


echo "post_hoc_0_2_0,fqzcomp_-Q_5,$TP, $FP, $FN, $precision, $recall, $f1_score" >> "/home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/metrics/compression_metrics_ERR103405_1.fastq.csv"

# Deactivate the conda environment
conda deactivate