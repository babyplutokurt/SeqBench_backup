#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N post_hoc_0_6_0
#PBS -q normal
#PBS -l nodes=1:ppn=24
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_6_0_output.log
#PBS -e /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_6_0_error.log

#PBS -W depend=afterok:1888929:1888945:1888930


# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Step 1: Index the reference fasta file (if not already indexed)
if [ ! -f /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna.bwt ]; then
    bwa index /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna
fi

# Step 2: Alignment
bwa mem -t 48 /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/work/work/DecompressedOutput/ERR103405_1.fastq_-Q_20.fqz.fastq /home/tus53997/work/work/DecompressedOutput/ERR103405_2.fastq_-Q_20.fqz.fastq > /home/tus53997/work/work/SAM/ERR103405_1.fastq_-Q_20.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /home/tus53997/work/work/SAM/ERR103405_1.fastq_-Q_20.sam | samtools sort -@ 48 -o /home/tus53997/work/work/BAM/ERR103405_1.fastq_-Q_20_sorted.bam
samtools index --threads 48 /home/tus53997/work/work/BAM/ERR103405_1.fastq_-Q_20_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 48 -f /home/tus53997/work/work/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/work/work/BAM/ERR103405_1.fastq_-Q_20_sorted.bam | bcftools call -mv -Ov -o /home/tus53997/work/work/VCF/ERR103405_1.fastq_-Q_20.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 48 -c /home/tus53997/work/work/VCF/ERR103405_1.fastq_-Q_20.vcf > /home/tus53997/work/work/VCF/ERR103405_1.fastq_-Q_20.vcf.gz
tabix --threads 48 -p vcf /home/tus53997/work/work/VCF/ERR103405_1.fastq_-Q_20.vcf.gz

# Step 6: Compare the new VCF file with the original VCF file
bcftools isec --threads 48  -p /home/tus53997/work/work/VCF/comparison/0_6_0 -Oz /home/tus53997/work/work/VCF/ERR103405_1.fastq.vcf.gz /home/tus53997/work/work/VCF/ERR103405_1.fastq_-Q_20.vcf.gz

# Step 7: Calculate precision, recall, and F1 score
TP=$(bcftools view /home/tus53997/work/work/VCF/comparison/0_6_0/0002.vcf.gz | grep -v '^#' | wc -l)
FP=$(bcftools view /home/tus53997/work/work/VCF/comparison/0_6_0/0001.vcf.gz | grep -v '^#' | wc -l)
FN=$(bcftools view /home/tus53997/work/work/VCF/comparison/0_6_0/0000.vcf.gz | grep -v '^#' | wc -l)

precision=$(echo "scale=6; $TP / ($TP + $FP)" | bc)
recall=$(echo "scale=6; $TP / ($TP + $FN)" | bc)
f1_score=$(echo "scale=6; 2 * ($precision * $recall) / ($precision + $recall)" | bc)

echo "Precision: $precision"
echo "Recall: $recall"
echo "F1 Score: $f1_score"


echo "post_hoc_0_6_0,fqzcomp_-Q_20,$TP, $FP, $FN, $precision, $recall, $f1_score" >> "/home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/metrics/compression_metrics_ERR103405_1.fastq.csv"

# Deactivate the conda environment
conda deactivate