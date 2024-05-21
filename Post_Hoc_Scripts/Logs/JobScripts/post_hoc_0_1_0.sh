#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N post_hoc_0_1_0
#PBS -q normal
#PBS -l nodes=1:ppn=24
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_1_0_output.log
#PBS -e /home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/logs/job_0_1_0_error.log
<<<<<<< HEAD
#PBS -W depend=afterok:72286:72285:72295
=======

#PBS -W depend=afterok:1889110:1889109:1889135
>>>>>>> 1f0bd7d1b255f712c9ab85d63df20a581d9db948


# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Step 1: Index the reference fasta file (if not already indexed)
<<<<<<< HEAD
if [ ! -f /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna.bwt ]; then
    bwa index /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna
fi

# Step 2: Alignment
bwa mem -t 24 /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/SeqBench/DecompressedOutput/ERR103405_1.fastq_-t_24.renano.fastq /home/tus53997/SeqBench/DecompressedOutput/ERR103405_2.fastq_-t_24.renano.fastq > /home/tus53997/SeqBench/SAM/ERR103405_1.fastq_-t_24.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /home/tus53997/SeqBench/SAM/ERR103405_1.fastq_-t_24.sam | samtools sort -@ 24 -o /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_-t_24_sorted.bam
samtools index --threads 24 /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_-t_24_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 24 -f /home/tus53997/SeqBench/RefSeq/GCF_000013425.1_ASM1342v1_genomic.fna /home/tus53997/SeqBench/BAM/ERR103405_1.fastq_-t_24_sorted.bam | bcftools call -mv -Ov -o /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-t_24.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 24 -c /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-t_24.vcf > /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-t_24.vcf.gz
tabix --threads 24 -p vcf /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-t_24.vcf.gz

# Step 6: Compare the new VCF file with the original VCF file
bcftools isec --threads 24  -p /home/tus53997/SeqBench/VCF/comparison/0_1_0 -Oz /home/tus53997/SeqBench/VCF/ERR103405_1.fastq.vcf.gz /home/tus53997/SeqBench/VCF/ERR103405_1.fastq_-t_24.vcf.gz
=======
if [ ! -f /work/tus53997/RefSeq/hg38.fa.bwt ]; then
    bwa index /work/tus53997/RefSeq/hg38.fa
fi

# Step 2: Alignment
bwa mem -t 48 /work/tus53997/RefSeq/hg38.fa /work/tus53997/DecompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.fqz.fastq /work/tus53997/DecompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2.fastq_-Q_2.fqz.fastq > /work/tus53997/SAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -bS /work/tus53997/SAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.sam | samtools sort -@ 48 -o /work/tus53997/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2_sorted.bam
samtools index --threads 48 /work/tus53997/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2_sorted.bam

# Step 4: Call variants using bcftools
bcftools mpileup --threads 48 -f /work/tus53997/RefSeq/hg38.fa /work/tus53997/BAM/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2_sorted.bam | bcftools call -mv -Ov -o /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.vcf

# Step 5: Compress and index the new VCF file
bgzip --threads 48 -c /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.vcf > /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.vcf.gz
tabix --threads 48 -p vcf /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.vcf.gz

# Step 6: Compare the new VCF file with the original VCF file
bcftools isec --threads 48  -p /work/tus53997/VCF/comparison/0_1_0 -Oz /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.vcf.gz /work/tus53997/VCF/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_2.vcf.gz
>>>>>>> 1f0bd7d1b255f712c9ab85d63df20a581d9db948

# Step 7: Calculate precision, recall, and F1 score
TP=$(bcftools view /work/tus53997/VCF/comparison/0_1_0/0002.vcf.gz | grep -v '^#' | wc -l)
FP=$(bcftools view /work/tus53997/VCF/comparison/0_1_0/0001.vcf.gz | grep -v '^#' | wc -l)
FN=$(bcftools view /work/tus53997/VCF/comparison/0_1_0/0000.vcf.gz | grep -v '^#' | wc -l)

precision=$(echo "scale=6; $TP / ($TP + $FP)" | bc)
recall=$(echo "scale=6; $TP / ($TP + $FN)" | bc)
f1_score=$(echo "scale=6; 2 * ($precision * $recall) / ($precision + $recall)" | bc)

echo "Precision: $precision"
echo "Recall: $recall"
echo "F1 Score: $f1_score"


echo "post_hoc_0_1_0,renano_-t_24,$TP, $FP, $FN, $precision, $recall, $f1_score" >> "/home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/metrics/compression_metrics_ERR103405_1.fastq.csv"

# Deactivate the conda environment
conda deactivate