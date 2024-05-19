#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -N job_0_6_0
#PBS -l nodes=1:ppn=8
#PBS -M @
#PBS -o /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_6_0_output.log
#PBS -e /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_6_0_error.log

# change to directory where 'qsub' was called
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run compression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/fqzcomp/fqzcomp -Q 20 /home/tus53997/SeqBench/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq /home/tus53997/SeqBench/CompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_20.fqz
END_TIME=$SECONDS
COMPRESSION_DURATION=$((END_TIME - START_TIME))

INPUT_SIZE=$(stat -c %s "/home/tus53997/SeqBench/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq")
OUTPUT_SIZE=$(stat -c %s "/home/tus53997/SeqBench/CompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_20.fqz")
RATIO=$(echo "scale=6; $INPUT_SIZE / $OUTPUT_SIZE" | bc)

sleep 10

# Run decompression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/fqzcomp/fqzcomp -d /home/tus53997/SeqBench/CompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_20.fqz /home/tus53997/SeqBench/DecompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq_-Q_20.fqz.fastq
END_TIME=$SECONDS
DECOMPRESSION_DURATION=$((END_TIME - START_TIME))

echo "job_0_6_0,fqzcomp_-Q_20,$COMPRESSION_DURATION,$DECOMPRESSION_DURATION,$RATIO" >> "/home/tus53997/SeqBench/Compression_Scripts/Logs/metrics/compression_metrics_HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1.fastq.csv"

conda deactivate