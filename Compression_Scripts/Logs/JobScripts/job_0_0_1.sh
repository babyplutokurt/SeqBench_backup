#!/bin/sh
#PBS -l walltime=168:00:00
#PBS -N job_0_0_1
#PBS -l nodes=1:ppn=6
#PBS -l mem=128gb
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_0_1_output.log
#PBS -e /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_0_1_error.log

# change to directory where 'qsub' was called
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run compression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/SZ3/bin/sz3 -f -1 37500000 -M REL 0.1 -i /home/tus53997/SeqBench/binFile/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2_1000000.bin -z /home/tus53997/SeqBench/CompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2_1000000.bin_-f_-1_37500000_-M_REL_0.1.sz
END_TIME=$SECONDS
COMPRESSION_DURATION=$((END_TIME - START_TIME))

INPUT_SIZE=$(stat -c %s "/home/tus53997/SeqBench/binFile/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2_1000000.bin")
OUTPUT_SIZE=$(stat -c %s "/home/tus53997/SeqBench/CompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2_1000000.bin_-f_-1_37500000_-M_REL_0.1.sz")
RATIO=$(echo "scale=2; $INPUT_SIZE / $OUTPUT_SIZE" | bc)

# Run decompression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/SZ3/bin/sz3 -f -1 37500000 -M REL 0.1 -z /home/tus53997/SeqBench/CompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2_1000000.bin_-f_-1_37500000_-M_REL_0.1.sz -o /home/tus53997/SeqBench/DecompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2_1000000.bin_-f_-1_37500000_-M_REL_0.1.sz.fastq
END_TIME=$SECONDS
DECOMPRESSION_DURATION=$((END_TIME - START_TIME))

echo "job_0_0_1,SZ3_-f_-1_37500000_-M_REL_0.1,$COMPRESSION_DURATION,$DECOMPRESSION_DURATION,$RATIO" >> "/home/tus53997/SeqBench/Compression_Scripts/Logs/metrics/compression_metrics_HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R2_1000000.fastq.csv"

conda deactivate