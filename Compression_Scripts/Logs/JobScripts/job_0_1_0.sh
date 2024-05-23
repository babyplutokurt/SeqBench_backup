#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N job_0_1_0
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_1_0_output.log
#PBS -e /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_1_0_error.log

# change to directory where 'qsub' was called
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run compression and capture metrics
START_TIME=$SECONDS

/home/tus53997/SeqBench/External_Dependencies/fqzcomp/fqzcomp -Q 8 /work/tus53997/FASTQ/ERR016162_1.fastq /work/tus53997/CompressedOutput/ERR016162_1.fastq_-Q_8.fqz
END_TIME=$SECONDS
COMPRESSION_DURATION=$((END_TIME - START_TIME))

INPUT_SIZE=$(stat -c %s "/work/tus53997/FASTQ/ERR016162_1.fastq")
OUTPUT_SIZE=$(stat -c %s "/work/tus53997/CompressedOutput/ERR016162_1.fastq_-Q_8.fqz")
RATIO=$(echo "scale=6; $INPUT_SIZE / $OUTPUT_SIZE" | bc)

sleep 10

# Run decompression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/fqzcomp/fqzcomp -d /work/tus53997/CompressedOutput/ERR016162_1.fastq_-Q_8.fqz /work/tus53997/DecompressedOutput/ERR016162_1.fastq_-Q_8.fqz.fastq
END_TIME=$SECONDS
DECOMPRESSION_DURATION=$((END_TIME - START_TIME))

echo "job_0_1_0,fqzcomp_-Q_8,$COMPRESSION_DURATION,$DECOMPRESSION_DURATION,$RATIO" >> "/home/tus53997/SeqBench/Compression_Scripts/Logs/metrics/compression_metrics_ERR016162_1.fastq.csv"

conda deactivate