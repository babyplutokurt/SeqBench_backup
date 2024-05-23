#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N job_0_7_0
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_7_0_output.log
#PBS -e /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_7_0_error.log

# change to directory where 'qsub' was called
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run compression and capture metrics
START_TIME=$SECONDS

/home/tus53997/SeqBench/External_Dependencies/SPRING/build/spring -c -l -q qvz 0.6 -i /work/tus53997/FASTQ/ERR016162_1.fastq -o /work/tus53997/CompressedOutput/ERR016162_1.fastq_-c_-l_-q_qvz_0.6.spring
END_TIME=$SECONDS
COMPRESSION_DURATION=$((END_TIME - START_TIME))

INPUT_SIZE=$(stat -c %s "/work/tus53997/FASTQ/ERR016162_1.fastq")
OUTPUT_SIZE=$(stat -c %s "/work/tus53997/CompressedOutput/ERR016162_1.fastq_-c_-l_-q_qvz_0.6.spring")
RATIO=$(echo "scale=6; $INPUT_SIZE / $OUTPUT_SIZE" | bc)

sleep 10

# Run decompression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/SPRING/build/spring -d -i /work/tus53997/CompressedOutput/ERR016162_1.fastq_-c_-l_-q_qvz_0.6.spring -o /work/tus53997/DecompressedOutput/ERR016162_1.fastq_-c_-l_-q_qvz_0.6.spring.fastq
END_TIME=$SECONDS
DECOMPRESSION_DURATION=$((END_TIME - START_TIME))

echo "job_0_7_0,spring_-c_-l_-q_qvz_0.6,$COMPRESSION_DURATION,$DECOMPRESSION_DURATION,$RATIO" >> "/home/tus53997/SeqBench/Compression_Scripts/Logs/metrics/compression_metrics_ERR016162_1.fastq.csv"

conda deactivate