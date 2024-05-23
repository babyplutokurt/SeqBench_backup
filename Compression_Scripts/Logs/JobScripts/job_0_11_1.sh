#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N job_0_11_1
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_11_1_output.log
#PBS -e /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_11_1_error.log

# change to directory where 'qsub' was called
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run compression and capture metrics
START_TIME=$SECONDS

/home/tus53997/SeqBench/External_Dependencies/ENANO/enano/enano -t 28 /work/tus53997/FASTQ/ERR016162_2.fastq /work/tus53997/CompressedOutput/ERR016162_2.fastq_-t_28.enano
END_TIME=$SECONDS
COMPRESSION_DURATION=$((END_TIME - START_TIME))

INPUT_SIZE=$(stat -c %s "/work/tus53997/FASTQ/ERR016162_2.fastq")
OUTPUT_SIZE=$(stat -c %s "/work/tus53997/CompressedOutput/ERR016162_2.fastq_-t_28.enano")
RATIO=$(echo "scale=6; $INPUT_SIZE / $OUTPUT_SIZE" | bc)

sleep 10

# Run decompression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/ENANO/enano/enano -d -t 28 /work/tus53997/CompressedOutput/ERR016162_2.fastq_-t_28.enano /work/tus53997/DecompressedOutput/ERR016162_2.fastq_-t_28.enano.fastq
END_TIME=$SECONDS
DECOMPRESSION_DURATION=$((END_TIME - START_TIME))

echo "job_0_11_1,enano_-t_28,$COMPRESSION_DURATION,$DECOMPRESSION_DURATION,$RATIO" >> "/home/tus53997/SeqBench/Compression_Scripts/Logs/metrics/compression_metrics_ERR016162_2.fastq.csv"

conda deactivate