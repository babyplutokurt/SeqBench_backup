#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N job_0_1_0
#PBS -l nodes=1:ppn=24
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_1_0_output.log
#PBS -e /home/tus53997/SeqBench/Compression_Scripts/Logs/logs/job_0_1_0_error.log

# change to directory where 'qsub' was called
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run compression and capture metrics
START_TIME=$SECONDS

/home/tus53997/SeqBench/External_Dependencies/RENANO/renano/renano -t 24 /home/tus53997/SeqBench/FASTQ/ERR103405_1.fastq /home/tus53997/SeqBench/CompressedOutput/ERR103405_1.fastq_-t_24.renano
END_TIME=$SECONDS
COMPRESSION_DURATION=$((END_TIME - START_TIME))

INPUT_SIZE=$(stat -c %s "/home/tus53997/SeqBench/FASTQ/ERR103405_1.fastq")
OUTPUT_SIZE=$(stat -c %s "/home/tus53997/SeqBench/CompressedOutput/ERR103405_1.fastq_-t_24.renano")
RATIO=$(echo "scale=6; $INPUT_SIZE / $OUTPUT_SIZE" | bc)

sleep 10

# Run decompression and capture metrics
START_TIME=$SECONDS
/home/tus53997/SeqBench/External_Dependencies/RENANO/renano/renano -d -t 24 /home/tus53997/SeqBench/CompressedOutput/ERR103405_1.fastq_-t_24.renano /home/tus53997/SeqBench/DecompressedOutput/ERR103405_1.fastq_-t_24.renano.fastq
END_TIME=$SECONDS
DECOMPRESSION_DURATION=$((END_TIME - START_TIME))

echo "job_0_1_0,renano_-t_24,$COMPRESSION_DURATION,$DECOMPRESSION_DURATION,$RATIO" >> "/home/tus53997/SeqBench/Compression_Scripts/Logs/metrics/compression_metrics_ERR103405_1.fastq.csv"

conda deactivate