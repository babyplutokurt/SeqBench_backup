#!/bin/sh
#PBS -l walltime={{ walltime }}
#PBS -N {{ job_name }}
#PBS -q normal
#PBS -l {{ nodes }}
#PBS -l mem={{ memory }}
#PBS -M {{ email }}
#PBS -o {{ output_log }}
#PBS -e {{ error_log }}

# change to directory where 'qsub' was called
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run compression and capture metrics
START_TIME=$SECONDS
{{ compression_command }}
END_TIME=$SECONDS
COMPRESSION_DURATION=$((END_TIME - START_TIME))

INPUT_SIZE=$(stat -c %s "{{ input_path }}")
OUTPUT_SIZE=$(stat -c %s "{{ compressed_output_path }}")
RATIO=$(echo "scale=2; $INPUT_SIZE / $OUTPUT_SIZE" | bc)

# Run decompression and capture metrics
START_TIME=$SECONDS
{{ decompression_command }}
END_TIME=$SECONDS
DECOMPRESSION_DURATION=$((END_TIME - START_TIME))

echo "{{ job_name }},{{ compressor_name }},$COMPRESSION_DURATION,$DECOMPRESSION_DURATION,$RATIO" >> "{{ metrics_csv_path }}"

conda deactivate