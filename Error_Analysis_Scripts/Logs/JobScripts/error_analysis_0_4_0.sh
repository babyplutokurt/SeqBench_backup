#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N error_analysis_0_4_0
#PBS -l nodes=1:ppn=24
#PBS -o /home/tus53997/SeqBench/Error_Analysis_Scripts/Logs/logs/job_0_4_0_output.log
#PBS -e /home/tus53997/SeqBench/Error_Analysis_Scripts/Logs/logs/job_0_4_0_error.log
#PBS -W depend=afterok:72217


# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run error analysis using the fastq_metrics API
mse_psnr_output=$(python3 -c "
import sys
sys.path.append('/home/tus53997/SeqBench/Error_Analysis_Scripts/build')
import fastq_metrics
original_file = '/home/tus53997/SeqBench/FASTQ/ERR103405_1.fastq'
decompressed_file = '/home/tus53997/SeqBench/DecompressedOutput/ERR103405_1.fastq_-Q_10.fqz.fastq'
threads = 24
mse_v3, psnr_v3 = fastq_metrics.calculate_mse_psnr_v3(original_file, decompressed_file, threads)
print(f'{mse_v3},{psnr_v3}')
")

# Extract MSE and PSNR values
mse_v3=$(echo $mse_psnr_output | cut -d',' -f1)
psnr_v3=$(echo $mse_psnr_output | cut -d',' -f2)

# Write the results to the CSV file
echo "error_analysis_0_4_0,fqzcomp_-Q_10,$mse_v3,$psnr_v3" >> "/home/tus53997/SeqBench/Error_Analysis_Scripts/Logs/metrics/error_analysis_metrics_ERR103405_1.fastq.csv"

conda deactivate