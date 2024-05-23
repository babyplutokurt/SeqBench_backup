#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N error_analysis_0_12_0
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -M taolue.yang@temple.edu
#PBS -o /home/tus53997/SeqBench/Error_Analysis_Scripts/Logs/logs/job_0_12_0_output.log
#PBS -e /home/tus53997/SeqBench/Error_Analysis_Scripts/Logs/logs/job_0_12_0_error.log

#PBS -W depend=afterok:1893698


# Change to the working directory
cd $PBS_O_WORKDIR

source /home/tus53997/miniconda3/bin/activate compression

# Run error analysis using the fastq_metrics API
mse_psnr_output=$(python3 -c "
import sys
sys.path.append('/home/tus53997/SeqBench/Error_Analysis_Scripts/build')
import fastq_metrics
original_file = '/work/tus53997/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq'
decompressed_file = '/work/tus53997/DecompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq_-t_28.enano.fastq'
threads = 28
mse_v3, psnr_v3 = fastq_metrics.calculate_mse_psnr_v3(original_file, decompressed_file, threads)
print(f'{mse_v3},{psnr_v3}')
")

# Extract MSE and PSNR values
mse_v3=$(echo $mse_psnr_output | cut -d',' -f1)
psnr_v3=$(echo $mse_psnr_output | cut -d',' -f2)

# Write the results to the CSV file
echo "error_analysis_0_12_0,enano_-t_28,$mse_v3,$psnr_v3" >> "/home/tus53997/SeqBench/Error_Analysis_Scripts/Logs/metrics/error_analysis_metrics_HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq.csv"

conda deactivate