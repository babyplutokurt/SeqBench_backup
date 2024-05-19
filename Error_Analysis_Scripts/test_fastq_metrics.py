import time
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), 'build'))
import fastq_metrics

def main():
    original_file = "/home/tus53997/SeqBench/FASTQ/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq"
    decompressed_file = "/home/tus53997/SeqBench/DecompressedOutput/HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq_-Q_5.fqz.fastq"
    threads = 8


    # Test original version
    original_quality_scores = fastq_metrics.read_quality_scores(original_file)
    decompressed_quality_scores = fastq_metrics.read_quality_scores(decompressed_file)

    start_time = time.time()
    mse, psnr = fastq_metrics.calculate_mse_psnr(original_quality_scores, decompressed_quality_scores)
    end_time = time.time()
    
    print(f"Original Version - MSE: {mse}, PSNR: {psnr}, Time: {end_time - start_time} seconds")

    # Test v2 version
    start_time = time.time()
    mse_v2, psnr_v2 = fastq_metrics.calculate_mse_psnr_v2(original_file, decompressed_file)
    end_time = time.time()

    print(f"v2 Version - MSE: {mse_v2}, PSNR: {psnr_v2}, Time: {end_time - start_time} seconds")


    # Test v3 multithreaded version
    start_time = time.time()
    mse_v3, psnr_v3 = fastq_metrics.calculate_mse_psnr_v3(original_file, decompressed_file, threads)
    end_time = time.time()

    print(f"v3 Multithreaded Version - MSE: {mse_v3}, PSNR: {psnr_v3}, Time: {end_time - start_time} seconds")


if __name__ == "__main__":
    main()