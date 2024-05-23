import pandas as pd


def merge_csv_files(file1, file2, file3, output_file):
    """
    Merges three CSV files based on the 'Compressor_Name' column and saves the result to an output file.

    Parameters:
    file1 (str): Path to the first CSV file.
    file2 (str): Path to the second CSV file.
    file3 (str): Path to the third CSV file.
    output_file (str): Path to save the merged CSV file.
    """
    # Load the first CSV file
    df1 = pd.read_csv(file1)

    # Load the second CSV file
    df2 = pd.read_csv(file2)

    # Load the third CSV file
    df3 = pd.read_csv(file3)

    # Drop the 'job_id' column from the third dataframe
    df3 = df3.drop(columns=['job_id'])

    # Rename columns in the first dataframe to have a common 'job_id' for consistency if needed
    df1 = df1.rename(columns={'job_id_compression': 'job_id'})

    # Merge the first two dataframes on 'Compressor_Name'
    merged_df_1_2 = pd.merge(df1, df2, on='Compressor_Name')

    # Merge the result with the third dataframe on 'Compressor_Name'
    final_merged_df = pd.merge(merged_df_1_2, df3, on='Compressor_Name')

    # Save the merged dataframe to a new CSV file
    final_merged_df.to_csv(output_file, index=False)

    print(f'Final merged file saved as {output_file}')


if __name__ == "__main__":
    file1 = '/home/tus53997/SeqBench/Compression_Scripts/Logs/metrics/compression_metrics_HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq.csv'
    file2 = '/home/tus53997/SeqBench/Error_Analysis_Scripts/Logs/metrics/error_analysis_metrics_HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq.csv'
    file3 = '/home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/metrics/post_hoc_metrics_HG00097_CCAAGTCT-AAGGATGA_HCLHLDSXX_L004_001.R1_1000000.fastq.csv'
    output_file = '/home/tus53997/SeqBench/Post_Hoc_Scripts/Logs/metrics/final_merged_file.csv'

    merge_csv_files(file1, file2, file3, output_file)
