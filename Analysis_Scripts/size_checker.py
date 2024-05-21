import os


def get_file_size(file_path):
    try:
        size = os.path.getsize(file_path)
        return size
    except OSError as e:
        print(f"Error: {e}")
        return None


# Example usage
if __name__ == "__main__":
    # Replace with the actual file path
    file_path = "/home/tus53997/SeqBench/DecompressedOutput/ERR1539081.fastq_-Q_5.fqz.fastq"
    size = get_file_size(file_path)
    if size is not None:
        print(f"The size of the file '{file_path}' is {size} bytes.")
    else:
        print("Could not retrieve the file size.")
