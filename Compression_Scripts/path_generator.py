import os
import json
from Analysis_Scripts.size_checker import get_file_size  # Ensure this path is correct based on your project structure


class PathGenerator:
    def __init__(self, config_path):
        self.config_path = config_path
        self.project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
        self.config = self.load_config()

    def load_config(self):
        try:
            with open(self.config_path) as json_file:
                return json.load(json_file)
        except FileNotFoundError:
            raise Exception(f"Configuration file not found: {self.config_path}")
        except json.JSONDecodeError:
            raise Exception(f"Error decoding JSON from the configuration file: {self.config_path}")

    def get_full_path(self, relative_path):
        if os.path.isabs(relative_path):
            return os.path.normpath(relative_path)
        # Normalize and resolve path to absolute
        full_path = os.path.abspath(os.path.join(self.project_base_dir, relative_path))
        return full_path

    def get_input_file_path(self, job_index, file_pair_index, file_index):
        job_name = self.config['jobs'][job_index]['name'].upper()
        file_set = "input_file_binary" if "SZ3" in job_name else "input_file"
        input_files = self.config[file_set][file_pair_index]
        return self.get_full_path(input_files[file_index])

    def get_compressed_output_path(self, job_index, file_pair_index, file_index):
        input_file_path = self.get_input_file_path(job_index, file_pair_index, file_index)
        job_name = self.config['jobs'][job_index]['name'].upper()
        option_str = self.config['jobs'][job_index]['options'][0]  # Get the first option

        if "{Binary_length}" in option_str:
            # Calculate binary length as required for SZ3
            file_size = get_file_size(input_file_path)
            binary_length = file_size // 4
            option_str = option_str.replace("{Binary_length}", str(binary_length))

        sanitized_option_str = option_str.replace(" ", "_").replace("/", "_").replace("{", "").replace("}", "").replace(
            "Binary_length", "len")

        suffix_mapper = {
            "SZ3": ".sz",
            "FQZCOMP": ".fqz",
            "SPRING": ".spring"
        }
        suffix = suffix_mapper.get(job_name, ".out")
        base_filename = os.path.basename(input_file_path)
        compressed_output_dir = os.path.abspath(os.path.join(self.project_base_dir, 'CompressedOutput'))
        compressed_file_name = f"{base_filename}_{sanitized_option_str}{suffix}"
        return os.path.join(compressed_output_dir, compressed_file_name)

    def get_decompressed_output_path(self, job_index, file_pair_index, file_index):
        compressed_path = self.get_compressed_output_path(job_index, file_pair_index, file_index)
        decompressed_output_dir = os.path.abspath(os.path.join(self.project_base_dir, 'DecompressedOutput'))
        decompressed_file_name = f"{os.path.basename(compressed_path)}.fastq"
        return os.path.join(decompressed_output_dir, decompressed_file_name)

    def get_compression_metric_path(self, file_pair_index, file_index):
        input_files = self.config['input_file'][file_pair_index]
        input_file_path = self.get_full_path(input_files[file_index])
        base_filename = os.path.basename(input_file_path)
        metrics_filename = f"compression_metrics_{base_filename}.csv"
        metrics_dir = os.path.abspath(os.path.join(self.project_base_dir, 'Compression_Scripts', 'Logs', 'metrics'))
        os.makedirs(metrics_dir, exist_ok=True)
        return os.path.join(metrics_dir, metrics_filename)

    def get_output_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Compression_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"job_{file_pair_index}_{job_index}_{file_index}_output.log")

    def get_error_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Compression_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"job_{file_pair_index}_{job_index}_{file_index}_error.log")

    def get_script_path(self, job_index, file_pair_index, file_index):
        scripts_dir = os.path.join(self.project_base_dir, 'Compression_Scripts/Logs/JobScripts')
        os.makedirs(scripts_dir, exist_ok=True)
        return os.path.join(scripts_dir, f"job_{file_pair_index}_{job_index}_{file_index}.sh")

    def get_compressor_name(self, job_index, file_pair_index, file_index):
        job = self.config['jobs'][job_index]
        compressor_name = job['name']
        options = job['options'][0]  # Get the first option for compression

        if "SZ3" in compressor_name:
            input_path = self.get_input_file_path(job_index, file_pair_index, file_index)
            file_size = get_file_size(input_path)
            binary_length = file_size // 4
            options = options.replace("{Binary_length}", str(binary_length))

        formatted_options = options.replace(" ", "_")
        return f"{compressor_name}_{formatted_options}"


if __name__ == "__main__":
    config_path = "/home/tus53997/SeqBench/Jobs/bench.json"  # Adjust the path as necessary
    pg = PathGenerator(config_path)
    job_index = 0
    file_pair_index = 1  # example file pair index
    file_index = 0  # example file index
    try:
        print(f"Compression Metric Path: {pg.get_compressor_name(job_index, file_pair_index, file_index)}")
    except Exception as e:
        print(f"Error: {e}")
