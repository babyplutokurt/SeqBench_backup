import os
import json
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
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
        full_path = os.path.abspath(os.path.join(self.project_base_dir, relative_path))
        return full_path

    def get_reference_file_path(self):
        reference_file = self.config['reference_file']
        return self.get_full_path(reference_file)

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

    def get_truth_sam_path(self, job_index, file_pair_index, file_index):
        input_file_path = self.get_input_file_path(job_index, file_pair_index, file_index)
        sam_filename = f"{os.path.basename(input_file_path)}.sam"
        sam_path = os.path.abspath(os.path.join(self.project_base_dir, 'SAM', sam_filename))
        return sam_path

    def get_truth_sorted_bam_path(self, job_index, file_pair_index, file_index):
        sam_path = self.get_truth_sam_path(job_index, file_pair_index, file_index)
        bam_filename = f"{os.path.basename(sam_path).replace('.sam', '')}_sorted.bam"
        bam_path = os.path.abspath(os.path.join(self.project_base_dir, 'BAM', bam_filename))
        return bam_path

    def get_truth_variant_path(self, job_index, file_pair_index, file_index):
        bam_path = self.get_truth_sorted_bam_path(job_index, file_pair_index, file_index)
        vcf_filename = f"{os.path.basename(bam_path).replace('_sorted.bam', '')}.vcf"
        vcf_path = os.path.abspath(os.path.join(self.project_base_dir, 'VCF', vcf_filename))
        return vcf_path

    def get_truth_compressed_variant_path(self, job_index, file_pair_index, file_index):
        vcf_path = self.get_truth_variant_path(job_index, file_pair_index, file_index)
        compressed_variant_path = vcf_path + '.gz'
        return compressed_variant_path

    def get_sam_path(self, job_index, file_pair_index, file_index):
        input_file_path = self.get_input_file_path(job_index, file_pair_index, file_index)
        job_options = self.config['jobs'][job_index]['options'][0]
        sanitized_options = job_options.replace(" ", "_").replace("/", "_")
        sam_filename = f"{os.path.basename(input_file_path)}_{sanitized_options}.sam"
        sam_path = os.path.abspath(os.path.join(self.project_base_dir, 'SAM', sam_filename))
        return sam_path

    def get_sorted_bam_path(self, job_index, file_pair_index, file_index):
        sam_path = self.get_sam_path(job_index, file_pair_index, file_index)
        bam_filename = f"{os.path.basename(sam_path).replace('.sam', '')}_sorted.bam"
        bam_path = os.path.abspath(os.path.join(self.project_base_dir, 'BAM', bam_filename))
        return bam_path

    def get_variant_path(self, job_index, file_pair_index, file_index):
        bam_path = self.get_sorted_bam_path(job_index, file_pair_index, file_index)
        vcf_filename = f"{os.path.basename(bam_path).replace('_sorted.bam', '')}.vcf"
        vcf_path = os.path.abspath(os.path.join(self.project_base_dir, 'VCF', vcf_filename))
        return vcf_path

    def get_compressed_variant_path(self, job_index, file_pair_index, file_index):
        vcf_path = self.get_variant_path(job_index, file_pair_index, file_index)
        compressed_variant_path = vcf_path + '.gz'
        return compressed_variant_path

    def get_truth_vcf_script_path(self, job_index, file_pair_index, file_index):
        scripts_dir = os.path.join(self.project_base_dir, 'Post_Hoc_Scripts/Logs/JobScripts')
        os.makedirs(scripts_dir, exist_ok=True)
        return os.path.join(scripts_dir, f"truth_vcf.sh")

    def get_lossy_vcf_script_path(self, job_index, file_pair_index, file_index):
        scripts_dir = os.path.join(self.project_base_dir, 'Post_Hoc_Scripts/Logs/JobScripts')
        os.makedirs(scripts_dir, exist_ok=True)
        return os.path.join(scripts_dir, f"post_hoc_{file_pair_index}_{job_index}_{file_index}.sh")

    def get_vcf_script_path(self, job_index, file_pair_index, file_index):
        scripts_dir = os.path.join(self.project_base_dir, 'Post_Hoc_Scripts/Logs/JobScripts')
        os.makedirs(scripts_dir, exist_ok=True)
        return os.path.join(scripts_dir, f"truth_vcf_{file_pair_index}_{job_index}.sh")

    def get_comparison_dir_path(self, job_index, file_pair_index, file_index):
        comparison_dir = os.path.join(self.project_base_dir, 'VCF', 'comparison',
                                      f"{file_pair_index}_{job_index}_{file_index}")
        os.makedirs(comparison_dir, exist_ok=True)
        return comparison_dir

    def get_post_hoc_metric_path(self, file_pair_index, file_index):
        input_files = self.config['input_file'][file_pair_index]
        input_file_path = self.get_full_path(input_files[0])
        base_filename = os.path.basename(input_file_path)
        metrics_filename = f"compression_metrics_{base_filename}.csv"
        metrics_dir = os.path.abspath(os.path.join(self.project_base_dir, 'Post_Hoc_Scripts', 'Logs', 'metrics'))
        os.makedirs(metrics_dir, exist_ok=True)
        return os.path.join(metrics_dir, metrics_filename)

    def get_post_hoc_output_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Post_Hoc_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"job_{file_pair_index}_{job_index}_{file_index}_output.log")

    def get_post_hoc_error_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Post_Hoc_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"job_{file_pair_index}_{job_index}_{file_index}_error.log")

    def get_post_hoc_truth_output_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Post_Hoc_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"truth_vcf_{file_pair_index}_output.log")

    def get_post_hoc_truth_error_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Post_Hoc_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"truth_vcf_{file_pair_index}_error.log")

    def get_error_analysis_metric_path(self, file_pair_index, file_index):
        input_files = self.config['input_file'][file_pair_index]
        input_file_path = self.get_full_path(input_files[file_index])
        base_filename = os.path.basename(input_file_path)
        metrics_filename = f"error_analysis_metrics_{base_filename}.csv"
        metrics_dir = os.path.abspath(os.path.join(self.project_base_dir, 'Error_Analysis_Scripts', 'Logs', 'metrics'))
        os.makedirs(metrics_dir, exist_ok=True)
        return os.path.join(metrics_dir, metrics_filename)

    def get_error_analysis_script_path(self, job_index, file_pair_index, file_index):
        scripts_dir = os.path.join(self.project_base_dir, 'Error_Analysis_Scripts', 'Logs', 'JobScripts')
        os.makedirs(scripts_dir, exist_ok=True)
        return os.path.join(scripts_dir, f"error_analysis_{file_pair_index}_{job_index}_{file_index}.sh")

    def get_build_cpp_path(self):
        build_path = os.path.join(self.project_base_dir, 'Error_Analysis_Scripts', 'build')
        return build_path

    def get_error_analysis_output_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Error_Analysis_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"job_{file_pair_index}_{job_index}_{file_index}_output.log")

    def get_error_analysis_error_log_path(self, job_index, file_pair_index, file_index):
        logs_dir = os.path.join(self.project_base_dir, 'Error_Analysis_Scripts/Logs/logs')
        os.makedirs(logs_dir, exist_ok=True)
        return os.path.join(logs_dir, f"job_{file_pair_index}_{job_index}_{file_index}_error.log")


if __name__ == "__main__":
    config_path = "/home/tus53997/SeqBench/Jobs/bench.json"  # Adjust the path as necessary
    pg = PathGenerator(config_path)
    job_index = 0
    file_pair_index = 0  # example file pair index
    file_index = 0  # example file index
    try:
        print({pg.get_input_file_path(job_index, file_pair_index, file_index)})
        # print({pg.get_sam_path(job_index, file_pair_index, file_index)})
        # print({pg.get_sorted_bam_path(job_index, file_pair_index, file_index)})
        # print({pg.get_variant_path(job_index, file_pair_index, file_index)})
        # print({pg.get_compressed_variant_path(job_index, file_pair_index, file_index)})
        # print({pg.get_truth_vcf_script_path(job_index, file_pair_index, file_index)})
        # print({pg.get_vcf_script_path(job_index, file_pair_index, file_index)})
    except Exception as e:
        print(f"Error: {e}")
