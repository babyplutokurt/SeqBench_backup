import os
import csv
import subprocess
import logging
from jinja2 import Template
from command_generator import CommandGeneratorFactory
from dependency_linker import DependencyLinker


class JobGenerator:
    def __init__(self, config_name, template_path="/home/tus53997/SeqBench/Scripts_Template/Nest.sh",
                 dependency_file="/home/tus53997/SeqBench/Compression_Scripts/Logs/job_dependencies.json"):
        self.factory = CommandGeneratorFactory(config_name)
        self.config = self.factory.config
        self.path_generator = self.factory.path_generator
        self.template_path = template_path
        self.dependency_linker = DependencyLinker(dependency_file)
        self.all_commands = self.factory.generate_all_commands()

    def create_compression_metrics_csv(self, file_pair_index, file_index):
        metrics_path = self.path_generator.get_compression_metric_path(file_pair_index, file_index)
        header = ['job_id_compression', 'Compressor_Name', 'Compression_Time', 'Decompression_Time', 'Ratio']
        os.makedirs(os.path.dirname(metrics_path), exist_ok=True)
        with open(metrics_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
        logging.info(f"Created CSV file for compression metrics: {metrics_path}")

    def create_job_script(self, job_name, commands, job_header, job_index, file_pair_index, file_index,
                          dependencies=None):
        job_script_path = self.path_generator.get_script_path(job_index, file_pair_index, file_index)

        with open(self.template_path) as f:
            template = Template(f.read())

        compression_command = commands[0]
        decompression_command = commands[1]
        input_path = self.path_generator.get_input_file_path(job_index, file_pair_index, file_index)
        compressed_output_path = self.path_generator.get_compressed_output_path(job_index, file_pair_index, file_index)
        metrics_csv_path = self.path_generator.get_compression_metric_path(file_pair_index, file_index)
        output_log = self.path_generator.get_output_log_path(job_index, file_pair_index, file_index)
        error_log = self.path_generator.get_error_log_path(job_index, file_pair_index, file_index)
        compressor_name = self.path_generator.get_compressor_name(job_index, file_pair_index, file_index)
        nodes = self.config.get('nodes', 1)  # Default to 1 node if not specified
        ppn = self.config.get('ppn', 8)  # Default to 8 processor per node if not specified
        conda_path = self.config.get('conda_path', '/home/tus53997/miniconda3/bin/activate')
        walltime = self.config.get('walltime', "24:00:00")
        email = self.config.get('email', "default@gamil.com")

        dependency_line = f"#PBS -W depend=afterok:{':'.join(dependencies)}\n" if dependencies else ""

        job_script_content = template.render(
            job_name=job_name,
            nodes=nodes,
            ppn=ppn,
            walltime=walltime,
            conda_path=conda_path,
            email=email,
            output_log=output_log,
            error_log=error_log,
            compression_command=compression_command,
            decompression_command=decompression_command,
            input_path=input_path,
            compressed_output_path=compressed_output_path,
            metrics_csv_path=metrics_csv_path,
            compressor_name=compressor_name,
            dependency_line=dependency_line
        )

        with open(job_script_path, 'w') as f:
            f.write(job_script_content)

        logging.info(f"Created job script: {job_script_path}")
        return job_script_path

    def submit_job(self, job_script_path, job_name, use_append=False):
        submit_command = f"qsub {job_script_path}"  # Change this if using a different scheduler like `sbatch` for SLURM
        try:
            result = subprocess.run(submit_command, check=True, shell=True, capture_output=True, text=True)
            job_id = result.stdout.strip().split('.')[0]  # Get the job ID from the output
            if use_append:
                self.dependency_linker.append_job_id(job_name, job_id)
            else:
                self.dependency_linker.add_job_id(job_name, job_id)
            logging.info(f"Job submitted: {job_script_path} with job ID: {job_id}")
            return job_id
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to submit job: {job_script_path}, Error: {e}")
            return None

    def generate_and_submit_jobs(self):
        # Step 1: Create CSV files for each input file
        for file_pair_index, file_pair in enumerate(self.config['input_file']):
            for file_index in range(len(file_pair)):
                self.create_compression_metrics_csv(file_pair_index, file_index)

        # Step 2: Create and submit job scripts with dependencies
        previous_job_ids = []
        for file_pair_index, file_pair in enumerate(self.config['input_file']):
            for job_index in range(len(self.config['jobs'])):
                for file_index in range(len(file_pair)):
                    job_name = f"job_{file_pair_index}_{job_index}_{file_index}"
                    # job_script_path = self.path_generator.get_script_path(job_index, file_pair_index, file_index)
                    job_script_path = self.create_job_script(
                        job_name,
                        self.all_commands[file_pair_index][job_index][file_index],
                        self.config['jobs_header'],
                        job_index,
                        file_pair_index,
                        file_index,
                        dependencies=previous_job_ids
                    )
                    job_id = self.submit_job(job_script_path, job_name)
                    if job_id:
                        self.dependency_linker.add_job_id(job_name, job_id)
                        previous_job_ids.append(job_id)  # Only append if the job was successfully submitted


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname=s - %(message=s')
    config_name = "/home/tus53997/SeqBench/Jobs/bench.json"
    job_generator = JobGenerator(config_name)
    job_generator.generate_and_submit_jobs()
