import os
import json

import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Analysis_Scripts.size_checker import get_file_size
from path_generator import PathGenerator
from compressor_paths import COMPRESSOR_PATHS
import logging


class CommandGenerator:
    def __init__(self, config_path, path_generator):
        self.config_path = config_path
        self.path_generator = path_generator
        self.config = self.path_generator.load_config()

    def generate_commands(self, job_index, file_pair_index, file_index):
        raise NotImplementedError("Subclasses should implement this method.")


class SZ3CommandGenerator(CommandGenerator):
    def generate_commands(self, job_index, file_pair_index, file_index):
        job = self.config['jobs'][job_index]
        if job['name'].upper() != "SZ3":
            return []

        executable_path = COMPRESSOR_PATHS.get("SZ3")
        if not executable_path:
            raise ValueError("Path for SZ3 compressor not found.")

        input_path = self.path_generator.get_input_file_path(job_index, file_pair_index, file_index)
        output_path = self.path_generator.get_compressed_output_path(job_index, file_pair_index, file_index)
        decompressed_path = self.path_generator.get_decompressed_output_path(job_index, file_pair_index, file_index)

        file_size = get_file_size(input_path)
        binary_length = file_size // 4  # Calculate binary length as required

        # Replace {Binary_length} in the command options
        compression_options = job['options'][0].replace("{Binary_length}", str(binary_length))
        decompression_options = job['options'][1].replace("{Binary_length}", str(binary_length))

        compression_command = f"{executable_path} {compression_options} -i {input_path} -z {output_path}"
        decompression_command = f"{executable_path} {decompression_options} -z {output_path} -o {decompressed_path}"

        return [compression_command, decompression_command]


class FQZCompCommandGenerator(CommandGenerator):
    def generate_commands(self, job_index, file_pair_index, file_index):
        job = self.config['jobs'][job_index]
        if job['name'].upper() != "FQZCOMP":
            return []

        executable_path = COMPRESSOR_PATHS.get("FQZCOMP")
        if not executable_path:
            raise ValueError("Path for FQZComp compressor not found.")

        input_path = self.path_generator.get_input_file_path(job_index, file_pair_index, file_index)
        output_path = self.path_generator.get_compressed_output_path(job_index, file_pair_index, file_index)
        decompressed_path = self.path_generator.get_decompressed_output_path(job_index, file_pair_index, file_index)

        commands = []
        for option in job['options']:
            if "-d" in option:  # Assuming "-d" indicates a decompression option
                # Modify the input/output paths for decompression
                decompression_input_path = output_path  # Decompress the previously compressed file
                command = f"{executable_path} {option} {decompression_input_path} {decompressed_path}"
            else:
                command = f"{executable_path} {option} {input_path} {output_path}"
            commands.append(command)

        return commands


class SpringCommandGenerator(CommandGenerator):
    def generate_commands(self, job_index, file_pair_index, file_index):
        job = self.config['jobs'][job_index]
        normalized_job_name = job['name'].upper()
        if normalized_job_name != "SPRING":
            return []

        executable_path = COMPRESSOR_PATHS.get("SPRING")
        if not executable_path:
            raise ValueError("Path for Spring compressor not found.")

        input_path = self.path_generator.get_input_file_path(job_index, file_pair_index, file_index)
        output_path = self.path_generator.get_compressed_output_path(job_index, file_pair_index, file_index)
        decompressed_path = self.path_generator.get_decompressed_output_path(job_index, file_pair_index, file_index)

        commands = []
        for option in job['options']:
            if "-d" in option:  # Assuming "-d" indicates a decompression option
                command = f"{executable_path} {option} -i {output_path} -o {decompressed_path}"
            else:
                command = f"{executable_path} {option} -i {input_path} -o {output_path}"
            commands.append(command)

        return commands


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class CommandGeneratorFactory:
    def __init__(self, config_name):
        try:
            self.config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'jobs', config_name)
            self.config = self.load_config()
            self.path_generator = PathGenerator(self.config_path)
            self.generators = {
                'SZ3': SZ3CommandGenerator(self.config_path, self.path_generator),
                'FQZCOMP': FQZCompCommandGenerator(self.config_path, self.path_generator),
                'SPRING': SpringCommandGenerator(self.config_path, self.path_generator)
            }
        except Exception as e:
            logging.error(f"Failed to initialize command generators: {e}")
            raise

    def load_config(self):
        try:
            with open(self.config_path) as json_file:
                return json.load(json_file)
        except FileNotFoundError:
            logging.error(f"Configuration file not found: {self.config_path}")
            raise
        except json.JSONDecodeError:
            logging.error("JSON decoding error occurred while reading the configuration.")
            raise

    def generate_all_commands(self):
        all_commands_for_files = []
        for file_pair_index, file_pair in enumerate(self.config['input_file']):
            commands_for_current_file_pair = []
            for job_index, job in enumerate(self.config['jobs']):
                command_for_current_job = []
                generator = self.generators.get(job['name'].upper())
                if generator:
                    for file_index in range(len(file_pair)):
                        try:
                            job_commands = generator.generate_commands(job_index, file_pair_index, file_index)
                            command_for_current_job.append(job_commands)
                        except Exception as e:
                            logging.warning(
                                f"Failed to generate commands for job {job['name']} at index {job_index}: {e}")
                else:
                    logging.info(f"No generator found for {job['name']}. Skipping...")
                commands_for_current_file_pair.append(command_for_current_job)
            all_commands_for_files.append(commands_for_current_file_pair)
        return all_commands_for_files


if __name__ == "__main__":

    config_path = "/home/tus53997/SeqBench/Jobs/bench.json"
    factory = CommandGeneratorFactory(config_path)
    all_commands = factory.generate_all_commands()
    for file_pair_index, file_commands in enumerate(all_commands):
        for job_index, command_set in enumerate(file_commands):
            for file_index, command in enumerate(command_set):
                print(f"File Pair: {file_pair_index}, File: {file_index}, Job: {job_index}")
                for cmd in command:
                    print(cmd)
                print('-' * 50)
