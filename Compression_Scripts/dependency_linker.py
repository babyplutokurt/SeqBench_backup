import os
import json

class DependencyLinker:
    def __init__(self, dependency_file):
        self.dependencies = None
        self.dependency_file = dependency_file
        self.load_dependencies()

    def load_dependencies(self):
        if os.path.exists(self.dependency_file):
            with open(self.dependency_file, 'r') as f:
                self.dependencies = json.load(f)
        else:
            self.dependencies = {}
            self.save_dependencies()

    def save_dependencies(self):
        with open(self.dependency_file, 'w') as f:
            json.dump(self.dependencies, f)

    def add_job_id(self, job_name, job_id):
        self.dependencies[job_name] = job_id
        self.save_dependencies()

    def get_dependency(self, job_name):
        return self.dependencies.get(job_name, None)
