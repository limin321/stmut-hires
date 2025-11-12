import os

class OutputDirManager:
    """Responsible for managing output dirs"""
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.txt_output_dir = None

    def create_output_dir(self):
        """Create output_dir to store grouped barcodes expr"""
        self.txt_output_dir = f"{self.output_dir}/txt"
        os.makedirs(self.txt_output_dir, exist_ok=True)
        return self.txt_output_dir

