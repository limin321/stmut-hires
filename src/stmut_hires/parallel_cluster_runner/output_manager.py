import os
from pathlib import Path
import shutil

class OutputDirManager:
    """Responsible for managing output dirs"""
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.txt_output_dir = None
        self.cnr_dir = None
        self.wtcnr_dir = None
        self.cdt_dir = None
        self.figures_dir = None
        self.cluster_summary = None

    def create_output_dir(self):
        """Create output_dir to store grouped barcodes expr"""
        self.txt_output_dir = f"{self.output_dir}/txt"
        os.makedirs(self.txt_output_dir, exist_ok=True)

        self.cnr_dir = f"{self.output_dir}/cnr"
        os.makedirs(self.cnr_dir, exist_ok=True)

        self.wtcnr_dir = f"{self.output_dir}/wtcnr"
        os.makedirs(self.wtcnr_dir, exist_ok=True)

        self.cdt_dir = f"{self.output_dir}/cdt"
        os.makedirs(self.cdt_dir, exist_ok=True)

        self.figures_dir = f"{self.output_dir}/figures"
        os.makedirs(self.figures_dir, exist_ok=True)

        self.cluster_summary = f"{self.output_dir}/cluster_summary"
        os.makedirs(self.cluster_summary, exist_ok=True)
        return self.txt_output_dir, self.cnr_dir, self.wtcnr_dir, self.cdt_dir,self.figures_dir,self.cluster_summary

    def move_csvs_to_summary(self):
        """Moves all CSV files from txt directory to cluster_summary folder"""
        # Convert to Path objects for easier globbing
        txt_path = Path(self.txt_output_dir)

        # Use .glob() to find matching files
        for csv_file in txt_path.glob("*.csv"):
            shutil.move(str(csv_file), self.cluster_summary)

        


