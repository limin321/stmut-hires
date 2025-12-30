import pandas as pd
import numpy as np
import os


class WtMedianWriter:
    """Handles writing weighted median of cnr"""
    def __init__(self, output_dir,df_arms, df_wmedian, cnr_name):
        self.output_dir = output_dir
        self.wtcnr_output_dir = None
        self.df_arms = df_arms
        self.df_wmedian = df_wmedian
        self.cnr_name = cnr_name
    
    def save_wt_cnr(self):
        self.wtcnr_output_dir = f"{self.output_dir}/wtcnr"
        os.makedirs(self.wtcnr_output_dir, exist_ok=True)
        print(f"save to {self.wtcnr_output_dir}")
        self.df_arms.to_csv(f"{self.wtcnr_output_dir}/summary.csv", index=False)
        self.df_wmedian.to_csv(f"{self.wtcnr_output_dir}/{self.cnr_name}", sep = "\t", index=False)

