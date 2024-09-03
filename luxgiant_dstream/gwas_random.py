"""
Module to perform a GWAS analysis using a random effect model.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from luxgiant_dstream.Helpers import shell_do, delete_temp_files

class GWASrandom:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, preps_path:str) -> None:

        # check if paths are set
        if input_path is None or output_path is None:
            raise ValueError("Values for input_path, output_path and dependables_path must be set upon initialization.")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.config_dict = config_dict
        self.preps_path  = preps_path

        self.files_to_keep = []

        # create results folder
        self.results_dir = os.path.join(output_path, 'gwas_random')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create figures folder
        self.plots_dir = os.path.join(output_path, 'plots_random')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

        pass

    def prepare_aux_files(self)->dict:
        
        input_path  = self.input_path
        input_name  = self.input_name
        results_dir  = self.results_dir
        output_name = self.output_name

        step = "prepare_aux_files"

        df_fam = pd.read_csv(os.path.join(input_path, input_name+'.fam'), sep='\t', header=None)

        df_pheno = df_fam[[0,1,5]].copy()

        # recode phenotype
        df_pheno[5] = df_pheno[5]-1

        df_pheno.to_csv(os.path.join(results_dir, output_name+'_pheno.phen'), sep='\t', header=False, index=False)

        df_sex = df_fam[[0,1,4]].copy()

        df_sex.to_csv(os.path.join(results_dir, output_name+'_sex.covar'), sep='\t', header=False, index=False)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def compute_grm(self)->dict:

        return dict()