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
        results_dir = self.results_dir
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
            'python_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def compute_grm(self)->dict:

        results_dir= self.results_dir
        prep_path  = self.preps_path
        output_name= self.output_name

        step = "compute_grm"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # gcta commands
        gcta_cmd1 = f"gcta64 --bfile {os.path.join(prep_path, output_name+'_LDpruned')} --make-grm --thread-num {max_threads} --out {os.path.join(results_dir, output_name+'_grm')}"

        gcta_cmd2 = f"gcta64 --grm {os.path.join(results_dir, output_name+'_grm')} --make-bK-sparse 0.05 --out {os.path.join(results_dir, output_name+'_sparse')}"

        # run gctag4 commands
        cmds = [gcta_cmd1, gcta_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'gcta_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def run_gwas_random(self)->dict:

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name
        config_dict = self.config_dict
        preps_path  = self.preps_path

        maf = config_dict['maf']

        step = "run_gwas_random"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # gcta commands
        gcta_cmd = f"gcta64 --bfile {os.path.join(input_path, input_name)} --fastGWA-mlm-binary --maf {maf}  --grm-sparse {os.path.join(results_dir, output_name+'_sparse')} --qcovar {os.path.join(preps_path, output_name+'_pca.eigenvec')} --covar {os.path.join(results_dir, output_name+'_sex.covar')} --pheno {os.path.join(results_dir, output_name+'_pheno.phen')} --out {os.path.join(results_dir,output_name+'_assocSparseCovar_pca_sex-mlm-binary')}--thread-num {max_threads}"

        # run gcta command
        shell_do(gcta_cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'gcta_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
