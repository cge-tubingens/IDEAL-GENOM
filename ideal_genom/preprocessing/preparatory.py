"""
Module to prepare data for the downstream analysis

The module provides a class to perform data preparation for downstream analysis.

Classes:
--------
PrepDS
    Class to perform data preparation for downstream analysis.
"""

import os

from ideal_genom.Helpers import shell_do, delete_temp_files

class Preparatory:

    """
    Class designed to perform data preparation for downstream analysis.
    """

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, dependables:str) -> None:

        """
        Initialize the PrepDS class.

        Parameters:
        -----------
        input_path : str
            Path to the input data.
        input_name : str
            Name of the input data.
        output_path : str
            Path to the output data.
        output_name : str
            Name of the output data.
        config_dict : dict
            Dictionary containing the configuration parameters.
        dependables_path : str
            Path to the dependables data.

        Returns:
        --------
        None
        """
        
        # check if paths are set
        if input_path is None or output_path is None or dependables is None:
            raise ValueError("Values for input_path, output_path and dependables_path must be set upon initialization.")
        
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input path does not exist: {input_path}")
        if not os.path.exists(dependables):
            raise FileNotFoundError(f"Dependables path does not exist: {dependables}")
        if not os.path.exists(output_path):
            raise FileNotFoundError(f"Output path does not exist: {output_path}")
        
        # check if input_name and output_name are set
        if input_name is None or output_name is None:
            raise ValueError("Values for input_name and output_name must be set upon initialization.")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str.")
        
        # check existence of PLINK files
        if not os.path.exists(os.path.join(input_path, input_name+'.bed')):
            raise FileNotFoundError(f"PLINK bed file was not found: {os.path.join(input_path, input_name+'.bed')}")
        if not os.path.exists(os.path.join(input_path, input_name+'.bim')):
            raise FileNotFoundError(f"PLINK bim file was not found: {os.path.join(input_path, input_name+'.bim')}")
        if not os.path.exists(os.path.join(input_path, input_name+'.fam')):
            raise FileNotFoundError(f"PLINK fam file was not found: {os.path.join(input_path, input_name+'.fam')}")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.dependables = dependables
        
        self.files_to_keep = []

        # create results folder
        self.results_dir = os.path.join(output_path, 'preparatory')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        pass

    def execute_ld_prunning(self, maf:float=0.05, geno:float=0.1, mind:float=0.1, hwe:float=5e-8, ind_pair:list=[50, 5, 0.2])->dict:

        """
        Method to exclude high LD regions and perform LD pruning.
        
        Returns:
        --------
        out_dict : dict
            Dictionary containing a report of the process.
        """

        input_path       = self.input_path
        input_name       = self.input_name
        results_dir      = self.results_dir
        output_name      = self.output_name
        dependables_path = self.dependables

        # Check type of maf
        if not isinstance(maf, float):
             raise TypeError("maf should be of type float.")

        # Check type of geno
        if not isinstance(geno, float):
            raise TypeError("geno should be of type float.")

        # Check type of mind
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check if maf is in range
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 0.5")
        
        # Check if geno is in range
        if geno < 0 or geno > 1:
            raise ValueError("geno should be between 0 and 1")
        
        # Check if mind is in range
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1")
        
        # Check if hwe is in range
        if hwe < 0 or hwe > 1:
            raise ValueError("hwe should be between 0 and 1")
        
        # check existence of high LD regions file
        high_ld_regions_file = os.path.join(dependables_path, 'high-LD-regions.txt')
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError(f"File with high LD region was not found: {high_ld_regions_file}")

        step = "ld_prune"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # plink command to exclude high LD regions
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --chr 1-22 --maf {maf} --geno {geno}  --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --make-bed --out {os.path.join(results_dir, output_name+'_prunning')}"

        # plink command to perform LD pruning
        plink_cmd2 = f"plink2 --bfile {os.path.join(results_dir, output_name+'_prunning')} --extract {os.path.join(results_dir, output_name+'_prunning.prune.in')} --make-bed --out {os.path.join(results_dir, output_name)} --threads {max_threads}"

        # execute plink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # self.files_to_keep.append(output_name.bed)
        # self.files_to_keep.append(output_name.bim)
        # self.files_to_keep.append(output_name.fam)

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
    
    def execute_pc_decomposition(self, pca:int=10)->dict:

        """
        Method to perform PCA decomposition.

        Returns:
        --------
        out_dict : dict
            Dictionary containing a report of the process
        """

        results_dir = self.results_dir
        input_name = self.input_name
        input_path = self.input_path
        output_name = self.output_name

        # Check type of pca and range
        if not isinstance(pca, int):
            raise TypeError("pca should be of type int.")
        if pca < 1:
            raise ValueError("pca should be greater than 0.")

        step = "pca_decomposition"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        if not os.path.exists(os.path.join(results_dir, output_name+'.bed')) or not os.path.exists(os.path.join(results_dir, output_name+'.bim')) or not os.path.exists(os.path.join(results_dir, output_name+'.fam')):
            raise FileNotFoundError(f"File with pruned data was not found: {os.path.join(results_dir, output_name)}")

        # plink command to perform PCA decomposition
        plink_cmd = f"plink --bfile {os.path.join(results_dir, output_name)} --pca {pca} --threads {max_threads} --out {os.path.join(input_path, input_name)}"

        # execute plink command
        shell_do(plink_cmd, log=True)

        self.files_to_keep.append(output_name+'.eigenvec')

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
