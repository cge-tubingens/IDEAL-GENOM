"""
Module to prepare data for the downstream analysis

The module provides a class to perform data preparation for downstream analysis.

Classes:
--------
PrepDS
    Class to perform data preparation for downstream analysis.
"""

import os
import psutil

from ideal_genom.Helpers import shell_do, delete_temp_files

class Preparatory:

    """
    Class designed to perform data preparation for downstream analysis.
    """

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, dependables:str) -> None:
        
        """
        Initializes the preparatory class with the given input and output paths, file names, and dependables path.
        
        Parameters:
        -----------
        input_path (str): 
            The path to the input directory.
        input_name (str): 
            The name of the input file (without extension).
        output_path (str): 
            The path to the output directory.
        output_name (str): 
            The name of the output file (without extension).
        dependables (str): 
            The path to the dependables directory.
        
        Raises:
        -------
        ValueError: 
            If any of input_path, output_path, or dependables is None.
        FileNotFoundError: 
            If input_path, output_path, or dependables does not exist.
        ValueError: 
            If input_name or output_name is None.
        TypeError: 
            If input_name or output_name is not a string.
        FileNotFoundError: 
            If any of the required PLINK files (.bed, .bim, .fam) are not found in the input_path.
        
        Attributes:
        -----------
        input_path (str): 
            The path to the input directory.
        output_path (str): 
            The path to the output directory.
        input_name (str): 
            The name of the input file (without extension).
        output_name (str): 
            The name of the output file (without extension).
        dependables (str): 
            The path to the dependables directory.
        results_dir (str): 
            The directory where results will be stored.
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

        # create results folder
        self.results_dir = os.path.join(output_path, 'preparatory')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        pass

    def execute_ld_prunning(self, maf:float=0.01, geno:float=0.1, hwe:float=5e-6, ind_pair:list=[50, 5, 0.2], memory:int=None)->dict:
        
        """
        Executes LD pruning using PLINK.
        
        This function performs linkage disequilibrium (LD) pruning on genotype data using PLINK. It filters SNPs based on minor allele frequency (MAF), genotype missingness, and Hardy-Weinberg equilibrium (HWE). It also excludes high LD regions and performs LD pruning using the specified parameters.

        Parameters:
        -----------
        maf : float, optional
            Minor allele frequency threshold (default is 0.01).
        geno : float, optional
            Genotype missingness threshold (default is 0.1).
        hwe : float, optional
            Hardy-Weinberg equilibrium threshold (default is 5e-6).
        ind_pair : list, optional
            Parameters for the --indep-pairwise option in PLINK (default is [50, 5, 0.2]).
        memory : int, optional
            Amount of memory (in MB) to allocate for PLINK (default is None, which calculates 2/3 of available memory).
        
        Returns:
        --------
        dict
            A dictionary containing the status of the process, the step name, and the output directory.
        
        Raises:
        -------
        TypeError
            If any of the parameters `maf`, `geno`, or `hwe` are not of type float.
        ValueError
            If any of the parameters `maf`, `geno`, or `hwe` are out of their respective valid ranges.
        FileNotFoundError
            If the file with high LD regions is not found.
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
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check if maf is in range
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 0.5")
        
        # Check if geno is in range
        if geno < 0 or geno > 1:
            raise ValueError("geno should be between 0 and 1")
        
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

        if memory is None:
            # get virtual memory details
            memory_info = psutil.virtual_memory()
            available_memory_mb = memory_info.available / (1024 * 1024)
            memory = round(2*available_memory_mb/3,0)

        # plink command to exclude high LD regions
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --chr 1-22 --maf {maf} --geno {geno}  --hwe {hwe} --exclude {high_ld_regions_file} --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --memory {memory} --make-bed --out {os.path.join(results_dir, output_name+'_prunning')}"

        # plink command to perform LD pruning
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, output_name+'_prunning')} --extract {os.path.join(results_dir, output_name+'_prunning.prune.in')} --make-bed --out {os.path.join(input_path, input_name+'-pruned')} --threads {max_threads}"

        # execute plink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

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
        Executes Principal Component Analysis (PCA) decomposition using PLINK.

        This method performs PCA decomposition on pruned genetic data using the PLINK software. It checks for the existence of necessary input files, constructs the PLINK command, and executes it. The results are stored in the specified output directory.

        Parameters:
        -----------
        pca (int, optional): 
            The number of principal components to compute. Default is 10.

        Returns:
        --------
        dict: 
            A dictionary containing the status of the process, the step name, and the output files.

        Raises:
        ------
        TypeError: 
            If `pca` is not of type int.
        ValueError: 
            If `pca` is less than 1.
        FileNotFoundError: 
            If the required pruned data files are not found.
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

        if not os.path.exists(os.path.join(input_path, input_name+'-pruned.bed')) or not os.path.exists(os.path.join(input_path, input_name+'-pruned.bim')) or not os.path.exists(os.path.join(input_path, input_name+'-pruned.fam')):
            raise FileNotFoundError(f"File with pruned data was not found: {os.path.join(input_path, input_name+'-pruned')}")

        # plink command to perform PCA decomposition
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name+'-pruned')} --pca {pca} --threads {max_threads} --out {os.path.join(input_path, input_name)}"

        # execute plink command
        shell_do(plink_cmd, log=True)

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
