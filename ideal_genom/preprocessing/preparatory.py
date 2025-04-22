"""
Module to prepare data for the GWAS and another advanced analysis.

The module provides a class to perform data preparation for advanced analysis.

Classes:
--------
Preparatory:
    Class to perform data preparation for downstream analysis.
"""

import os
import psutil
import logging

from ideal_genom.Helpers import shell_do
from ideal_genom.get_references import FetcherLDRegions

from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class Preparatory:

    def __init__(self, input_path: str | Path, input_name: str, output_path: str | Path, output_name: str, high_ld_file: str | Path, built: str = '38') -> None:
        
        """
        Initialize the preparatory class for genomic data processing.
        
        This class handles the initialization and validation of paths and files needed for genomic data preprocessing.
        It checks for the existence of required PLINK files and sets up the directory structure for output results.
        Parameters:
        ----------
        input_path (str | Path): 
            Directory path containing input PLINK files
        input_name (str): 
            Base name of input PLINK files (without extension)
        output_path (str | Path): 
            Directory path for output files
        output_name (str): 
            Base name for output files (without extension)
        high_ld_file (str | Path): 
            Path to the high LD regions file
        built (str, optional): 
            Genome build version. Must be either '38' or '37'. Defaults to '38'
        
        Raises:
        -------
            ValueError: If input_path, output_path, input_name, or output_name are not set,
                       or if built is not '38' or '37'
            TypeError: If input_path, output_path are not str or Path objects,
                       if input_name, output_name are not str objects,
                       or if built is not a str object
            FileNotFoundError: If input_path or output_path don't exist,
                              or if required PLINK files (.bed, .bim, .fam) are not found
        
        Note:
        -----
            If high_ld_file is not found, it will be automatically fetched from the package
            using the FetcherLDRegions class.
        """
        
        # check if paths are set
        if input_path is None or output_path:
            raise ValueError("Values for input_path and output_path must be set upon initialization.")
        
        if not isinstance(input_path, (str, Path)) or not isinstance(output_path, (str, Path)):
            raise TypeError("input_path and output_path should be of type str or Path.")
        
        input_path  = Path(input_path)
        output_path = Path(output_path)
        high_ld_file = Path(high_ld_file)
        
        if not input_path.exists() or not input_path.is_dir():
            raise FileNotFoundError(f"Input directory path does not exist: {input_path}")
        if not output_path.exists() or not output_path.is_dir():
            raise FileNotFoundError(f"Output path does not exist: {output_path}")
        
        # check if input_name and output_name are set
        if input_name is None or output_name is None:
            raise ValueError("Values for input_name and output_name must be set upon initialization.")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str.")
        
        if not isinstance(built, str):
            raise TypeError("built should be of type str.")
        if built not in ['38', '37']:
            raise ValueError("built should be either '38' or '37'.")
        
        # check existence of PLINK files
        if not (input_path / f"{input_name}.bed").exists():
            raise FileNotFoundError(f"PLINK bed file was not found: {input_path / f"{input_name}.bed"}")
        if not (input_path / f"{input_name}.bim").exists():
            raise FileNotFoundError(f"PLINK bim file was not found: {input_path / f"{input_name}.bim"}")
        if not (input_path / f"{input_name}.fam").exists():
            raise FileNotFoundError(f"PLINK fam file was not found: {input_path / f"{input_name}.fam"}")
        if not high_ld_file.is_file():
            logger.info(f"High LD file not found at {high_ld_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions(built=built)
            ld_fetcher.get_ld_regions()

            high_ld_file = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.built       = built
        self.high_ld_file = high_ld_file

        # create results folder
        self.results_dir = self.output_path / 'preparatory'
        self.results_dir.mkdir(parents=True, exist_ok=True)

        pass

    def execute_ld_prunning(self, mind: float = 0.2, maf: float = 0.01, geno: float = 0.1, hwe: float = 5e-6, ind_pair: list = [50, 5, 0.2], memory: int = None) -> None:
        """
        Execute LD (Linkage Disequilibrium) pruning on genetic data using PLINK.

        This method performs LD pruning in two steps:
        1. Excludes high LD regions and identifies independent SNPs
        2. Extracts the identified independent SNPs
        
        Parameters
        ----------
        mind : float, optional (default=0.2)
            Missing rate per individual threshold. Excludes individuals with missing rate higher than threshold.
        maf : float, optional (default=0.01) 
            Minor allele frequency threshold. Must be between 0 and 0.5.
        geno : float, optional (default=0.1)
            Missing rate per SNP threshold. Must be between 0 and 1.
        hwe : float, optional (default=5e-6)
            Hardy-Weinberg equilibrium exact test p-value threshold. Must be between 0 and 1.
        ind_pair : list, optional (default=[50, 5, 0.2])
            Parameters for pairwise pruning: [window size(variants), step size(variants), r^2 threshold]
        memory : int, optional (default=None)
            Memory in MB to allocate. If None, uses 2/3 of available system memory.
        
        Returns
        -------
        None
            The results are saved to disk and the pruned file path is stored in self.pruned_file
        
        Raises
        ------
        TypeError
            If mind, maf, geno, or hwe are not float
        ValueError
            If maf is not between 0 and 0.5
            If geno is not between 0 and 1
            If hwe is not between 0 and 1
        
        Notes
        -----
        Uses PLINK software for the pruning operations.
        Operates on chromosomes 1-22 only.
        Automatically determines optimal thread count based on system CPU cores.
        """
        
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")

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
        
        logger.info("STEP: LD pruning")

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
        plink_cmd1 = f"plink --bfile {self.input_path / self.input_name} --chr 1-22 --mind {mind} --maf {maf} --geno {geno}  --hwe {hwe} --exclude {self.high_ld_file} --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --keep-allele-order --threads {max_threads} --memory {memory} --make-bed --out {self.results_dir / (self.output_name+'-prunning')}"

        prune_in_file = (self.results_dir / (self.input_name+'-prunning')).with_suffix('.prune.in')

        # plink command to perform LD pruning
        plink_cmd2 = f"plink --bfile {self.results_dir / (self.output_name+'-prunning')} --extract {prune_in_file} --keep-allele-order --make-bed --out {self.input_path / (self.input_name+'-pruned')} --threads {max_threads}"

        # execute plink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.pruned_file = self.input_path / (self.input_name+'-pruned')

        return
    
    def execute_pc_decomposition(self, pca: int = 10) -> None:
        """
        Execute PCA decomposition on pruned PLINK binary files.

        This method performs Principal Component Analysis (PCA) on the pruned genotype data
        using PLINK software. It requires the existence of pruned binary PLINK files
        (.bed, .bim, .fam) and generates PCA eigenvectors and eigenvalues.

        Parameters
        ----------
        pca : int, default=10
            Number of principal components to compute. Must be greater than 0.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If pca parameter is not an integer.
        ValueError
            If pca parameter is less than 1.
        FileNotFoundError
            If any of the required pruned PLINK files (.bed, .bim, .fam) are not found.

        Notes
        -----
        The method automatically determines the optimal number of threads to use based on
        CPU count, reserving 2 cores for other processes. If CPU count cannot be determined,
        it defaults to 10 threads.

        The output files will be created in the same directory as the input files, using
        the input name as prefix with extensions .eigenvec and .eigenval.
        """

        # Check type of pca and range
        if not isinstance(pca, int):
            raise TypeError("pca should be of type int.")
        if pca < 1:
            raise ValueError("pca should be greater than 0.")

        logger.info("STEP: PCA decomposition")

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        if not self.pruned_file.with_suffix('.bed').exists():
            raise FileNotFoundError(f"bed file with pruned data was not found: {self.pruned_file.with_suffix('.bed')}")
        if not self.pruned_file.with_suffix('.bim').exists():
            raise FileNotFoundError(f"bim file with pruned data was not found: {self.pruned_file.with_suffix('.bim')}")
        if not self.pruned_file.with_suffix('.fam').exists():
            raise FileNotFoundError(f"fam file with pruned data was not found: {self.pruned_file.with_suffix('.fam')}")

        # plink command to perform PCA decomposition
        plink_cmd = f"plink --bfile {self.input_path / (self.input_name+'-pruned')} --pca {pca} --threads {max_threads} --out {self.input_path / self.input_name}"

        # execute plink command
        shell_do(plink_cmd, log=True)

        return
