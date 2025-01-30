"""
Module to perform a GWAS analysis using a fixed model.
"""

import gzip
import os
import shutil

import pandas as pd

from ideal_genom.Helpers import shell_do

from gwaslab.bd_download import download_file
from gwaslab.util_in_get_sig import annogene
from gwaslab.g_Log import Log

class GWASfixed:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, recompute:bool=True) -> None:
        """
        Initialize the GWAS fixed model analysis.

        Parameters:
        -----------
        input_path (str): 
            Path to the input directory containing PLINK files.
        input_name (str): 
            Base name of the input PLINK files (without extensions).
        output_path (str): 
            Path to the output directory where results will be saved.
        output_name (str): 
            Base name for the output files.
        recompute (bool): 
            Flag indicating whether to recompute the analysis if results already exist. Default is True.
        
        Raises:
        -------
        ValueError: 
            If input_path, output_path, input_name, or output_name are not provided.
        FileNotFoundError: 
            If the specified input_path or output_path does not exist.
        FileNotFoundError: 
            If the required PLINK files (.bed, .bim, .fam) are not found in the input_path.
        TypeError: 
            If input_name or output_name are not strings, or if recompute is not a boolean.
        
        Attributes:
        ------------
        input_path (str): 
            Path to the input directory.
        output_path (str): 
            Path to the output directory.
        input_name (str): 
            Base name of the input PLINK files.
        output_name (str): 
            Base name for the output files.
        recompute (bool): 
            Flag indicating whether to recompute the analysis.
        results_dir (str): 
            Directory where the results will be saved.
        """

    
        # check if paths are set
        if input_path is None or output_path is None:
            raise ValueError("Values for input_path, output_path and dependables_path must be set upon initialization.")
        
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input path does not exist: {input_path}")
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
        
        if not isinstance(recompute, bool):
            raise TypeError("recompute should be of type bool.")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.recompute   = recompute

        # create results folder
        self.results_dir = os.path.join(output_path, 'gwas_fixed')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        print("\033[1;32mAnalysis of GWAS data using a fixed model initialized.\033[0m")

        pass

    def fixed_model_association_analysis(self, maf:float, mind:float, hwe:float, ci:float)->dict:
        """
        Perform fixed model association analysis using PLINK2.
        
        This method performs a fixed model association analysis on genomic data using PLINK2. It checks the validity of the input parameters, ensures necessary files exist, and executes the PLINK2 command to perform the analysis.

        Parameters:
        ------------
        maf (float): 
            Minor allele frequency threshold. Must be between 0 and 0.5.
        mind (float): 
            Individual missingness threshold. Must be between 0 and 1.
        hwe (float): 
            Hardy-Weinberg equilibrium threshold. Must be between 0 and 1.
        ci (float): 
            Confidence interval threshold. Must be between 0 and 1.
        
        Returns:
        -------
        dict: 
            A dictionary containing the status of the process, the step name, and the output directory.
        
        Raises:
        -------
        TypeError: 
            If any of the input parameters are not of type float.
        ValueError: 
            If any of the input parameters are out of their respective valid ranges.
        FileNotFoundError: 
            If the required PCA file is not found.
        """

        output_name= self.output_name
        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        recompute  = self.recompute

        step = "association_analysis"

        # Check type of maf
        if not isinstance(maf, float):
             raise TypeError("maf should be of type float.")

        # Check type of mind
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check type of ci
        if not isinstance(ci, float):
            raise TypeError("ci should be of type float.")
        
        # Check if maf is in range
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 0.5")
        
        # Check if mind is in range
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1")
        
        # Check if hwe is in range
        if hwe < 0 or hwe > 1:
            raise ValueError("hwe should be between 0 and 1")
        
        # Check if ci is in range
        if ci <= 0 or ci >= 1:
            raise ValueError("ci should be between 0 and 1")
        
        # check if the PCA file exists
        if not os.path.exists(os.path.join(input_path, input_name+'.eigenvec')):
            raise FileNotFoundError(f"PCA file was not found: {os.path.join(input_path, input_name+'.eigenvec')}")

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2  # use all available cores
        else:
            max_threads = 10

        if recompute:

            # plink2 command to perform association analysis
            plink2_cmd = f"plink2 --bfile {os.path.join(input_path, input_name)} --adjust --ci {ci} --maf {maf} --mind {mind} --hwe {hwe} --covar {os.path.join(input_path, input_name+'.eigenvec')} --glm hide-covar omit-ref sex cols=+a1freq,+beta --out {os.path.join(results_dir, output_name+'_glm')} --threads {max_threads}"

            # execute plink command
            shell_do(plink2_cmd, log=True)

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

    def get_top_hits(self, maf:float)->dict:

        """
        Get the top hits from the association analysis with GCTA.
        """

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name
        recompute   = self.recompute

        # check type and range of maf
        if not isinstance(maf, float):
            raise TypeError("maf should be of type float.")
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 0.5")

        step = "get_top_hits"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # load results of association analysis and rename columns according GCTA requirements
        df = pd.read_csv(os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), sep="\t")
        rename = {
            '#CHROM'          : 'CHR',	
            'POS'             : 'POS',	
            'ID'              : 'SNP',
            'REF'             : 'A2',	
            'ALT'             : 'ALT',	
            'PROVISIONAL_REF?': 'PROVISIONAL_REF',	
            'A1'              : 'A1',	
            'OMITTED'         : 'OMITTED',	
            'A1_FREQ'         : 'freq',	
            'FIRTH?'          : 'FIRTH',	
            'TEST'            : 'TEST',	
            'OBS_CT'          : 'N',	
            'BETA'            : 'b',	
            'SE'              : 'se',	
            'L95'             : 'L95',	
            'U95'             : 'U95',	
            'Z_STAT'          : 'Z_STAT',	
            'P'               : 'p',	
            'ERRCODE'         : 'ERRCODE'
        }
        df = df.rename(columns=rename)

        # prepare .ma file
        df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

        df.to_csv(os.path.join(results_dir, 'cojo_file.ma'), sep="\t", index=False)

        del df

        if recompute:
            # gcta command
            gcta_cmd = f"gcta64 --bfile {os.path.join(input_path, input_name)} --maf {maf} --cojo-slct --cojo-file {os.path.join(results_dir, 'cojo_file.ma')}   --out {os.path.join(results_dir, 'cojo_file')} --thread-num {max_threads}"

            # execute gcta command
            shell_do(gcta_cmd, log=True)

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
    
    def annotate_top_hits(self, gtf_path:str=None, build:str='38')->dict:

        """
        Annotate the top hits from the association analysis.
        """

        import time

        results_dir = self.results_dir

        step = "annotate_hits"

        # load the data
        if os.path.exists(os.path.join(results_dir, 'cojo_file.jma.cojo')):
            df_hits = pd.read_csv(os.path.join(results_dir, 'cojo_file.jma.cojo'), sep="\t")
        else:
            FileExistsError("File cojo_file.jma not found in the results directory.")

        df_hits = df_hits[['Chr', 'SNP', 'bp']].copy()

        if gtf_path is None:
            gtf_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
            path_to_gz = os.path.join(os.path.abspath('..'), 'GCF_000001405.40_GRCh38.p14_genomic.gtf.gz')
            path_to_gtf= os.path.join(os.path.abspath('..'), 'GCF_000001405.40_GRCh38.p14_genomic.gtf')
            
            if os.path.exists(path_to_gz) is not True or os.path.exists(path_to_gtf) is not True:

                download_file(gtf_url, path_to_gz)
                
                with gzip.open(path_to_gz, 'rb') as f_in:
                     with open(path_to_gtf, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            gtf_path = path_to_gtf

        if (df_hits.empty is not True):
            variants_toanno = annogene(
                variants_toanno,
                id     ='SNP',
                chrom  ='Chr',
                pos    ='bp',
                log    =Log(),
                build  =build,
                source ="refseq",
                verbose=False,
                gtf_path=gtf_path
            ).rename(columns={"GENE":"GENENAME"})

        df_hits.to_csv(os.path.join(results_dir, 'top_hits_annotated.tsv'), sep="\t", index=False)

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
