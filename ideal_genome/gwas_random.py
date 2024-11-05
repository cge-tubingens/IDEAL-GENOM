"""
Module to perform a GWAS analysis using a random effect model.
"""

import os
import pandas as pd

from ideal_genome.Helpers import shell_do, delete_temp_files

from gwaslab.bd_download import download_file
from gwaslab.util_in_get_sig import annogene
from gwaslab.g_Log import Log

class GWASrandom:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, preps_path:str, recompute:str=True) -> None:

        """
        Initialize the GWASrandom class.

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
        preps_path : str
            Path to the preparatory data.

        Returns:
        --------
        None
        """
        
       # check if paths are set
        if input_path is None or output_path is None or preps_path is None:
            raise ValueError("Values for input_path, output_path and dependables_path must be set upon initialization.")
        
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input path does not exist: {input_path}")
        if not os.path.exists(preps_path):
            raise FileNotFoundError(f"Dependables path does not exist: {preps_path}")
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
        
        if not isinstance(config_dict, dict):
            raise TypeError("config_dict should be of type dict.")
        if not isinstance(recompute, bool):
            raise TypeError("recompute should be of type bool.")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.config_dict = config_dict
        self.preps_path  = preps_path
        self.recompute   = recompute

        self.files_to_keep = []

        # create results folder
        self.results_dir = os.path.join(output_path, 'gwas_random')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create figures folder
        self.plots_dir = os.path.join(output_path, 'plots_random')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

        print("\033[1;32mAnalysis of GWAS data using a random effect model initialized.\033[0m")

        pass

    def prepare_aux_files(self)->dict:

        """
        Prepare auxiliary files for the GWAS analysis.
        """
        
        input_path  = self.input_path
        input_name  = self.input_name
        results_dir = self.results_dir
        output_name = self.output_name

        step = "prepare_aux_files"

        df_fam = pd.read_csv(os.path.join(input_path, input_name+'.fam'), sep=' ', header=None)

        df_pheno = df_fam[[df_fam.columns[0], df_fam.columns[1], df_fam.columns[5]]].copy()

        # recode phenotype
        df_pheno[5] = df_pheno[5]-1

        df_pheno.to_csv(os.path.join(results_dir, output_name+'_pheno.phen'), sep='\t', header=False, index=False)

        # recode sex
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

        """
        Compute the genetic relationship matrix (GRM) for the GWAS analysis using GCTA.
        """

        results_dir= self.results_dir
        prep_path  = self.preps_path
        output_name= self.output_name
        recompute  = self.recompute

        step = "compute_grm"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # gcta commands
        gcta_cmd1 = f"gcta64 --bfile {os.path.join(prep_path, output_name+'_LDpruned')} --make-grm --thread-num {max_threads} --out {os.path.join(results_dir, output_name+'_grm')}"

        gcta_cmd2 = f"gcta64 --grm {os.path.join(results_dir, output_name+'_grm')} --make-bK-sparse 0.05 --out {os.path.join(results_dir, output_name+'_sparse')}"

        # run gcta commands
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

        """
        Method to run the GWAS analysis using a random effect model.
        """

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name
        config_dict = self.config_dict
        preps_path  = self.preps_path
        recompute   = self.recompute

        maf = config_dict['maf']

        if not isinstance(maf, float):
            raise TypeError("maf should be of type float.")
        if maf < 0 or maf > 1:
            raise ValueError("maf should be between 0 and 1.")

        step = "run_gwas_random"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # gcta command
        gcta_cmd = f"gcta64 --bfile {os.path.join(input_path, input_name)} --fastGWA-mlm-binary --maf {maf}  --grm-sparse {os.path.join(results_dir, output_name+'_sparse')} --qcovar {os.path.join(preps_path, output_name+'_pca.eigenvec')} --covar {os.path.join(results_dir, output_name+'_sex.covar')} --pheno {os.path.join(results_dir, output_name+'_pheno.phen')} --out {os.path.join(results_dir,output_name+'_assocSparseCovar_pca_sex-mlm-binary')}--thread-num {max_threads}"

        # run gcta command
        shell_do(gcta_cmd, log=True)

        # rename columns for later use with GCTA
        df = pd.read_csv(os.path.join(results_dir, output_name+'_assocSparseCovar_pca_sex-mlm-binary--thread-num.fastGWA'), sep="\t")
        rename = {
            'CHR'     :'CHR',	
            'SNP'     :'SNP',
            'POS'     :'POS',	
            'A1'      :'A1', 
            'A2'      :'A2', 
            'N'       :'N', 
            'AF1'     :'freq', 
            'T'       :'T', 
            'SE_T'    :'SE_T', 
            'P_noSPA' :'P_noSPA',
            'BETA'    :'b', 
            'SE'      :'se', 
            'P'       :'p', 
            'CONVERGE':'CONVERGE'
        }
        df = df.rename(columns=rename)
        df.to_csv(os.path.join(results_dir, output_name+'_assocSparseCovar_pca_sex-mlm-binary--thread-num.fastGWA'), sep="\t", index=False)

        self.files_to_keep.append(output_name+'_assocSparseCovar_pca_sex-mlm-binary--thread-num.fastGWA')

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
    
    def get_top_hits(self)->dict:

        """
        Method to extract the top hits from the GWAS analysis.
        """

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name
        recompute   = self.recompute

        maf = self.config_dict['maf']

        if not isinstance(maf, float):
            raise TypeError("maf should be of type float.")
        if maf < 0 or maf > 1:
            raise ValueError("maf should be between 0 and 1.")

        step = "get_top_hits"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # load results of association analysis and rename columns according to GCTA requirements
        df = pd.read_csv(os.path.join(results_dir, output_name+'_assocSparseCovar_pca_sex-mlm-binary--thread-num.fastGWA'), sep="\t")
        rename = {
            'CHR'     :'CHR',	
            'SNP'     :'SNP',
            'POS'     :'POS',	
            'A1'      :'A1', 
            'A2'      :'A2', 
            'N'       :'N', 
            'AF1'     :'freq', 
            'T'       :'T', 
            'SE_T'    :'SE_T', 
            'P_noSPA' :'P_noSPA',
            'BETA'    :'b', 
            'SE'      :'se', 
            'P'       :'p', 
            'CONVERGE':'CONVERGE'
        }
        df = df.rename(columns=rename)        

        # prepare .ma file
        df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

        df.to_csv(os.path.join(results_dir, 'cojo_file.ma'), sep="\t", index=False)

        del df

        if recompute:
            # gcta command
            gcta_cmd = f"gcta64 --bfile {os.path.join(input_path, input_name)} --maf {maf} --cojo-slct --cojo-file {os.path.join(results_dir, 'cojo_file.ma')}   --out {os.path.join(results_dir, output_name+'_assocSparseCovar_pca_sex-mlm-binary-cojo')} --thread-num {max_threads}"

            # execute gcta command
            shell_do(gcta_cmd, log=True)

        self.files_to_keep.append(output_name+'_assocSparseCovar_pca_sex-mlm-binary-cojo.jma.cojo')
        self.files_to_keep.append(output_name+'_assocSparseCovar_pca_sex-mlm-binary-cojo.ldr.cojo')

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

        results_dir = self.results_dir
        output_name = self.output_name

        step = "annotate_hits"

        # load the data
        cojo_file_path = os.path.join(results_dir, output_name+'_assocSparseCovar_pca_sex-mlm-binary-cojo.jma.cojo')

        # check if .jma file exists
        if os.path.exists(cojo_file_path):
            df_hits = pd.read_csv(cojo_file_path, sep="\t")
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

        self.files_to_keep.append('top_hits_annotated.tsv')

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
