"""
Module to perform a GWAS analysis using a fixed model.
"""

import gzip
import os
import shutil

import pandas as pd

from ideal_genome.Helpers import shell_do, delete_temp_files

from gwaslab.bd_download import download_file
from gwaslab.util_in_get_sig import annogene
from gwaslab.g_Log import Log

class GWASfixed:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, preps_path:str, dependables:str) -> None:

        """
        Initialize the GWASfixed class.

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
        if input_path is None or output_path is None:
            raise ValueError("Values for input_path, output_path and dependables_path must be set upon initialization.")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.config_dict = config_dict
        self.preps_path  = preps_path
        self.dependables = dependables

        self.files_to_keep = []
        self.compare_gwas_fixed_file_name = None
        self.compare_gwas_fixed_highlights= None

        # create results folder
        self.results_dir = os.path.join(output_path, 'gwas_fixed')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create figures folder
        self.plots_dir = os.path.join(output_path, 'plots_fixed')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

        print("\033[1;32mAnalysis of GWAS data using a fixed model initialized.\033[0m")

        pass

    def fixed_model_association_analysis(self)->dict:

        """
        Perform association analysis using a fixed model with PLINK 2.0.
        """

        output_name= self.output_name
        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        preps_dir  = self.preps_path

        maf = self.config_dict['maf']
        mind= self.config_dict['mind']
        hwe = self.config_dict['hwe']
        ci  = self.config_dict['ci']

        step = "association_analysis"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2  # use all available cores
        else:
            max_threads = 10

        # plink2 command to perform association analysis
        plink2_cmd = f"plink2 --bfile {os.path.join(input_path, input_name)} --adjust --ci {ci} --maf {maf} --mind {mind} --hwe {hwe} --covar {os.path.join(preps_dir, output_name+'_pca.eigenvec')} --glm hide-covar omit-ref sex cols=+a1freq,+beta --out {os.path.join(results_dir, output_name+'_glm')} --threads {max_threads}"

        # execute plink command
        shell_do(plink2_cmd, log=True)

        # rename columns for later use with GCTA
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
        df.to_csv(os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), sep="\t", index=False)

        self.files_to_keep.append(output_name+'_glm.PHENO1.glm.logistic.hybrid')
        self.files_to_keep.append(output_name+'_glm.PHENO1.glm.logistic.hybrid.adjusted')

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

    def get_top_hits(self)->dict:

        """
        Get the top hits from the association analysis with GCTA.
        """

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name

        maf = self.config_dict['maf']

        step = "get_top_hits"

        # compute the number of threads to use
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # load results of association analysis
        df = pd.read_csv(os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), sep="\t")

        # prepare .ma file
        df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

        df.to_csv(os.path.join(results_dir, 'cojo_file.ma'), sep="\t", index=False)

        del df

        # gcta command
        gcta_cmd = f"gcta64 --bfile {os.path.join(input_path, input_name)} --maf {maf} --cojo-slct --cojo-file {os.path.join(results_dir, 'cojo_file.ma')}   --out {os.path.join(results_dir, 'cojo_file')} --thread-num {max_threads}"

        # execute gcta command
        shell_do(gcta_cmd, log=True)

        self.files_to_keep.append('cojo_file.jma.cojo')
        self.files_to_keep.append('cojo_file.ldr.cojo')

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
    
    def annotate_top_hits(self)->dict:

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

        for k in range(df_hits.shape[0]):
            # get variant context
            chr = df_hits.loc[k, 'Chr']
            pos = df_hits.loc[k, 'bp']

            context = get_variant_context(chr, pos)

            if context is None:
                context = 'NA'
            df_hits.loc[k, 'GENE'] = context[0]

            # sleep for 1.5 seconds
            time.sleep(1.5)

        df_hits = df_hits[['SNP', 'GENE']].copy()

        # save the annotated data
        df_hits.to_csv(os.path.join(results_dir, 'snps_annotated.csv'), sep="\t", index=False)

        self.files_to_keep.append('snps_annotated.csv')

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

    def plot_drawings(self)->dict:

        """
        Draw Manhattan plot and QQ plot for the GWAS analysis.
        """

        plots_dir  = self.plots_dir
        results_dir= self.results_dir
        output_name= self.output_name

        annotate = self.config_dict['annotate']

        step = "draw_plots"

        # load GWAS results
        df_gwas = pd.read_csv(
            os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), 
            sep="\t",
            usecols=['SNP', 'CHR', 'p']
        )

        # draw manhattan plot
        if annotate:
            df_annot = pd.read_csv(os.path.join(results_dir, "snps_annotated.csv"), sep="\t")

            mann_plot = manhattan_plot(
                df_gwas    =df_gwas,
                df_annot   =df_annot,
                plots_dir  =plots_dir,
                annotate   =True
            )
        else:
            mann_plot = manhattan_plot(
                df_gwas    =df_gwas,
                df_annot   =df_annot,
                plots_dir  =plots_dir,
                annotate   =False
            )

        # draw QQ plot
        QQ_plot = qq_plot(df_gwas, plots_dir)

        delete_temp_files(self.files_to_keep, results_dir)

        # report
        process_complete = (mann_plot & QQ_plot)

        outfiles_dict = {
            'plink_out': plots_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }
        
        return out_dict

    def one_cohort_comparison(self)->None:

        """
        Compare the results of the GWAS analysis with other datasets.
        """

        output_name  = self.output_name
        plots_dir    = self.plots_dir
        results_dir  = self.results_dir
        dependables  = self.dependables
        gwas_compare = self.compare_gwas_fixed_file_name
        gwas_highlts = self.compare_gwas_fixed_highlights

        if gwas_compare is None or gwas_compare == '':
            print("\033[1;31mNo comparison file provided. Skipping comparison.\033[0m")
            pass
        if gwas_highlts is None or gwas_highlts == '':
            df_bottom_annot=None
            print("\033[1;31mNo highlights file provided\033[0m")
        else:
            df_bottom_annot = pd.read_csv(
            os.path.join(dependables, gwas_highlts),
            sep="\t"
        )

        # load the data
        # it is expected that the data comes with columns 'SNP', 'CHR', 'bp', 'p'
        df_gwas_bottom = pd.read_csv(
            os.path.join(dependables, gwas_compare),
            sep="\t",
            usecols=['SNP', 'CHR', 'POS', 'p']
        )

        df_gwas_top = pd.read_csv(
            os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'),
            sep="\t",
            usecols=['SNP', 'CHR', 'POS', 'p']
        )
        if os.path.isfile(os.path.join(results_dir, 'snps_annotated.csv')):
            df_top_annot = pd.read_csv(os.path.join(results_dir, 'snps_annotated.csv'), sep="\t")
        else:
            df_top_annot = None

        miami = miami_plot(
            df_top       =df_gwas_top,
            top_higlights=df_top_annot,
            df_bottom    =df_gwas_bottom,
            bottom_higlights=df_bottom_annot,
            plots_dir    =plots_dir
        )

        return miami
    
    def create_trumpet_plot(self)->None:

        """
        Create a trumpet plot for the GWAS analysis.
        """

        input_path = self.input_path
        input_name = self.input_name
        output_name= self.output_name
        plots_dir  = self.plots_dir
        results_dir= self.results_dir
        preps_path = self.preps_path

        prevalence = self.config_dict['prevalence']
        maf = self.config_dict['maf']

        # plink command to compute MAF
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --freq --maf {maf} --out {os.path.join(results_dir, output_name)}"

        # execute plink command
        shell_do(plink_cmd, log=True)

        # load gwas data
        df_gwas = pd.read_csv(
            os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), 
            sep="\t",
        )
        # ensure column 'CHR' is of type int
        df_gwas['CHR'] = df_gwas['CHR'].astype(int)

        # load file with 'MAF'
        df_freq = pd.read_csv(
            os.path.join(results_dir, output_name+'.frq'),
            sep="\s+",
            usecols=['SNP', 'MAF']
        )

        # load .fam file to compute the number of cases and controls
        df_fam = pd.read_csv(
            os.path.join(preps_path, output_name+'_LDpruned.fam'),
            sep="\s+",
            header=None
        )

        counts  = df_fam[5].value_counts()
        ncase   = counts[2]
        ncontrol= counts[1]

        # estimate prevalence from ncase and ncontrol
        if prevalence is None:
            prevalence = ncase / (ncase + ncontrol)

        # load SNPs to highlight and annotate
        if os.path.isfile(os.path.join(results_dir, 'snps_annotated.csv')):
            df_annot = pd.read_csv(os.path.join(results_dir, 'snps_annotated.csv'), sep="\t")
        else:
            df_annot = None
        
        # merge the data
        df = pd.merge(df_gwas, df_freq, on='SNP', how='inner')

        del df_gwas, df_freq, df_fam

        # rename columns to match GWASlab format
        rename_cols = {
            'SNP': 'rsID',
            'CHR': 'CHR',
            'POS': 'POS',
            'A1' : 'EA',
            'A2' : 'NEA',
            'freq': 'EAF',
            'b'  : 'BETA',
            'se' : 'SE',
            'Z_STAT': 'Z',
            'p'  : 'P',
            'MAF': 'MAF'
        }

        df = df.rename(columns=rename_cols)

        # create a column with the SNP ID CHR:POS:EA:NEA
        df['SNPID'] = df['CHR'].astype(str) +':'+ df['POS'].astype(str) +':'+ df['EA'].astype(str) +':'+ df['NEA'].astype(str)
        
        # find SNP IDs to highlight
        if df_annot is not None:
            df_annot = df_annot.merge(
                df[['rsID', 'SNPID', 'CHR', 'POS']], 
                left_on='SNP', 
                right_on='rsID', 
                how='inner'
            )
            df_annot = df_annot.drop(columns=['SNP'])

        # power curves thresholds
        ts=[0.2,0.4,0.6,0.8]

        trumpet = draw_trumpet_plot(
            df_gwas    =df,
            plot_dir   =plots_dir,
            mode       ="b",
            xscale     ='nonlog',
            n_matrix   =3000,
            ts         =ts,
            ncase      =ncase, 
            ncontrol   =ncontrol, 
            prevalence =prevalence,
            figargs    ={"figsize": (20, 16), "dpi": 400},
            font_family='DejaVu Sans',
            highlight  =df_annot['SNPID'].tolist(),
            highlight_windowkb = 0.01,
            anno       ="GENENAME",
            build      ="38",
            anno_set   =df_annot['SNPID'].tolist(),
            anno_style ="expand",
            ylim       =(-df['BETA'].abs().max()*1.5, df['BETA'].abs().max()*1.5),
            xlim       =(df['MAF'].min()*0.5,0.52),
            anno_args  ={'bbox':dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.8)}
        )

        return trumpet
