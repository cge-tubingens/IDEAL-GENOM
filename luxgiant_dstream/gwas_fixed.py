"""
Module to perform a GWAS analysis using a fixed model.
"""

import os
import pandas as pd

from luxgiant_dstream.Helpers import shell_do, delete_temp_files
from luxgiant_dstream.plots import manhattan_plot, qq_plot
from luxgiant_dstream.annotate_tools import get_variant_context

class GWASfixed:

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
        self.results_dir = os.path.join(output_path, 'gwas_fixed')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create figures folder
        self.plots_dir = os.path.join(output_path, 'plots')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

        pass

    def fixed_model_association_analysis(self)->dict:

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

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2  # use all available cores
        else:
            max_threads = 10

        # Run plink2 to perform association analysis
        plink2_cmd = f"plink2 --bfile {os.path.join(input_path, input_name)} --adjust --ci {ci} --maf {maf} --mind {mind} --hwe {hwe} --covar {os.path.join(preps_dir, output_name+'_pca.eigenvec')} --glm hide-covar omit-ref sex cols=+a1freq,+beta --out {os.path.join(results_dir, output_name+'_glm')} --threads {max_threads}"

        shell_do(plink2_cmd, log=True)

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

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name

        maf = self.config_dict['maf']

        step = "get_top_hits"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # CHANGE FOR FINAL VERSION
        # load results of association analysis
        df = pd.read_csv(os.path.join(results_dir, output_name+'_glm1.PHENO1.glm.logistic.hybrid'), sep="\t")

        # delete for final version
        self.files_to_keep.append(output_name+'_glm1.PHENO1.glm.logistic.hybrid')

        # prepare .ma file
        df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

        df.to_csv(os.path.join(results_dir, 'cojo_file.ma'), sep="\t", index=False)

        del df

        # gcta command
        gcta_cmd = f"gcta64 --bfile {os.path.join(input_path, input_name)} --maf {maf} --cojo-slct --cojo-file {os.path.join(results_dir, 'cojo_file.ma')}   --out {os.path.join(results_dir, 'cojo_file')} --thread-num {max_threads}"

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

            time.sleep(1.5)

        df_hits = df_hits[['SNP', 'GENE']].copy()

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

        plots_dir  = self.plots_dir
        results_dir= self.results_dir
        output_name= self.output_name

        annotate = self.config_dict['annotate']

        step = "draw_plots"

        df_gwas = pd.read_csv(
            os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), 
            sep="\t",
            usecols=['SNP', 'CHR', 'p']
        )

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
