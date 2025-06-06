import os
import json

from ideal_genom.Helpers import arg_parser

from ideal_genom.preprocessing.preparatory import Preparatory
from ideal_genom.preprocessing.PostImputation import ProcessVCF, GetPLINK
from ideal_genom.gwas.gwas_glmm import GWASrandom

def analysis_pipe(params_dict: dict, data_dict: dict, steps_dict: dict) -> None:

    post_imputation_params = params_dict['post_imputation']
    preparatory_params     = params_dict['preparatory']
    gwas_glm_params        = params_dict['gwas_glm']
    gwas_glmm_params       = params_dict['gwas_glmm']

    if steps_dict['pos_imputation']:

        # instantiate the PostImputation class
        post_imp = PostImputation(
            input_path = data_dict['input_directory'],
            output_path= data_dict['output_directory'],
            output_name= data_dict['output_prefix'],
            dependables= data_dict['dependables'],
        )

        # pipeline steps
        post_imp_steps = {
            'unzip_chrom'  : (post_imp.execute_unzip_chromosome_files, {'password': post_imputation_params['zip_password']}),
            'filter_by_R2' : (post_imp.execute_filter_variants, {'r2_threshold': post_imputation_params['r2_threshold']}),
            'normalize'    : (post_imp.execute_normalize_vcf, {}),
            'normalize_ref': (post_imp.execute_normalize_with_reference, {'reference_genome': post_imputation_params['ref_vcf']}),
            'index'        : (post_imp.execute_index_vcf, {}),
            'annotate'     : (post_imp.execute_annotate_vcf, {'annotations_file': post_imputation_params['ref_annotation']}),
            'concatenate'  : (post_imp.execute_concat_vcf, {}),
            'get_plink'    : (post_imp.get_plink_files, {}),
            'clean_up'     : (post_imp.cleanup, {})
        }

        step_description = {
            'unzip_chrom' : 'Unzip chromosome files',
            'filter_by_R2': 'Filter imputed variants by R2',
            'normalize'   : 'Normalize VCF files',
            'normalize_ref': 'Normalize VCF files with reference genome',
            'index'       : 'Index VCF files',
            'annotate'    : 'Annotate VCF files',
            'concatenate' : 'Concatenate VCF files',
            'get_plink'   : 'Get PLINK files',
            'clean_up'    : 'Remove intermediate files'
        }

        for name, (func, params) in post_imp_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            func(*params)

        print("\033[1mPost-imputation steps completed.\033[0m")

    if steps_dict['preparatory']:
        # instantiate the Preparatory class
        preps = Preparatory(
            input_path  =os.path.join(data_dict['output_directory'], 'post_imputation', 'analysis_ready'),
            input_name  =data_dict['input_prefix'],
            output_path =data_dict['output_directory'],
            output_name =data_dict['output_prefix'],
            dependables =data_dict['dependables'],
        )

        # pipeline steps
        prep_steps = {
            'ld_prune': (preps.execute_ld_prunning, {
                'maf'     : preparatory_params['maf'], 
                'geno'    : preparatory_params['geno'],
                'hwe'     : preparatory_params['hwe'], 
                'ind_pair': preparatory_params['ind_pair'],
            }),
            'pca': (preps.execute_pc_decomposition, {
                'pca': preparatory_params['pca']
            }),
        }

        step_description = {
            'ld_prune': 'Linkage Disequilibrium Prunning',
            'pca'     : 'Principal Component Analysis'
        }

        for name, (func, params) in prep_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            func(**params)

        print("\033[1mPreprocessing steps completed.\033[0m")

    if steps_dict['gwas_glm']:
        # instantiate the GWASfixed class
        gwas_glm = GWASfixed(
            input_path = os.path.join(data_dict['output_directory'], 'post_imputation', 'analysis_ready'),
            input_name = data_dict['input_prefix'],
            output_path= data_dict['output_directory'],
            output_name= data_dict['output_prefix'],
        )

        # pipeline steps
        gwas_steps = {
            'train_model': (gwas_glm.fixed_model_association_analysis, {
                'maf' :gwas_glm_params['maf'], 
                'mind':gwas_glm_params['mind'], 
                'hwe' :gwas_glm_params['hwe'], 
                'ci'  :gwas_glm_params['ci']
            }),
            'top_hits'   : (gwas_glm.get_top_hits, {'maf':gwas_glm_params['maf']}),
        }
        
        step_description = {
            'train_model': 'Train the model',
            'top_hits'   : 'Get top hits'
        }
        
        for name, (func, params) in gwas_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            func(**params)

        print("\033[1mGLM GWAS steps completed.\033[0m")

    if steps_dict['gwas_glmm']:
        # instantiate the GWASrandom class
        gwas_glmm = GWASrandom(
            input_path = os.path.join(data_dict['output_directory'], 'post_imputation', 'analysis_ready'),
            input_name = data_dict['input_prefix'],
            output_path= data_dict['output_directory'],
            output_name= data_dict['output_prefix'],
        )

        # pipeline steps
        gwas_steps = {
            'aux_files'  : (gwas_glmm.prepare_aux_files, {}),
            'compute_grm': (gwas_glmm.compute_grm, {}),
            'run_gwas'   : (gwas_glmm.run_gwas_random, {'maf' :gwas_glmm_params['maf']}),
            'top_hits'   : (gwas_glmm.get_top_hits, {'maf' :gwas_glmm_params['maf']})
        }

        step_description = {
            'aux_files'  : 'Prepare auxiliary files',
            'compute_grm': 'Compute genetic relationship matrix',
            'run_gwas'   : 'Run GWAS',
            'top_hits'   : 'Get top hits'
        }

        for name, (func, params) in gwas_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            func(**params)

        print("\033[1mGLLM GWAS steps completed.\033[0m")

    pass

def execute_main()->str:

    args = arg_parser()
    args_dict = vars(args)

    params_path = args_dict['path_params']
    data_path   = args_dict['file_folders']
    steps_path  = args_dict['steps']

    # check path to config files
    if not os.path.exists(data_path):
        raise FileNotFoundError("Configuration file with path to data and analysis results cannot be found.")
    
    if not os.path.exists(params_path):
        raise FileNotFoundError("Configuration file with pipeline parameters cannot be found.")
    
    if not os.path.exists(steps_path):
        raise FileNotFoundError("Configuration file with pipeline steps cannot be found.")
    
    # open config file
    with open(data_path, 'r') as file:
        data_dict = json.load(file)

    with open(params_path, 'r') as file:
        params_dict = json.load(file)

    with open(steps_path, 'r') as file:
        steps_dict = json.load(file)

    analysis_pipe(params_dict, data_dict, steps_dict)

    return "Analysis pipeline completed."

if __name__ == '__main__':
    execute_main()
