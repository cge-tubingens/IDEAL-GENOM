import os
import json

from ideal_genom.Helpers import arg_parser

from ideal_genom.preprocessing.preparatory import Preparatory
from ideal_genom.preprocessing.post_imputation import PostImputation
from ideal_genom.gwas.gwas_fixed import GWASfixed
from ideal_genom.gwas.gwas_random import GWASrandom

def analysis_pipe(params_dict:dict, data_dict:dict, steps_dict:dict)->None:

    if steps_dict['prep_ds']:
        # instantiate the PrepDS class
        prep = PrepDS(
            input_path      =data_dict['input_directory'],
            input_name      =data_dict['input_prefix'],
            output_path     =data_dict['output_directory'],
            output_name     =data_dict['output_prefix'],
            config_dict     =params_dict,
            dependables_path=data_dict['dependables_directory'],
        )

        # pipeline steps
        prep_steps = {
            'ld_prune': prep.exclude_high_ld_hla,
            'pca'     : prep.pca_decomposition
        }

        for step in prep_steps.keys():
            prep_steps[step]()

        print("\033[1mPreprocessing steps completed.\033[0m")

    if steps_dict['gwas_fixed']:
        # instantiate the GWASfixed class
        gwas_f = GWASfixed(
            input_path = data_dict['input_directory'],
            input_name = data_dict['input_prefix'],
            output_path= data_dict['output_directory'],
            output_name= data_dict['output_prefix'],
            config_dict= params_dict,
            preps_path = os.path.join(data_dict['output_directory'], 'preparatory'),
        )

        # pipeline steps
        gwas_f_steps = {
            'train_model': gwas_f.fixed_model_association_analysis,
            'top_hits'   : gwas_f.get_top_hits,
            'annotate'   : gwas_f.annotate_top_hits,
        }

        for step in gwas_f_steps.keys():
            gwas_f_steps[step]()

        print("\033[1mFixed model GWAS steps completed.\033[0m")

    if steps_dict['gwas_random']:
        # instantiate the GWASrandom class
        gwas_r = GWASrandom(
            input_path = data_dict['input_directory'],
            input_name = data_dict['input_prefix'],
            output_path= data_dict['output_directory'],
            output_name= data_dict['output_prefix'],
            config_dict= params_dict,
            preps_path = os.path.join(data_dict['output_directory'], 'preparatory'),
        )

        # pipeline steps
        gwas_r_steps = {
            'preparatory': gwas_r.prepare_aux_files,
            'grm'        : gwas_r.compute_grm,
            'random'     : gwas_r.run_gwas_random,
            'top_hits'   : gwas_r.get_top_hits,
            'annotate'   : gwas_r.annotate_top_hits,
        }
        for step in gwas_r_steps.keys():
            gwas_r_steps[step]()

        print("\033[1mRandom model GWAS steps completed.\033[0m")

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
