import os
import json
import pandas as pd

from luxgiant_dstream.Helpers import arg_parser

from luxgiant_dstream.prep_ds import PrepDS
from luxgiant_dstream.gwas_fixed import GWASfixed
from luxgiant_dstream.gwas_random import GWASrandom

def analysis_pipe(params_dict:dict, data_dict:dict, steps_dict)->None:

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

        print("Preprocessing steps completed.")