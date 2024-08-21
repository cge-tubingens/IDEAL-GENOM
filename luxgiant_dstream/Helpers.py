import subprocess
import sys
import argparse
import os
import shutil

def shell_do(command, print_cmd=False, log=False, return_log=False, err=False):

    """
    From GenoTools
    """
    
    if print_cmd:
        print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res = subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')
    
def arg_parser()->dict:

    # define parser
    parser = argparse.ArgumentParser(description='Adresses to configuration files')

    # parameters of quality control
    parser.add_argument('--path_params', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing genotype quality control parameters.')

    # path to data and names of files
    parser.add_argument('--file_folders', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing folder names and locations for genotype quality control data.')

    # path to steps of the pipeline to be executed
    parser.add_argument('--steps', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing the pipeline steps to be executed.')

    # order of the pipeline
    parser.add_argument('--pca-first', type=str, nargs='?', default=None, const=None, help='boolean that determines if the pca analysis is done before or after the sample qc.')

    # execute kingship analysis instead of ibd
    parser.add_argument('--use-kingship', type=str, nargs='?', default=None, const=None, help='boolean that determines if the relatives and duplicates prune is run using kingship or IBD.')

    # parse args and turn into dict
    args = parser.parse_args()

    return args

def delete_temp_files(files_to_keep:list, path_to_folder:str)->None:

    """
    Function to delete temporary files that were created during the pipeline execution. Moreover, it creates a directory called 'log_files' to save al `.log` files originated from the pipeline execution.

    Parameters
    ----------
    files_to_keep: list
        list of strings where its elements are the names of files and folders that should be kept.
    path_to_folder: str
        full path to the folder where the temporary files are located.
    """

    for file in os.listdir(path_to_folder):
        file_split = file.split('.')
        if file_split[-1]!='log' and file not in files_to_keep and file_split[-1]!='hh':
            if os.path.isfile(os.path.join(path_to_folder, file)):
                os.remove(
                    os.path.join(path_to_folder, file)
                )
        
    # create log folder for dependables
    logs_dir = os.path.join(path_to_folder, 'log_files')
    if not os.path.exists(logs_dir):
        os.mkdir(logs_dir)

    for file in os.listdir(path_to_folder):
        if file.split('.')[-1]=='log' or file.split('.')[-1]=='hh':
            shutil.move(
                os.path.join(path_to_folder, file),
                os.path.join(logs_dir, file)
            )

    pass
