import subprocess
import sys
import argparse
import requests
import os
import shutil
import zipfile
import gzip

from pathlib import Path

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
    
def arg_parser() -> dict:

    # define parser
    parser = argparse.ArgumentParser(description='Adresses to configuration files')

    # parameters of quality control
    parser.add_argument('--path-params', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing genotype quality control parameters.')

    # path to data and names of files
    parser.add_argument('--file-folders', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing folder names and locations for genotype quality control data.')

    # path to steps of the pipeline to be executed
    parser.add_argument('--steps', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing the pipeline steps to be executed.')

    # parse args and turn into dict
    args = parser.parse_args()

    return vars(args)

def download_file(url:str, local_filename: Path) -> None:

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    return

def unzip_file_flat(in_file: Path, target_file: str, out_dir: Path, remove_zip: bool = False) -> Path:
    """Extracts a specific file from a ZIP archive, decompresses it if it's a .gz file, and optionally deletes original files.

    Args:
        in_file (str): Path to the ZIP file.
        target_file (str): The file inside the ZIP to extract.
        out_dir (str): Directory where the extracted file will be saved.
        remove_zip (bool): If True, delete the original ZIP file after extraction.
        remove_gz (bool): If True, delete the .gz file after decompression.

    Returns:
        Path: Path to the final extracted file.
    """
    in_file = Path(in_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

    extracted_gz_path = out_dir / Path(target_file).name  # Target extracted .gz file

    try:
        with zipfile.ZipFile(in_file, "r") as zip_ref:
            if target_file in zip_ref.namelist():
                # Extract the .gz file from ZIP
                with zip_ref.open(target_file) as source, open(extracted_gz_path, "wb") as dest:
                    dest.write(source.read())
                print(f"Extracted: {extracted_gz_path}")
            else:
                print(f"File {target_file} not found in the archive.")
                return Path()

        # Optionally delete the ZIP file
        if remove_zip:
            in_file.unlink()
            print(f"Deleted ZIP file: {in_file}")

        return extracted_gz_path

    except zipfile.BadZipFile:
        print(f"Error: {in_file} is not a valid ZIP file.")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return Path()  # Return None if extraction fails

def extract_gz_file(gz_file: Path, out_dir: Path, remove_gz: bool = False) -> Path:
    """Extracts a .gz file and saves the decompressed content in the same directory.

    Args:
        gz_file (str): Path to the .gz file.
        out_dir (str): Directory where the decompressed file will be saved.
        remove_gz (bool): If True, delete the .gz file after extraction.

    Returns:
        Path: Path to the extracted file.
    """
    gz_file = Path(gz_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

    # Define the decompressed file path (removes .gz extension)
    decompressed_file = out_dir / gz_file.stem

    try:
        with gzip.open(gz_file, "rb") as f_in, open(decompressed_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)  # Copy content from .gz to uncompressed file
        print(f"Decompressed: {decompressed_file}")

        if remove_gz:
            gz_file.unlink()  # Delete the .gz file
            print(f"Removed original .gz file: {gz_file}")

    except Exception as e:
        print(f"Error extracting {gz_file}: {e}")

    return decompressed_file


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
