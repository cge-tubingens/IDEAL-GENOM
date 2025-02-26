import logging
import requests

import pandas as pd

from pathlib import Path

from ideal_genom.Helpers import download_file, unzip_file_flat, extract_gz_file

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def get_trumpet_quantitative_example() -> Path:

    library_path = Path(__file__).resolve().parent.parent

    URL = "https://humandbs.dbcls.jp/files/hum0014/hum0014.v15.ht.v1.zip"
    FILENAME = "hum0014.v15.ht.v1.zip"
    LOCAL_PATH = library_path / "data" / "sumstats"

    LOCAL_PATH.mkdir(parents=True, exist_ok=True)

    local_filename = LOCAL_PATH / FILENAME

    example_path = LOCAL_PATH / "2019_BBJ_Height_autosomes_BOLT.txt"

    if example_path.exists():
        logger.info(f"File already exists: {example_path}")
        return example_path
    
    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz = unzip_file_flat(local_filename, "hum0014.v15.ht.v1/2019_BBJ_Height_autosomes_BOLT.txt.gz", LOCAL_PATH, remove_zip=True)

    logger.info(f"Decompressing file: {extracted_gz}")
    uncompressed_file = extract_gz_file(extracted_gz, LOCAL_PATH, remove_gz=True)

    return uncompressed_file

def get_top_loci_trumpet_quantitative() -> Path:

    library_path = Path(__file__).resolve().parent.parent

    url = r"https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12276-5/MediaObjects/41467_2019_12276_MOESM5_ESM.xlsx"
    output_filename = "2019_BBJ_Height_autosomes_BOLT_loci.xlsx"

    output_path = library_path / "data" / "sumstats" / output_filename

    output_csv = output_path.with_suffix('.csv')
    if output_csv.exists():
        logger.info(f"File already exists: {output_csv}")
        return output_csv

    response = requests.get(url, stream=True)  # Stream to handle large files
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
                file.write(chunk)
        logger.info(f"Downloaded file: {output_path}")
    else:
        logger.info(f"Failed to download file. Status code: {response.status_code}")

    df_top = pd.read_excel(output_path, engine='openpyxl')

    df_top = df_top[['Unnamed: 7', 'Unnamed: 9']].copy()
    df_top = df_top.rename(columns={'Unnamed: 7': 'Variants', 'Unnamed: 9': 'New_Locus'}).dropna().iloc[1:, :].copy()

    mask_autosome = df_top['Variants'].str.contains('X')

    df_top = df_top[~mask_autosome].reset_index(drop=True)

    df_top.to_csv(output_path.with_suffix('.csv'), index=False, sep='\t')
    logger.info(f"Saved top hits to: {output_path.with_suffix('.csv')}")

    output_path.unlink()

    return output_csv

def get_top_cond_trumpet_quantitative() -> Path:

    library_path = Path(__file__).resolve().parent.parent

    url = r"https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12276-5/MediaObjects/41467_2019_12276_MOESM6_ESM.xlsx"
    output_filename = "2019_BBJ_Height_autosomes_BOLT_cond.xlsx"

    output_path = library_path / "data" / "sumstats" / output_filename

    output_csv = output_path.with_suffix('.csv')
    if output_csv.exists():
        logger.info(f"File already exists: {output_csv}")
        return output_csv

    response = requests.get(url, stream=True)  # Stream to handle large files
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
                file.write(chunk)
        logger.info(f"Downloaded file: {output_path}")
    else:
        logger.info(f"Failed to download file. Status code: {response.status_code}")

    df_top = pd.read_excel(output_path, engine='openpyxl', header=1)

    df_top = df_top[['rsID', 'Candiate gene(s)', 'CHRa', 'POSa']].copy()
    df_top = df_top.rename(columns={'Candiate gene(s)': 'Gene', 'CHRa': 'CHR', 'POSa': 'POS'}).dropna()

    mask_autosome = (df_top['CHR']!='X')

    df_top = df_top[mask_autosome].reset_index(drop=True)

    df_top.to_csv(output_path.with_suffix('.csv'), index=False, sep='\t')
    logger.info(f"Saved top hits to: {output_path.with_suffix('.csv')}")

    output_path.unlink()

    return output_csv

def get_trumpet_binary_example() -> Path:

    library_path = Path(__file__).resolve().parent.parent

    URL = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278092/harmonised/GCST90278092.h.tsv.gz"
    FILENAME = "GCST90278092.h.tsv.gz"
    LOCAL_PATH = library_path / "data" / "sumstats"

    LOCAL_PATH.mkdir(parents=True, exist_ok=True)

    local_filename = LOCAL_PATH / FILENAME

    example_path = LOCAL_PATH / "GCST90278092.h.tsv"

    if example_path.exists():
        logger.info(f"File already exists: {example_path}")
        return example_path
    
    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

#
    logger.info(f"Decompressing file: {local_filename}")
    uncompressed_file = extract_gz_file(local_filename, LOCAL_PATH, remove_gz=True)

    return uncompressed_file

def get_bmi_japanese_gwas() -> tuple:

    library_path = Path(__file__).resolve().parent.parent

    URL = r"https://humandbs.dbcls.jp/files/hum0014/hum0014.v6.158k.v1.zip"
    FILENAME = "hum0014.v6.158k.v1.zip"
    LOCAL_PATH = library_path / "data" / "sumstats"

    LOCAL_PATH.mkdir(parents=True, exist_ok=True)

    local_filename = LOCAL_PATH / FILENAME

    female_path = LOCAL_PATH / "Female_2017_BMI_BBJ_autosome.txt"
    male_path = LOCAL_PATH / "Male_2017_BMI_BBJ_autosome.txt"

    if female_path.exists() and male_path.exists():
        logger.info(f"Files {female_path} and {male_path} already exist")
        return female_path, male_path

    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz_f = unzip_file_flat(local_filename, "hum0014.v6.158k.v1/Female_2017_BMI_BBJ_autosome.txt.gz", LOCAL_PATH, remove_zip=False)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz_m = unzip_file_flat(local_filename, "hum0014.v6.158k.v1/Male_2017_BMI_BBJ_autosome.txt.gz", LOCAL_PATH, remove_zip=True)

    logger.info(f"Decompressing file: {extracted_gz_f}")
    uncompressed_file_f = extract_gz_file(extracted_gz_f, LOCAL_PATH, remove_gz=True)

    logger.info(f"Decompressing file: {extracted_gz_m}")
    uncompressed_file_m = extract_gz_file(extracted_gz_m, LOCAL_PATH, remove_gz=True)

    return uncompressed_file_f, uncompressed_file_m