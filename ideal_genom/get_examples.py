import logging

from pathlib import Path

from ideal_genom.Helpers import download_file, unzip_file_flat, extract_gz_file

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def get_trumpet_quantitative_example() -> Path:

    library_path = Path(__file__).resolve().parent.parent

    URL = "https://humandbs.dbcls.jp/files/hum0014/hum0014.v15.ht.v1.zip"
    FILENAME = "hum0014.v15.ht.v1.zip"
    LOCAL_PATH = library_path / "data" / "sumstats"

    local_filename = LOCAL_PATH / FILENAME

    example_path = LOCAL_PATH / "2019_BBJ_Height_autosomes_BOLT.txt.gz"

    if example_path.exists():
        logger.info(f"File already exists: {example_path}")
        return example_path
    
    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz = unzip_file_flat(local_filename, "hum0014.v15.ht.v1/2019_BBJ_Height_autosomes_BOLT.txt.gz", LOCAL_PATH, remove_zip=True)

    #logger.info(f"Decompressing file: {extracted_gz}")
    #extract_gz_file(extracted_gz, LOCAL_PATH, remove_gz=True)

    return extracted_gz