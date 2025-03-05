import requests
import os
import gzip
import shutil
import re
import logging

import pandas as pd

from bs4 import BeautifulSoup
from gtfparse import read_gtf
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
        
class ReferenceDataFetcher:

    def __init__(self, base_url: str, build: str, source: str, destination_folder: str=None) -> None:

        self.build = build
        self.source = source
        self.base_url = base_url
        self.destination_folder = destination_folder

        self.latest_url = None
        self.gz_file = None
        self.gtf_file = None

        pass

    def get_latest_release(self) -> None:
        """Determine the specific URL for fetching data."""
        raise NotImplementedError("Subclasses must implement this method.")

    def download_latest(self) -> str:
        """
        Downloads the latest file from `self.latest_url` to `self.destination_folder`.

        Raises:
        -------
            - AttributeError: If `self.latest_url` is not set.
            - requests.exceptions.RequestException: If the HTTP request fails.
        """

        if not getattr(self, 'latest_url', None):
            raise AttributeError("`self.latest_url` is not set. Call `get_latest_release` first.")

        self.destination_folder = self.get_destination_folder()

        file_path = self.destination_folder / Path(self.latest_url).name

        if file_path.exists():

            self.gz_file = str(file_path)
            logger.info(f"File already exists: {file_path}")
            
            return str(file_path)

        self._download_file(self.latest_url, file_path)

        self.gz_file = str(file_path)

        return str(file_path)

    def get_destination_folder(self) -> Path:

        """Determine the destination folder for downloads."""

        if self.destination_folder:
            destination = Path(self.destination_folder)
        else:
            # Determine project root and default `data` directory
            project_root = Path(__file__).resolve().parent.parent
            destination = project_root / "data" / f"{self.source}_latest"

        destination.mkdir(parents=True, exist_ok=True)

        return destination

    def _download_file(self, url: str, file_path: Path) -> None:

        """Download a file from the given URL and save it to `file_path`."""

        with requests.get(url, stream=True) as response:
            response.raise_for_status()

            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

        logger.info(f"Downloaded file to: {file_path}")
    
    def unzip_latest(self) -> str:
        """Unzips the latest downloaded file and stores it as a GTF file."""

        if not getattr(self, 'latest_url', None):
            raise AttributeError("`self.latest_url` is not set. Call `get_latest_release` first.")

        self.destination_folder = self.get_destination_folder()

        if not hasattr(self, 'gz_file') or not Path(self.gz_file).is_file():
            raise FileNotFoundError("Reference file not found")

        gtf_file = self.destination_folder / (Path(self.gz_file).stem)  # Removes .gz extension

        if gtf_file.exists():
            self.gtf_file = str(gtf_file)
            logger.info(f"File already exists: {gtf_file}")
            return str(gtf_file)

        try:
            with gzip.open(self.gz_file, 'rb') as f_in:
                with open(gtf_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info(f"Successfully unzipped file to: {gtf_file}")
        except OSError as e:
            logger.info(f"Error occurred while unzipping the file: {e}")
            raise

        self.gtf_file = str(gtf_file)

        return str(gtf_file)
    
    def get_all_genes(self) -> str:

        if not hasattr(self, 'gtf_file') or not os.path.isfile(self.gtf_file):
            raise FileNotFoundError("Reference file not found")
        
        if os.path.isfile(self.gtf_file[:-4]+"-all_genes.gtf.gz"):

            self.all_genes_path = self.gtf_file[:-4] + "-all_genes.gtf.gz"
            logger.info(f"File already exists: {self.all_genes_path}")

            return

        gtf = read_gtf(self.gtf_file, usecols=["feature", "gene_biotype", "gene_id", "gene_name"], result_type='pandas')

        gene_list = gtf.loc[gtf["feature"]=="gene", "gene_id"].values

        gtf_raw = pd.read_csv(
            self.gtf_file, 
            sep="\t", 
            header=None, 
            comment="#", 
            dtype="string"
        )

        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[ gtf_raw["_gene_id"].isin(gene_list) ,:]
        gtf_raw = gtf_raw.drop("_gene_id",axis=1)

        all_genes_path = self.gtf_file[:-4]+"-all_genes.gtf.gz"

        gtf_raw.to_csv(all_genes_path, header=None, index=None, sep="\t")
        logger.info(f"Saved all genes to: {all_genes_path}")

        self.all_genes_path = all_genes_path

        return all_genes_path
    
class Ensembl38Fetcher(ReferenceDataFetcher):

    def __init__(self, destination_folder = None):
        
        super().__init__(
            base_url = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/",
            build ='38',
            source = 'ensembl',
            destination_folder = destination_folder
        )

    def get_latest_release(self) -> None:
        """
        Fetches the latest GTF file dynamically from the specified base URL.

        This method sends a GET request to the base URL, parses the HTML response to find the latest GTF file link, and sets the `latest_url` attribute to the full URL of the latest GTF file.
        
        Raises:
        ------
        FileNotFoundError: If no GTF file is found in the HTML response.
        
        Returns:
        --------
        None
        """

        # Get the latest file dynamically
        response = requests.get(self.base_url)

        if response.status_code != 200:
            raise Exception(f"Failed to access {self.base_url}")
        
        soup = BeautifulSoup(response.text, "html.parser")

        # Find the latest GTF file
        latest_gtf = None
        for link in soup.find_all("a"):
            href = link.get("href")
            if href and "Homo_sapiens.GRCh38" in href and href.endswith(".chr.gtf.gz"):
                latest_gtf = href
                break  # Assuming the first match is the latest
            
        if latest_gtf:
            latest_url = self.base_url + latest_gtf
            logger.info(f"Latest GTF file: {latest_gtf}")
            logger.info(f"Download URL: {latest_url}")
            self.latest_url = latest_url
        else:
            raise FileNotFoundError("GTF file not found")
        
        pass

class Ensembl37Fetcher(ReferenceDataFetcher):

    def __init__(self, destination_folder = None):
        
        super().__init__(
            base_url = 'https://ftp.ensembl.org/pub/grch37/', 
            build ='37', 
            source = 'ensembl', 
            destination_folder = destination_folder
        )

    def get_latest_release(self) -> None:

        response = requests.get(self.base_url)

        if response.status_code != 200:
            raise Exception(f"Failed to access {self.base_url}")

        # Parse the HTML content
        soup = BeautifulSoup(response.text, "html.parser")

        # Find all folder names matching 'release-*'
        releases = []
        
        for link in soup.find_all("a"):
            
            href = link.get("href")
            match = re.match(r"release-(\d+)", href)
            
            if match:
                releases.append(int(match.group(1)))  # Extract the release number as integer

        if not releases:
            raise Exception("No release folders found.")

        latest_release = max(releases)  # Get the highest release number
        latest_folder = self.base_url + f"release-{latest_release}/" + 'gtf/homo_sapiens/'

        response = requests.get(latest_folder)

        if response.status_code != 200:
            raise Exception(f"Failed to access {latest_folder}")
        
        # Parse the HTML content
        soup = BeautifulSoup(response.text, "html.parser")

        latest_gtf = None

        for link in soup.find_all("a"):
            href = link.get("href")
            if href and "Homo_sapiens.GRCh37" in href and href.endswith(".chr.gtf.gz"):
                latest_gtf = href
                break  # Assuming the first match is the latest
            
        if latest_gtf:
            latest_url = latest_folder + latest_gtf
            logger.info(f"Latest GTF file: {latest_gtf}")
            logger.info(f"Download URL: {latest_url}")
            self.latest_url = latest_url
        else:
            raise FileNotFoundError("GTF file not found")
        
        pass

class RefSeqFetcher(ReferenceDataFetcher):

    def __init__(self, build: str, destination_folder: str = None):

        super().__init__(
            base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/", 
            build = build, 
            source = 'refseq', 
            destination_folder = destination_folder
        )

    def get_latest_release(self) -> None:
        """
        Fetches the latest GTF file dynamically from the specified base URL.

        This method sends a GET request to the base URL, parses the HTML response to find the latest GTF file link, and sets the `latest_url` attribute to the full URL of the latest GTF file.
        
        Raises:
        ------
        FileNotFoundError: If no GTF file is found in the HTML response.
        
        Returns:
        --------
        None
        """

        # Get the latest file dynamically
        response = requests.get(self.base_url)

        if response.status_code != 200:
            raise Exception(f"Failed to access {self.base_url}")
        
        soup = BeautifulSoup(response.text, "html.parser")

        if self.build == "38":
            version_name = "GRCh38"
        elif self.build == "37":
            version_name = "GRCh37"

        # Find all folder names matching 'release-*'
        latest_release = ''
        latest_release_num = 0

        for link in soup.find_all("a"):
            
            href = link.get("href")
            if version_name in href:
                version = href.split('.')[-1][1:-1]
                
                if version.isdigit():
                    version_num = int(version)
                    if version_num > latest_release_num:
                        latest_release_num = version_num
                        latest_release = href

        if len(latest_release)==0:
            raise Exception("No release folders found.")

        latest_folder = self.base_url + latest_release

        response = requests.get(latest_folder)

        if response.status_code != 200:
            raise Exception(f"Failed to access {latest_folder}")
        
        # Parse the HTML content
        soup = BeautifulSoup(response.text, "html.parser")

        latest_gtf = None

        for link in soup.find_all("a"):
            href = link.get("href")
            if href and version_name in href and href.endswith(".gtf.gz"):
                latest_gtf = href
                break  # Assuming the first match is the latest
            
        if latest_gtf:
            latest_url = latest_folder + latest_gtf
            logger.info(f"Latest GTF file: {latest_gtf}")
            logger.info(f"Download URL: {latest_url}")
            self.latest_url = latest_url
        else:
            raise FileNotFoundError("GTF file not found")
        
        return
