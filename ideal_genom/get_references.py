import requests
import os
import gzip
import shutil
import re

import pandas as pd

from bs4 import BeautifulSoup
from gtfparse import read_gtf

class Ensembl38:

    def __init__(self, base_url: str = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/") -> None:
        
        # URL of the GTF directory
        self.base_url = base_url

        pass

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
            print(f"Latest GTF file: {latest_gtf}")
            print(f"Download URL: {latest_url}")
            self.latest_url = latest_url
        else:
            raise FileNotFoundError("GTF file not found")
        
        pass
        
    def download_latest(self, destination_folder: str = None) -> None:
        
        """
        Downloads the latest file from the specified URL to the given destination folder.

        Note: `self.latest_url` must be set by calling `get_latest_release` before calling this method.

        If no destination folder is provided, the file will be saved in the `data/ensembl_latest`
        directory within the project root.

        Parameters:
        -----------
        destination_folder (str, optional): 
            The folder where the downloaded file will be saved. If not provided, defaults to `data/ensembl_latest` within the project root.

        Returns:
        --------
        None

        Raises:
        -------
            requests.exceptions.RequestException: If there is an issue with the HTTP request.

        Side Effects:
        -------------
            - Creates the destination folder if it does not exist.
            - Downloads the file from `self.latest_url` and saves it to the destination folder.
            - Sets `self.gz_file` to the path of the downloaded file.
        """

        if not getattr(self, 'latest_url', None):
            raise AttributeError("`self.latest_url` is not set. Call `get_latest_release` before calling this method.")

        if destination_folder is None:

            library_path = os.path.dirname(os.path.abspath(__file__))

            # Go up one level to reach the project root
            project_root = os.path.dirname(library_path)

            # Define the path to the `data` directory
            data_dir = os.path.join(project_root, "data", "ensembl_latest")

            # Ensure the `data` directory exists
            os.makedirs(data_dir, exist_ok=True)

            # Set the destination folder to the library path if not provided
            destination_folder = data_dir

        # Ensure the destination folder exists
        os.makedirs(destination_folder, exist_ok=True)

        # Define the path to save the downloaded file
        file_name = os.path.join(destination_folder, os.path.basename(self.latest_url))

        if os.path.exists(file_name):
            self.gz_file = file_name
            print(f"File already exists: {file_name}")
            return

        # Download the file
        with requests.get(self.latest_url, stream=True) as response:
            response.raise_for_status()
            with open(file_name, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
            self.gz_file = file_name

        print(f"Downloaded file to: {file_name}")

        return

    def unzip_latest(self, origin_folder: str = None, destination_folder: str = None) -> None:

        if not getattr(self, 'latest_url', None):
            raise AttributeError("`self.latest_url` is not set. Call `get_latest_release` before calling this method.")
            
        if origin_folder is None:

            library_path = os.path.dirname(os.path.abspath(__file__))

            # Go up one level to reach the project root
            project_root = os.path.dirname(library_path)

            # Define the path to the `data` directory
            data_dir = os.path.join(project_root, "data", "ensembl_latest")

            # Set the destination folder to the library path if not provided
            origin_folder = data_dir

        if destination_folder is None:

            destination_folder = origin_folder

        if not hasattr(self, 'gz_file') or not os.path.isfile(self.gz_file):
            raise FileNotFoundError("Reference file not found")
        
        gtf_file = os.path.join(destination_folder, os.path.basename(self.latest_url)[:-3])

        try:
            with gzip.open(self.gz_file, 'rb') as f_in:
                with open(gtf_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                self.gtf_file = gtf_file
        except OSError as e:
            print(f"Error occurred while copying file: {e}")
            raise

        return
    
    def get_all_genes(self) -> None:

        if not hasattr(self, 'gtf_file') or not os.path.isfile(self.gtf_file):
            raise FileNotFoundError("Reference file not found")
        
        if os.path.isfile(self.gtf_file[:-5]+"-all_genes.gtf.gz"):
            self.all_genes_path = self.gtf_file[:-5]+"-all_genes.gtf.gz"
            print(f"File already exists: {self.all_genes_path}")
            return

        gtf = read_gtf(self.gtf_file, usecols=["feature","gene_biotype","gene_id","gene_name"], result_type='pandas')

        gene_list = gtf.loc[gtf["feature"]=="gene","gene_id"].values

        gtf_raw = pd.read_csv(self.gtf_file, sep="\t", header=None, comment="#", dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[ gtf_raw["_gene_id"].isin(gene_list) ,:]
        gtf_raw = gtf_raw.drop("_gene_id",axis=1)

        all_genes_path = self.gtf_file[:-5]+"-all_genes.gtf.gz"

        gtf_raw.to_csv(all_genes_path, header=None, index=None, sep="\t")

        self.all_genes_path = all_genes_path

        return
    
class Ensembl37:

    def __init__(self, base_url: str = 'https://ftp.ensembl.org/pub/grch37/') -> None:

        self.base_url = base_url

        pass

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
            print(f"Latest GTF file: {latest_gtf}")
            print(f"Download URL: {latest_url}")
            self.latest_url = latest_url
        else:
            raise FileNotFoundError("GTF file not found")
        
        pass

    def download_latest(self, destination_folder: str = None) -> None:
        
        """
        Downloads the latest file from the specified URL to the given destination folder.

        Note: `self.latest_url` must be set by calling `get_latest_release` before calling this method.

        If no destination folder is provided, the file will be saved in the `data/ensembl_latest`
        directory within the project root.

        Parameters:
        -----------
        destination_folder (str, optional): 
            The folder where the downloaded file will be saved. If not provided, defaults to `data/ensembl_latest` within the project root.

        Returns:
        --------
        None

        Raises:
        -------
            requests.exceptions.RequestException: If there is an issue with the HTTP request.

        Side Effects:
        -------------
            - Creates the destination folder if it does not exist.
            - Downloads the file from `self.latest_url` and saves it to the destination folder.
            - Sets `self.gz_file` to the path of the downloaded file.
        """

        if not getattr(self, 'latest_url', None):
            raise AttributeError("`self.latest_url` is not set. Call `get_latest_release` before calling this method.")

        if destination_folder is None:

            library_path = os.path.dirname(os.path.abspath(__file__))

            # Go up one level to reach the project root
            project_root = os.path.dirname(library_path)

            # Define the path to the `data` directory
            data_dir = os.path.join(project_root, "data", "ensembl_latest")

            # Ensure the `data` directory exists
            os.makedirs(data_dir, exist_ok=True)

            # Set the destination folder to the library path if not provided
            destination_folder = data_dir

        # Ensure the destination folder exists
        os.makedirs(destination_folder, exist_ok=True)

        # Define the path to save the downloaded file
        file_name = os.path.join(destination_folder, os.path.basename(self.latest_url))

        if os.path.exists(file_name):
            self.gz_file = file_name
            print(f"File already exists: {file_name}")
            return

        # Download the file
        with requests.get(self.latest_url, stream=True) as response:
            response.raise_for_status()
            with open(file_name, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
            self.gz_file = file_name

        print(f"Downloaded file to: {file_name}")

        return
    
    def unzip_latest(self, origin_folder: str = None, destination_folder: str = None) -> None:

        if not getattr(self, 'latest_url', None):
            raise AttributeError("`self.latest_url` is not set. Call `get_latest_release` before calling this method.")
            
        if origin_folder is None:

            library_path = os.path.dirname(os.path.abspath(__file__))

            # Go up one level to reach the project root
            project_root = os.path.dirname(library_path)

            # Define the path to the `data` directory
            data_dir = os.path.join(project_root, "data", "ensembl_latest")

            # Set the destination folder to the library path if not provided
            origin_folder = data_dir

        if destination_folder is None:

            destination_folder = origin_folder

        if not hasattr(self, 'gz_file') or not os.path.isfile(self.gz_file):
            raise FileNotFoundError("Reference file not found")
        
        gtf_file = os.path.join(destination_folder, os.path.basename(self.latest_url)[:-2])

        try:
            with gzip.open(self.gz_file, 'rb') as f_in:
                with open(gtf_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                self.gtf_file = gtf_file
        except OSError as e:
            print(f"Error occurred while copying file: {e}")
            raise

        return
    
    def get_all_genes(self) -> None:

        if not hasattr(self, 'gtf_file') or not os.path.isfile(self.gtf_file):
            raise FileNotFoundError("Reference file not found")
        
        if os.path.isfile(self.gtf_file[:-5]+"protein_coding.gtf.gz"):
            self.protein_coding_path = self.gtf_file[:-3]+"protein_coding.gtf"
            print(f"File already exists: {self.protein_coding_path}")
            return

        gtf = read_gtf(self.gtf_file, usecols=["feature","gene_biotype","gene_id","gene_name"])

        gene_list = gene_list = gtf.loc[gtf["feature"]=="gene","gene_id"].values

        gtf_raw = pd.read_csv(self.gtf_file, sep="\t", header=None, comment="#", dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[ gtf_raw["_gene_id"].isin(gene_list) ,:]
        gtf_raw = gtf_raw.drop("_gene_id",axis=1)

        protein_coding_path = self.gtf_file[:-6]+"protein_coding.gtf.gz"

        gtf_raw.to_csv(protein_coding_path, header=None, index=None, sep="\t")

        self.protein_coding_path = protein_coding_path

        return
