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
from typing import Optional

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
        
class ReferenceDataFetcher:
    """A class for fetching, downloading, and processing reference genome data.
    
    This class provides a framework for retrieving genomic reference data from various
    sources. It handles downloading compressed files, unzipping them, and extracting
    gene information from GTF files.

    Attributes
    ----------
    build : str
        The genome build (e.g., 'hg38', 'GRCh38').
    source : str
        The data source (e.g., 'ensembl', 'ucsc').
    base_url : str
        The base URL to fetch data from.
    destination_folder : Optional[str]
        The directory to save downloaded files. If None, defaults to project_root/data/{source}_latest.
    latest_url : Optional[str]
        The URL of the latest release after calling get_latest_release().
    gz_file : Optional[str]
        Path to the downloaded compressed file.
    gtf_file : Optional[str]
        Path to the uncompressed GTF file.


    Notes
    -----
    This is an abstract base class that requires subclasses to implement the
    get_latest_release() method for specific data sources.
    """

    def __init__(self, base_url: str, build: str, source: str, destination_folder: Optional[str] = None) -> None:

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
        """Downloads the latest file from `self.latest_url` to `self.destination_folder`.

        Raises
        ------
        AttributeError 
            If `self.latest_url` is not set.
        requests.exceptions.RequestException 
            If the HTTP request fails.
        """

        if not self.latest_url:
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

        if not self.latest_url:
            raise AttributeError("`self.latest_url` is not set. Call `get_latest_release` first.")

        self.destination_folder = self.get_destination_folder()

        if not hasattr(self, 'gz_file') or self.gz_file is None or not Path(self.gz_file).is_file():
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
        """Extract all genes from the GTF file and save them to a new compressed file.
        
        This method reads the GTF file specified in self.gtf_file, filters for gene features,
        and creates a new GTF file containing only the gene entries. If the output file 
        already exists, it will return the path without reprocessing.
        
        Returns
        -------
        str: 
            Path to the compressed GTF file containing all genes
        
        Raises
        ------
        FileNotFoundError 
            If the reference GTF file (self.gtf_file) is not found
        TypeError 
            If read_gtf does not return a pandas DataFrame
        
        Note
        ----
            The output file will be named based on the input GTF file with "-all_genes.gtf.gz" suffix
        """

        if not hasattr(self, 'gtf_file') or self.gtf_file is None or not os.path.isfile(self.gtf_file):
            raise FileNotFoundError("Reference file not found")
        
        if os.path.isfile(self.gtf_file[:-4]+"-all_genes.gtf.gz"):

            self.all_genes_path = self.gtf_file[:-4] + "-all_genes.gtf.gz"
            logger.info(f"File already exists: {self.all_genes_path}")

            return self.all_genes_path

        gtf = read_gtf(
            self.gtf_file, 
            usecols    =["feature", "gene_biotype", "gene_id", "gene_name"], 
            result_type='pandas'
        )
        if not isinstance(gtf, pd.DataFrame):
            raise TypeError("read_gtf did not return a pandas DataFrame")

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

        gtf_raw.to_csv(all_genes_path, header=False, index=False, sep="\t")
        logger.info(f"Saved all genes to: {all_genes_path}")

        self.all_genes_path = all_genes_path

        return all_genes_path
    
class Ensembl38Fetcher(ReferenceDataFetcher):
    """A class for fetching human genome reference data from Ensembl based on GRCh38 build.
    
    This class extends ReferenceDataFetcher to specifically handle Ensembl's human
    genome data with build 38. It provides functionality to find and retrieve the
    latest GTF file from Ensembl's FTP server.
    
    Attributes
    ----------
    base_url : str
        Base URL for Ensembl FTP server where GTF files are stored
    build : str
        Genome build version ('38')
    source : str
        Data source ('ensembl')
    destination_folder : str 
        Local folder to store downloaded files
    latest_url : str
        URL of the latest GTF file after calling get_latest_release()
    
    Raises
    ------
    Exception
        If the Ensembl FTP server cannot be accessed
    FileNotFoundError
        If no matching GTF file is found
    """


    def __init__(self, destination_folder = None):
        
        super().__init__(
            base_url = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/",
            build ='38',
            source = 'ensembl',
            destination_folder = destination_folder
        )

    def get_latest_release(self) -> None:
        """Retrieves the URL of the latest GTF file for human genome (GRCh38) from the base URL.
        
        This method scrapes the base URL to find the most recent Homo_sapiens GRCh38 GTF file
        available for download. Upon finding the file, it constructs the complete URL and
        stores it in the instance variable `latest_url`.
        
        Returns
        -------
            None

        Raises
        ------
        Exception
            If the base URL cannot be accessed (non-200 response)
        FileNotFoundError
            If no GTF file matching the criteria is found
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
    """A class for fetching reference genome data from Ensembl's GRCh37 (hg19) repository.
    
    This class specializes the ReferenceDataFetcher to work specifically with 
    Ensembl's GRCh37 human genome build. It provides functionality to automatically
    detect and download the latest available GTF file for Homo sapiens from the
    Ensembl GRCh37 archive.
    The fetcher connects to Ensembl's FTP server, identifies the most recent release
    available for GRCh37, and locates the chromosome GTF file for human genome data.

    Attributes
    ----------
    base_url : str
        The base URL for Ensembl's GRCh37 repository
    build : str
        The genome build identifier ('37')
    source : str
        The data source identifier ('ensembl')
    latest_url : str
        The complete URL to the latest GTF file, populated after calling get_latest_release()
    """

    def __init__(self, destination_folder = None):
        """
        Initialize a reference genome downloader for Ensembl GRCh37.
        This constructor configures the downloader to retrieve data from Ensembl's GRCh37 
        repository.
        
        Parameters
        ----------
        destination_folder : str, optional
            The folder where downloaded files will be stored. If None, a default 
            location will be used based on the parent class implementation.
        """
        
        super().__init__(
            base_url = 'https://ftp.ensembl.org/pub/grch37/', 
            build ='37', 
            source = 'ensembl', 
            destination_folder = destination_folder
        )

    def get_latest_release(self) -> None:
        """Fetches the URL of the latest GTF file for Homo sapiens GRCh37 from Ensembl.
        
        This method:
        1. Connects to the base URL and identifies all available release folders
        2. Determines the latest release by finding the highest release number
        3. Navigates to the GTF directory for that release
        4. Locates the Homo sapiens GRCh37 chromosome GTF file
        5. Stores the complete download URL in self.latest_url
        
        Raises
        ------
        Exception
            If the base URL cannot be accessed
        Exception
            If no release folders are found
        Exception
            If the latest release folder cannot be accessed
        FileNotFoundError
            If the GTF file is not found in the latest release
        
        Returns
        -------
        None
        """

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
    """A class for fetching and downloading reference genome data from NCBI's RefSeq repository.
    
    This class extends ReferenceDataFetcher to specifically handle downloading human
    genome reference files from the RefSeq database. It supports different genome 
    builds (e.g., 'GRCh37', 'GRCh38') and automatically identifies the latest version
    available for the specified build.
    The class handles navigating the NCBI FTP directory structure, finding the
    appropriate GTF files for the requested genome build, and managing the
    download process.

    Attributes
    ----------
    base_url : str
        The base URL for the NCBI RefSeq FTP server directory.
    build : str
        The genome build version ('37' for GRCh37, '38' for GRCh38).
    source : str
        The source of the reference data (set to 'refseq').
    latest_url : str
        URL to the latest GTF file, set after calling get_latest_release().
    """

    def __init__(self, build: str, destination_folder: Optional[str] = None):

        super().__init__(
            base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/", 
            build = build, 
            source = 'refseq', 
            destination_folder = destination_folder
        )

    def get_latest_release(self) -> None:
        """Fetches the latest GTF file dynamically from the specified base URL.

        This method sends a GET request to the base URL, parses the HTML response 
        to find the latest GTF file link, and sets the `latest_url` attribute 
        to the full URL of the latest GTF file.
        
        Raises
        ------
        FileNotFoundError
            If no GTF file is found in the HTML response.
        
        Returns
        -------
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
        else:
            raise ValueError("Unsupported build version. Only 'GRCh37' and 'GRCh38' are supported.")

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

class AssemblyReferenceFetcher():
    """A class for fetching and preparing genomic reference files from online repositories.
    
    This class handles the process of:
    1. Finding the appropriate reference file URL based on build parameters
    2. Downloading the reference file 
    3. Unzipping compressed reference files if necessary
    
    Parameters
    ----------
    base_url : str
        The base URL where reference files are hosted
    build : str
        The genome build identifier (e.g., 'GRCh38', 'hg19')
    extension : str
        File extension to look for (e.g., '.gtf.gz', '.fa.gz')
    destination_folder : Optional[str], default=None
        Path where files should be downloaded. If None, uses project_root/data/assembly_references
    avoid_substring : str, default='extra'
        Substring to avoid when selecting reference files
    
    Attributes
    ----------
    reference_url : str or None
        URL of the identified reference file
    reference_file : str or None
        Filename of the identified reference file
    file_path : Path or None
        Local path to the downloaded reference file
    
    Raises
    ------
    Exception
        If the base URL cannot be accessed
    FileNotFoundError
        If no matching reference file is found
    AttributeError
        If methods are called out of sequence
    ValueError
        If required attributes are None when needed
    """


    def __init__(self, base_url: str, build: str, extension: str, destination_folder: Optional[str] = None, avoid_substring: str = 'extra') -> None:

        self.base_url = base_url
        self.build = build
        self.destination_folder = destination_folder
        self.extension = extension
        self.avoid_substring = avoid_substring

        self.file_path = None
        self.reference_url = None
        self.reference_file = None
        
        pass

    def get_reference_url(self) -> str:
        """Retrieves the URL for the reference file from the base URL.
        
        This method performs an HTTP GET request to the base URL, parses the HTML content,
        and searches for links matching specific criteria:
        - Contains the build version string
        - Ends with the specified extension
        - Does not contain the specified substring to avoid
        The first matching link is considered the reference file.
        
        Returns
        -------
            str: The complete URL to the reference file
        
        Raises
        ------
        Exception
            If the base URL cannot be accessed
        FileNotFoundError
            If no matching reference file is found
        
        Notes
        -----
            - Sets self.reference_file to the name of the found file
            - Sets self.reference_url to the complete URL
            - Logs information about the found file and URL
        """


        response = requests.get(self.base_url)

        if response.status_code != 200:
            raise Exception(f"Failed to access {self.base_url}")
        
        soup = BeautifulSoup(response.text, "html.parser")

        reference_file = None
        for link in soup.find_all("a"):
            href = link.get("href")
            if href and f"{self.build}" in href and href.endswith(f"{self.extension}") and self.avoid_substring not in href:
                reference_file = href
                break  # Assuming the first match is the latest

        self.reference_file = reference_file

        if reference_file:
            reference_url = self.base_url + reference_file
            logger.info(f"Latest GTF file: {reference_file}")
            logger.info(f"Download URL: {reference_url}")
            self.reference_url = reference_url
        else:
            raise FileNotFoundError("Reference file not found")
        
        return str(reference_url)

    def download_reference_file(self) -> str:
        """
        Downloads a reference file from the specified URL to the destination folder.
        
        This method first checks whether the reference file already exists locally.
        If not found, it also looks for an alternative version with a '.fa' extension.
        If neither is present, it downloads the file from the given URL.
        
        Raises
        ------
        AttributeError
            If `self.reference_url` or `self.reference_file` are not set.
        ValueError
            If `self.reference_url` or `self.reference_file` are set to None.
        
        Returns
        -------
        str
            The path to the downloaded or existing reference file.
        
        Note
        ----
        `self.reference_url` and `self.reference_file` must be set by calling `get_reference_url()` before using this method.
        """

        if not getattr(self, 'reference_url', None):
            raise AttributeError("`self.reference_url` is not set. Call `get_reference_url` first.")
        
        if not getattr(self, 'reference_file', None):
            raise AttributeError("`self.reference_file` is not set. Call `get_reference_url` first.")
        
        if self.reference_file is None:
            raise ValueError("reference_file is None. Cannot construct file path.")
        if self.reference_url is None:
            raise ValueError("reference_url is None. Cannot download file.")

        self.destination_folder = self.get_destination_folder()

        file_path = self.destination_folder / Path(self.reference_file).name

        if file_path.exists():
            logger.info(f"File already exists: {file_path}")
            self.file_path = file_path
            return str(file_path)
        
        fa_file = file_path.with_suffix('.fa')
        
        if fa_file.exists():
            logger.info(f"File already exists: {fa_file}")
            self.file_path = file_path
            return str(fa_file)
        
        self._download_file(self.reference_url, file_path)
        self.file_path = file_path

        return str(file_path)
    
    def unzip_reference_file(self) -> str:
        """Unzips a reference genome file (typically .fa.gz to .fa) and returns the path to the unzipped file.
        
        This method checks if the file is already unzipped, and if not, unzips it using gzip.
        After successful unzipping, the original compressed file is deleted.
        
        Returns
        -------
        str
            Path to the unzipped reference file (.fa)
        
        Raises
        ------
        AttributeError
            If self.reference_file is not set (get_reference_url should be called first)
        AttributeError
            If self.file_path is not set or None (download_reference_file should be called first)
        OSError
            If an error occurs during the unzipping process
        """

        if not getattr(self, 'reference_file', None):
            raise AttributeError("`self.reference_file` is not set. Call `get_reference_url` first.")
        
        if not getattr(self, 'file_path', None) or self.file_path is None:
            raise AttributeError("`self.file_path` is not set. Call `download_reference_file` first.")
        
        if self.file_path.suffix == '.fa':
            logger.info(f"File already unzipped: {self.file_path}")
            return str(self.file_path)
        
        fa_file = self.get_destination_folder() / Path(self.file_path.stem)
        if fa_file.exists():
            logger.info(f"File already exists: {fa_file}")
            return str(fa_file)
        
        try:
            with gzip.open(self.file_path, 'rb') as f_in:
                with open(fa_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info(f"Successfully unzipped file to: {fa_file}")
        except OSError as e:
            logger.info(f"Error occurred while unzipping the file: {e}")
            raise

        self.file_path.unlink()
        self.file_path = fa_file

        return str(fa_file)
    
    def get_destination_folder(self) -> Path:
        """Determines and creates (if necessary) the destination folder for reference files.

        If a destination folder was provided during initialization, it uses that path.
        Otherwise, it defaults to a 'data/assembly_references' directory in the project root.

        Returns
        -------
        Path 
            The path to the destination folder where reference files will be stored.
        """

        if self.destination_folder:
            destination = Path(self.destination_folder)
        else:
            # Determine project root and default `data` directory
            project_root = Path(__file__).resolve().parent.parent
            destination = project_root / "data" / "assembly_references"

        destination.mkdir(parents=True, exist_ok=True)

        return destination
    
    def _download_file(self, url: str, file_path: Path) -> None:
        """Download a file from a given URL and save it to the specified path.
        
        Parameters
        ----------
        url : str
            The URL to download the file from.
        file_path : Path
            The path where the downloaded file will be saved.
            
        Returns
        -------
        None
            
        Raises
        ------
        HTTPError
            If the HTTP request returns an unsuccessful status code.
        """

        with requests.get(url, stream=True) as response:
            response.raise_for_status()

            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

        logger.info(f"Downloaded file to: {file_path}")

class FetcherLDRegions:
    """A class for fetching high Linkage Disequilibrium (LD) regions files.
    
    This class handles downloading or creating files containing genomic regions
    with high LD for different genome builds (37 or 38). These regions are often
    excluded in GWAS analyses to avoid confounding effects.
    
    Parameters
    ----------
    destination : Path, optional
        Directory path where the LD regions files will be stored.
        Default is "../data/ld_regions_files" relative to the module location.
    built : str, default '38'
        Genome build version. Must be either '37' or '38'.
    
    Attributes
    ----------
    destination : Path
        Directory where LD regions files are stored.
    built : str
        Genome build version being used.
    ld_regions : Path or None
        Path to the LD regions file once retrieved, None otherwise.
    """

    def __init__(self, destination: Optional[Path] = None, built: str = '38'):

        if not destination:
            destination = Path(__file__).resolve().parent.parent / "data" / "ld_regions_files"

        self.destination = destination
        self.built = built

        self.ld_regions = None
        
        pass

    def get_ld_regions(self)-> Path:
        """Downloads or creates high LD regions file based on genome build version.
        
        This method handles the retrieval of high Linkage Disequilibrium (LD) regions for
        different genome builds (37 or 38). For build 37, it downloads the regions from a
        GitHub repository. For build 38, it creates the file from predefined coordinates.
        
        Returns
        -------
        Path
            Path to the created/downloaded LD regions file. Returns empty Path if download fails for build 37.
        
        Notes
        -----
        - For build 37: Downloads from genepi-freiburg/gwas repository
        - For build 38: Creates file from hardcoded coordinates from GWAS-pipeline
        - Files are named as 'high-LD-regions_GRCh{build}.txt'
        - Creates destination directory if it doesn't exist
        """

        self.destination.mkdir(parents=True, exist_ok=True)

        out_dir = self.destination

        if self.built == '37':
            url_ld_regions = r"https://raw.githubusercontent.com/genepi-freiburg/gwas/refs/heads/master/single-pca/high-LD-regions.txt"
        
            ld = requests.get(url_ld_regions)


            if ld.status_code == 200:
                with open((out_dir / f"high-LD-regions_GRCh{self.built}.txt"), "wb") as f:
                    f.write(ld.content)
                logger.info(f"LD regions file for built {self.built} downloaded successfully to {out_dir}")

                self.ld_regions = out_dir / f"high-LD-regions_GRCh{self.built}.txt"
                return out_dir / f"high-LD-regions_GRCh{self.built}.txt"
            else:
                logger.info(f"Failed to download .bim file: {ld.status_code}")

                return Path()

        elif self.built == '38':
            # extracted from
            # https://github.com/neurogenetics/GWAS-pipeline
            data = [
                (1, 47534328, 51534328, "r1"),
                (2, 133742429, 137242430, "r2"),
                (2, 182135273, 189135274, "r3"),
                (3, 47458510, 49962567, "r4"),
                (3, 83450849, 86950850, "r5"),
                (5, 98664296, 101164296, "r6"),
                (5, 129664307, 132664308, "r7"),
                (5, 136164311, 139164311, "r8"),
                (6, 24999772, 35032223, "r9"),
                (6, 139678863, 142178863, "r10"),
                (8, 7142478, 13142491, "r11"),
                (8, 110987771, 113987771, "r12"),
                (11, 87789108, 90766832, "r13"),
                (12, 109062195, 111562196, "r14"),
                (20, 33412194, 35912078, "r15")
            ]

            with open(out_dir / f'high-LD-regions_GRCH{self.built}.txt', 'w') as file:
                for line in data:
                    file.write(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\n")
            self.ld_regions = out_dir / f"high-LD-regions_GRCH{self.built}.txt"
            return out_dir / f'high-LD-regions_GRCH{self.built}.txt'
        else:
            raise ValueError("Unsupported genome build. Only '37' and '38' are supported.")
