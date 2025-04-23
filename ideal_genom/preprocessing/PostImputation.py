import time
import logging
import os
import threading
import zipfile
import subprocess

from pathlib import Path
from typing import Callable, List, Dict, Any, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

from ideal_genom.get_references import AssemblyReferenceFetcher

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class ParallelTaskRunner:
    """
    A base class for running parallel tasks on files.
    
    This class provides the basic infrastructure for parallel processing of files
    using ThreadPoolExecutor. It handles file collection and parallel task execution
    while providing progress monitoring and logging.
    
    Attributes:
    -----------
        input_path (Path): Directory path where input files are located.
        output_path (Path): Directory path where output files will be saved.
        max_workers (int): Maximum number of worker threads to use. Defaults to min(8, CPU count).
        files (List[Path]): List of files to be processed.
    
    Example:
    --------
        class MyTaskRunner(ParallelTaskRunner):
                def task_fn(file):
                    # Process file
                self._file_collector("*.txt")
                self._run_task(task_fn, {})
   
    Raises:
    -------
        TypeError: If input_path or output_path are not Path objects.
        FileNotFoundError: If input_path or output_path don't exist.
        NotADirectoryError: If input_path or output_path are not directories.
    """
    
    def __init__(self, input_path: Path, output_path: Path, max_workers: Optional[int] = None) -> None:
        """
        Initialize a ParallelTask processor.

        This constructor sets up paths for processing multiple VCF in parallel and validates
        that both input and output directories exist and are properly formatted.

        Parameters:
        -----------
        input_path (Path): 
            Path to the directory containing input files to process.
        output_path (Path): 
            Path to the directory where processed files will be saved.
        max_workers (Optional[int], optional): 
            Maximum number of concurrent workers for parallel processing. If None, defaults to min(8, os.cpu_count()).

        Raises:
        -------
        TypeError: If input_path or output_path is not a Path object.
        FileNotFoundError: If input_path or output_path does not exist.
        NotADirectoryError: If input_path or output_path is not a directory.
        """

        if not isinstance(input_path, Path):
            raise TypeError(f"input_path should be of type Path, got {type(input_path)}")
        if not isinstance(output_path, Path):
            raise TypeError(f"output_path should be of type Path, got {type(output_path)}")
        if not input_path.exists():
            raise FileNotFoundError(f"Input path {input_path} does not exist.")
        if not output_path.exists():
            raise FileNotFoundError(f"Output path {output_path} does not exist.")
        if not input_path.is_dir():
            raise NotADirectoryError(f"Output path {input_path} is not a directory.")
        if not output_path.is_dir():
            raise NotADirectoryError(f"Output path {output_path} is not a directory.")

        self.input_path = input_path
        self.output_path = output_path

        self.max_workers = max_workers or min(8, os.cpu_count())

        self.files = []

        pass

    def execute_task(self) -> None:
        """
        Execute the specific post-imputation processing task.

        This abstract method should be implemented by all subclasses to perform
        their specific post-imputation processing operations. Implementations
        should handle the execution logic for the particular task the subclass
        is designed to perform.

        Returns:
        --------
        None

        Raises:
        -------
        NotImplementedError: If the subclass does not implement this method.
        """

        raise NotImplementedError("Subclasses must implement this method.")

    def _file_collector(self, filename_pattern: str) -> List[Path]:
        """
        Collect files matching a given pattern from the input directory.
        This method finds all files matching the specified glob pattern within the
        input directory, sorts them, and stores the resulting list as an instance attribute.
        
        Parameters
        ----------
        filename_pattern : str
            A glob pattern string to match files (e.g., "*.vcf.gz").
        
        Returns
        -------
        List[Path]
            A sorted list of Path objects for the files matching the pattern.
        
        Raises
        ------
        TypeError
            If filename_pattern is not a string.
        FileNotFoundError
            If no files match the given pattern in the input directory.
        
        Notes
        -----
        The matched files are also stored in the instance attribute `files`.
        """

        if not isinstance(filename_pattern, str):
            raise TypeError(f"filename_pattern should be of type str, got {type(filename_pattern)}")

        files = list(self.input_path.glob(filename_pattern))
        files.sort()

        self.files = files
        logger.info(f"Found {len(files)} files matching pattern {filename_pattern} in {self.input_path}")
        
        if not files:
            raise FileNotFoundError(f"No files found matching pattern {filename_pattern} in {self.input_path}")

        return files
    
    def _run_task(self, task_fn: Callable, task_args: Dict[str, Any], desc: str = "Running tasks") -> None:
        """
        Execute a task function across all files using parallel processing with ThreadPoolExecutor.
        This method applies the given task function to each file in self.files concurrently,
        managing thread allocation, progress tracking, and error handling.
        
        Parameters
        ----------
        task_fn : Callable
            The function to execute for each file. First argument should accept a file,
            and it should accept **kwargs for additional arguments.
        task_args : Dict[str, Any]
            Dictionary of keyword arguments to pass to the task function.
        desc : str, optional
            Description for the progress bar and logging, by default "Running tasks".
        
        Returns
        -------
        None
        
        Notes
        -----
        - Uses ThreadPoolExecutor with max_workers defined in class initialization
        - Provides progress tracking via tqdm
        - Logs timing information and any exceptions that occur
        - Does not raise exceptions from individual tasks but logs them instead
        """
        
        start_time = time.time()

        logger.info(f"Active threads before: {threading.active_count()}")

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(task_fn, file, **task_args): file for file in self.files
            }
            logger.info(f"Active threads after submission: {threading.active_count()}")

            with tqdm(total=len(futures), desc=desc) as pbar:
                for future in as_completed(futures):
                    args = futures[future]
                    try:
                        future.result()
                    except Exception as e:
                        logger.error(f"Task failed with args {args}: {e}")
                    pbar.update(1)

        elapsed = time.time() - start_time
        logger.info(f"{desc} finished in {elapsed:.2f} seconds.")

class UnzipVCF(ParallelTaskRunner):
    """
    A class for unzipping VCF (Variant Call Format) files after imputation, with support for parallel processing.
    This class extends ParallelTaskRunner to efficiently extract VCF files from zip archives,
    including password-protected ones. It collects all zip files in the working directory
    and extracts their contents to the output directory.
    
    Attributes
    ----------
    Inherits all attributes from ParallelTaskRunner
    
    Methods
    -------
    execute_task(password=None)
        Execute the unzipping process for all zip files in the working directory
    unzip_files(zip_path, password)
        Extract files from a specific zip archive, handling password protection if needed

    Notes
    -----
    - VCF files are commonly used in genomics for storing gene sequence variations
    - The class only extracts files (not directories) from the zip archives
    - All extracted files are placed directly in the output directory without preserving paths
    - This class is designed for post-imputation processing in genetic data pipelines
    """

    def execute_task(self, password: Optional[str] = None) -> None:
        """
        Execute the post-imputation unzipping task on VCF files.

        This method performs the following steps:
        1. Collects all zip files in the working directory
        2. Unzips the VCF files, using the provided password if necessary

        Parameters
        ----------
        password : Optional[str]
            Password to decrypt zip files if they are password-protected.
            Default is None.

        Returns
        -------
        None
            This method doesn't return any value.
        """

        task_args ={'password': password}

        self._file_collector("*.zip")

        self._run_task(
            self.unzip_files,
            task_args=task_args,
            desc="Unzipping VCF files"
        )

        return
    
    def unzip_files(self, zip_path: Path, password: str) -> None:
        """
        Extract files from a password-protected zip archive.
        This method extracts all non-directory files from the specified zip archive
        to the class's output_path directory. If the zip file is password-protected,
        provide the password as a parameter.
        
        Parameters
        ----------
        zip_path : Path
            Path to the zip file to be extracted
        password : str
            Password for the zip file, None if the file is not password-protected
        
        Returns
        -------
        None
        
        Raises
        ------
        No exceptions are explicitly raised, but logging errors occur if the zip file is corrupted
        
        Notes
        -----
        Files are extracted to the output_path directory of the class instance.
        Only files (not directories) are extracted from the archive.
        File paths are not preserved - all files are placed directly in output_path.
        """
        
        with zipfile.ZipFile(zip_path, 'r') as zf:

            if password is not None:
                password_bytes = bytes(password, 'utf-8')
            else:
                password_bytes = None

            zf.setpassword(password_bytes)
            
            if zf.testzip() is not None:
                logger.error(f"Corrupted ZIP file: {zip_path}")
                return
            else:
                logger.info(f"Extracting {zip_path}...")
            
            for member in zf.namelist():
                
                if zf.getinfo(member).is_dir():
                    continue
                
                target_path = os.path.join(self.output_path, os.path.basename(member))
                
                with zf.open(member) as source, open(target_path, 'wb') as target:
                    target.write(source.read())

        logger.info(f"Successfully extracted {zip_path}")
        pass

class FilterVariants(ParallelTaskRunner):
    """
    A class for filtering genetic variants in VCF/BCF files based on imputation quality (R² statistic).
    This class extends ParallelTaskRunner to provide parallel processing capabilities for filtering
    variants across multiple VCF files. It identifies variants with imputation quality below a specified
    R² threshold and removes them from the output files.
    
    Attributes
    ----------
    Inherits all attributes from ParallelTaskRunner
    
    Methods
    -------
    execute_task(r2_threshold, output_prefix='filtered-')
        Set up and execute parallel filtering of variants based on R² threshold
    filter_variants(input_file, r2_threshold, output_prefix='filtered-')
        Filter a single VCF/BCF file based on the R² imputation quality threshold
    
    Dependencies
    -----------
    - bcftools must be installed and available in the system path
    - ParallelTaskRunner parent class for handling parallel execution
    
    Notes
    -----
    The class searches for files matching the pattern '*dose.vcf.gz' in the input directory
    and processes them in parallel. The filtered output files will be saved in the output
    directory with the specified prefix added to their original filenames.
    """

    def execute_task(self, r2_threshold: float, output_prefix: str = 'filtered-') -> None:
        """
        Execute the task of filtering variants based on an R² threshold.

        This method collects the necessary files with the pattern '*dose.vcf.gz' and runs 
        the filtering task with the specified parameters.

        Parameters
        ----------
        r2_threshold : float
            The threshold value for the R² statistic. Variants with an R² value below this threshold 
            will be filtered out.
        output_prefix : str, optional
            The prefix to be added to output filenames. Default is 'filtered-'.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If r2_threshold is not a float or output_prefix is not a string.

        Notes
        -----
        The method uses internal methods _file_collector and _run_task to perform the filtering operation.
        """

        task_args = {'r2_threshold': r2_threshold, 'output_prefix': output_prefix}
        if not isinstance(r2_threshold, float):
            raise TypeError(f"r2_threshold should be of type float, got {type(r2_threshold)}")
        if not isinstance(output_prefix, str):
            raise TypeError(f"prefix should be of type str, got {type(output_prefix)}")

        self._file_collector('*dose.vcf.gz')

        self._run_task(
            self.filter_variants,
            task_args=task_args,
            desc="Filter variants"
        )

        return
    
    def filter_variants(self, input_file: Path, r2_threshold: float, output_prefix: str = 'filtered-') -> None:
        """
        Filter variants from a VCF/BCF file based on R2 imputation quality threshold.
        This method takes an imputed VCF/BCF file and filters out variants with 
        imputation quality (R2) below the specified threshold. The filtered output 
        is saved as a compressed VCF.

        Parameters:
        -----------
        input_file (Path): 
            Path to the input VCF/BCF file to be filtered
        r2_threshold (float): 
            Minimum R2 imputation quality threshold (variants with R2 <= threshold will be removed)
        output_prefix (str, optional): 
            Prefix to add to the output filename. Defaults to 'filtered-'.
        
        Returns:
        --------
            None: The method outputs a filtered VCF file but doesn't return a value.
        
        Raises:
        -------
            FileExistsError: If the input file does not exist
            IsADirectoryError: If the input path is a directory, not a file
            TypeError: If r2_threshold is not a float or output_prefix is not a string
        
        Notes:
        ------
            The output file will be saved in the instance's output_path directory with
            the name constructed as: output_prefix + input_file.name
            This method requires bcftools to be installed and available in the system path.
        """

        if not input_file.exists():
            raise FileExistsError(f"Input file {input_file} does not exist")
        if not input_file.is_file():
            raise IsADirectoryError(f"Input file {input_file} is not a file")
        if not isinstance(r2_threshold, float):
            raise TypeError(f"r2_threshold should be of type float, got {type(r2_threshold)}")
        if not isinstance(output_prefix, str):
            raise TypeError(f"output_prefix should be of type str, got {type(output_prefix)}")
        
        output_file = self.output_path / (output_prefix + input_file.name)

        logger.info(f"Filtering {input_file} with R2 > {r2_threshold}")
        logger.info(f"Output file: {output_file}")

        # bcftools command
        bcf_cmd = [
            "bcftools", "view", "-Oz", "-i", f"R2>{r2_threshold}",
            input_file, "-o", output_file
        ]
        file_name = input_file.name
        
        try:
            # execute bcftools command
            subprocess.run(bcf_cmd, check=True)
            logger.info(f"Chromosome {file_name}: Completed")
        except subprocess.CalledProcessError as e:
            logger.error(f"Chromosome {file_name}: Failed with error {e}")
        except FileNotFoundError:
            logger.error(f"Chromosome {file_name}: Input file not found")
        pass

class NormalizeVCF(ParallelTaskRunner):
    """
    A class for normalizing VCF files post-imputation in parallel.
    This class provides functionality to process VCF files by normalizing them using
    bcftools. It's specifically designed to handle post-imputation VCF files and
    split multiallelic variants into separate entries.
    The class inherits from ParallelTaskRunner to enable parallel processing of
    multiple VCF files, which improves performance for large-scale genomic datasets.
    
    Attributes:
    ----------
        Inherits all attributes from ParallelTaskRunner

    Methods:
    -------
    execute_task(output_prefix='uncompressed-')
        Execute the normalization task on VCF files in parallel.
    normalize_vcf(input_file, output_prefix='uncompressed-')
        Normalize a single VCF file using bcftools norm with the -m -any option.
    """

    def execute_task(self, output_prefix: str = 'uncompressed-') -> None:
        """
        Execute the post-imputation normalization task on VCF files.

        This method collects filtered dose VCF files matching the pattern 'filtered-*dose.vcf.gz'
        and runs the normalization process on them. The normalized files will be prefixed with 
        the provided output_prefix.

        Parameters:
        -----------
            output_prefix (str, optional): Prefix to add to the output files. Defaults to 'uncompressed-'.

        Raises:
        -------
            TypeError: If output_prefix is not a string.

        Returns:
        --------
            None
        """

        task_args = {'output_prefix': output_prefix}
        if not isinstance(output_prefix, str):
            raise TypeError(f"prefix should be of type str, got {type(output_prefix)}")

        self._file_collector('filtered-*dose.vcf.gz')

        self._run_task(
            self.normalize_vcf,
            task_args=task_args,
            desc="Normalizing VCF files"
        )

        return
    
    def normalize_vcf(self, input_file: Path, output_prefix: str = 'uncompressed-') -> None:
        """
        Normalizes a VCF file using bcftools norm with the -m -any option.
        This method takes a VCF file, performs normalization using bcftools to split 
        multiallelic variants into separate entries, and outputs the normalized file 
        with the specified prefix.
        
        Parameters:
        -----------
            input_file (Path): Path to the input VCF file to be normalized
            output_prefix (str, optional): Prefix for the output file name. Defaults to 'uncompressed-'
        
        Returns:
        --------
            None

        Raises:
        -------
            FileExistsError: If the input file does not exist
            IsADirectoryError: If the input file path points to a directory
            TypeError: If output_prefix is not a string
        
        Note:
        -----
            The output file will be saved in the output_path directory with the naming 
            convention: output_prefix + base_name, where base_name is derived from the input file.
        """

        if not input_file.exists():
            raise FileExistsError(f"Input file {input_file} does not exist")
        if not input_file.is_file():
            raise IsADirectoryError(f"Input file {input_file} is not a file")
        if not isinstance(output_prefix, str):
            raise TypeError(f"output_prefix should be of type str, got {type(output_prefix)}")
        
        base_name = input_file.name.split('-')[1]

        output_file = self.output_path / (output_prefix + base_name)

        # Normalize with `bcftools norm -Ou -m -any`
        bcf_cmd = [
            "bcftools", "norm", "-Ou", "-o", output_file,"-m", "-any", input_file
        ]

        chr_number = os.path.basename(input_file)
        try:
            # execute bcftools command
            subprocess.run(bcf_cmd, check=True)
            logger.info(f"Chromosome {chr_number}: Completed")
        except subprocess.CalledProcessError as e:
            logger.error(f"Chromosome {chr_number}: Failed with error {e}")
        except FileNotFoundError:
            logger.error(f"Chromosome {chr_number}: Input file not found")
        pass

class ReferenceNormalizeVCF(ParallelTaskRunner):
    """
    A class for normalizing VCF files using a reference genome in parallel.
    This class extends ParallelTaskRunner to process multiple VCF files concurrently,
    normalizing them against a reference genome using bcftools. If a reference file
    is not provided, it will automatically download the appropriate reference genome
    based on the specified build.
    
    Attributes:
    -----------
        reference_file (str): Path to the reference genome file used for normalization.
    
    Methods:
    --------
        execute_task: Sets up normalization task and processes files in parallel.
        normalize_with_reference: Normalizes a single VCF file using bcftools.
    
    Example:
    --------
        normalizer = ReferenceNormalizeVCF()
        normalizer.execute_task(build='38', output_prefix='normalized-')
    """

    def execute_task(self, build: str = '38', output_prefix: str = 'normalized-', reference_file: Path = None) -> None:
        """
        Execute the post-imputation normalization task with reference genome.
        This method normalizes VCF files using a reference genome. If no reference file is provided,
        it automatically downloads the appropriate reference genome based on the build parameter.

        Parameters:
        -----------
            build (str, optional): Genome build version, either '37' or '38'. Defaults to '38'.
            output_prefix (str, optional): Prefix to add to output files. Defaults to 'normalized-'.
            reference_file (Path, optional): Path to the reference genome file. If None or the file doesn't exist,
                                            the reference will be downloaded. Defaults to None.
        Returns:
        --------
            None

        Raises:
        -------
            TypeError: If output_prefix is not a string.
        
        Notes:
        ------
            This method collects uncompressed dose VCF files using a pattern match and normalizes them
            against the reference genome. The downloaded reference genomes come from the 1000 Genomes Project.
        """

        task_args = {'output_prefix': output_prefix}
        if not isinstance(output_prefix, str):
            raise TypeError(f"prefix should be of type str, got {type(output_prefix)}")
        
        if not reference_file or not reference_file.exists():
            if build == '37':

                assemb37 = AssemblyReferenceFetcher(
                        base_url='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/',
                        build='37',
                        extension='.fa.gz'
                )
                assemb37.get_reference_url()
                assemb37.download_reference_file()
                assemb37.unzip_reference_file()
                reference_file = str(assemb37.file_path)

            elif build == '38':

                assemb38 = AssemblyReferenceFetcher(
                        base_url='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/',
                        build='38',
                        extension='.fa'
                )
                assemb38.get_reference_url()
                assemb38.download_reference_file()
                assemb38.unzip_reference_file()
                reference_file = str(assemb38.file_path)

        self._file_collector('uncompressed-*dose.vcf.gz')
        self.reference_file = reference_file

        self._run_task(
            self.normalize_with_reference,
            task_args=task_args,
            desc="Normalizing VCF files with Reference"
        )

        return
    
    def normalize_with_reference(self, input_file: Path, output_prefix: str = 'normalized-') -> None:
        """
        Normalize a VCF file with a reference genome using bcftools.
        This method takes an input VCF file and normalizes it against a reference genome
        using bcftools norm. The normalized output is compressed with gzip (-Oz).
        
        Parameters
        ----------
        input_file : Path
            Path to the input VCF file to be normalized.
        output_prefix : str, default='normalized-'
            Prefix to add to the output filename.
        
        Returns
        -------
        None
            The method doesn't return a value but creates a normalized VCF file
            at the output_path location.
        
        Raises
        ------
        TypeError
            If output_prefix is not a string.
        subprocess.CalledProcessError
            If the bcftools command fails.
        FileNotFoundError
            If the input file cannot be found.
        
        Notes
        -----
        The output filename is constructed from the output_prefix and the
        base name extracted from the input filename (after the first hyphen).
        """

        if not isinstance(output_prefix, str):
            raise TypeError(f"output_prefix should be of type str, got {type(output_prefix)}")
        
        base_name = input_file.name.split('-')[1]

        output_file = self.output_path / (output_prefix + base_name)

        # Normalize with `bcftools norm -Ou -m -any`
        bcf_cmd = [
            "bcftools", "norm", "-Oz", "-f", self.reference_file, "-o", output_file, input_file
        ]

        chr_number = os.path.basename(input_file)
        try:
            # execute bcftools command
            subprocess.run(bcf_cmd, check=True)
            logger.info(f"Chromosome {chr_number}: Completed")
        except subprocess.CalledProcessError as e:
            logger.error(f"Chromosome {chr_number}: Failed with error {e}")
        except FileNotFoundError:
            logger.error(f"Chromosome {chr_number}: Input file not found")
        pass

class IndexVCF(ParallelTaskRunner):

    def execute_task(self, pattern: str) -> None:

        task_args = dict()

        self._file_collector(pattern)
        self._run_task(
            self.index_vcf,
            task_args=task_args,
            desc="Indexing VCF files"
        )
        return
    
    def index_vcf(self, input_file: Path) -> None:

        if not os.path.isfile(input_file):
            raise FileExistsError(f"Input file {input_file} does not exist")

        # Build the bcftools index command
        command = ["bcftools", "index", input_file]

        # Execute the command
        try:
            subprocess.run(command, check=True)
            print(f"Successfully indexed: {input_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error indexing {input_file}: {e}")

        pass

class AnnotateVCF(ParallelTaskRunner):

    def execute_task(self, ref_annotation: Path, output_prefix: str = 'annotated-') -> None:

        task_args = {'ref_annotation': ref_annotation, 'output_prefix': output_prefix}
        if not isinstance(ref_annotation, Path):
            raise TypeError(f"ref_annotation should be of type Path, got {type(ref_annotation)}")
        if not isinstance(output_prefix, str):
            raise TypeError(f"prefix should be of type str, got {type(output_prefix)}")
        if not ref_annotation.exists():
            raise FileNotFoundError(f"Reference annotation file {ref_annotation} does not exist.")
        
        self._file_collector('normalized-*dose.vcf.gz')
        
        self._run_task(
            self.annotate_vcf,
            task_args=task_args,
            desc="Annotating VCF files"
        )

        return
    
    def annotate_vcf(self, input_file: Path, ref_annotation: Path, output_prefix: str = 'annotated-') -> None:

        if not input_file.exists():
            raise FileExistsError(f"Input file {input_file} does not exist")
        if not input_file.is_file():
            raise IsADirectoryError(f"Input file {input_file} is not a file")
        if not isinstance(ref_annotation, Path):
            raise TypeError(f"ref_annotation should be of type Path, got {type(ref_annotation)}")
        if not isinstance(output_prefix, str):
            raise TypeError(f"output_prefix should be of type str, got {type(output_prefix)}")
        
        base_name = input_file.name.split('-')[-1]

        output_file = self.output_path / (output_prefix + base_name)

        # Annotate with `bcftools annotate -a`
        bcf_cmd = [
            "bcftools", "annotate",
            "--annotations", ref_annotation,
            "--columns", "ID",
            "--output", output_file,
            input_file
        ]

        chr_number = os.path.basename(input_file)
        try:
            # execute bcftools command
            subprocess.run(bcf_cmd, check=True)
            logger.info(f"Chromosome {chr_number}: Completed")
        except subprocess.CalledProcessError as e:
            logger.error(f"Chromosome {chr_number}: Failed with error {e}")
        except FileNotFoundError:
            logger.error(f"Chromosome {chr_number}: Input file not found")
        pass

class ProcessVCF:
    """
    Class to run the entire post-imputation pipeline.
    """

    def __init__(self, input_path: Path, output_path: Path) -> None:
        
        if not isinstance(input_path, Path):
            raise TypeError(f"input_path should be of type Path, got {type(input_path)}")
        if not isinstance(output_path, Path):
            raise TypeError(f"output_path should be of type Path, got {type(output_path)}")
        if not input_path.exists():
            raise FileNotFoundError(f"Input path {input_path} does not exist.")
        if not output_path.exists():
            raise FileNotFoundError(f"Output path {output_path} does not exist.")
        if not input_path.is_dir():
            raise NotADirectoryError(f"Output path {input_path} is not a directory.")
        if not output_path.is_dir():
            raise NotADirectoryError(f"Output path {output_path} is not a directory.")
        
        self.input_path = input_path
        self.output_path= output_path

        self.process_vcf = self.input_path / 'process_vcf'
        self.process_vcf.mkdir(parents=True, exist_ok=True)

        pass

    def execute_unzip(self, password: str = None) -> None:
        
        unzipper = UnzipVCF(
            input_path = self.input_path,
            output_path=self.process_vcf
        )

        unzipper.execute_task(password=password)

        return

    def execute_filter(self, r2_threshold: float = 0.3) -> None:

        filter = FilterVariants(
            input_path = self.process_vcf,
            output_path=self.process_vcf
        )
        filter.execute_task(r2_threshold=r2_threshold)

        return
    
    def execute_normalize(self) -> None:
        
        normalizer = NormalizeVCF(
            input_path = self.process_vcf,
            output_path=self.process_vcf
        )
        normalizer.execute_task()

        return
    
    def execute_reference_normalize(self, build: str = '38', reference_file: Path = None) -> None:
        
        reference_normalizer = ReferenceNormalizeVCF(
            input_path = self.process_vcf,
            output_path=self.process_vcf
        )
        reference_normalizer.execute_task(build=build, reference_file=reference_file)

        return
    
    def execute_index(self, pattern: str = 'normalized-*dose.vcf.gz') -> None:

        indexer = IndexVCF(
            input_path = self.process_vcf,
            output_path=self.process_vcf
        )
        indexer.execute_task(pattern=pattern)

        return
    
    def execute_annotate(self, ref_annotation: Path, output_prefix: str = 'annotated-') -> None:
        
        annotator = AnnotateVCF(
            input_path = self.process_vcf,
            output_path=self.process_vcf
        )
        annotator.execute_task(ref_annotation=ref_annotation, output_prefix=output_prefix)

        return
    
    def execute_concatenate(self, output_name: str = 'final.vcf.gz', max_threads: int = None) -> None:
        """
        Concatenate all VCF files in the output directory into a single VCF file.
        """
        if not isinstance(output_name, str):
            raise TypeError(f"output_file should be of type str, got {type(output_file)}")
        
        output_path = self.output_path / output_name

        input_files = list(self.process_vcf.glob('annotated*.vcf.gz'))
        input_files.sort()
        if not input_files:
            raise FileNotFoundError(f"No VCF files annotated VCF files found in {self.process_vcf}")
        
        if not max_threads:
            max_threads = min(8, os.cpu_count())
        if max_threads < 1:
            raise ValueError(f"max_threads should be at least 1, got {max_threads}")

        # Concatenate with `bcftools concat`
        command = [
            "bcftools", "concat",
            *input_files,  # List of input files
            "--threads", str(max_threads),
            "-Oz",
            "-o", output_path
        ]

        try:
            subprocess.run(command, check=True)
            logger.info(f"Successfully concatenated and outputted to: {output_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error concatenating VCF files: {e}")
        pass