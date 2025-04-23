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
    Base class to run tasks in parallel using threads.
    """

    def __init__(self, input_path: Path, output_path: Path, max_workers: Optional[int] = None) -> None:

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

    def execute_task(self):
        """
        Execute the task. This method should be overridden by subclasses.
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def _file_collector(self, filename_pattern: str) -> List[Path]:
        """
        Collect files matching the filename pattern from the input path.
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
        Run a list of tasks in parallel using threads.

        Args:
            task_fn: A callable to run, e.g., a processing function.
            task_args: A list of tuples, each containing the arguments for one task.
            desc: Description for the progress bar.
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

    def execute_task(self, password: Optional[str] = None) -> None:

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
        Unzips a password-protected zip file into a specified output folder.
        Parameters:
        -----------
            zip_path (str): Path to the zip file
            output_folder (str): Folder where files will be extracted
            password (str): Password to decrypt the zip file
        Raises:
        -------
            zipfile.BadZipFile: If zip file is corrupted
            RuntimeError: If password is incorrect
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

    def execute_task(self, r2_threshold: float, output_prefix: str = 'filtered-') -> None:

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

    def execute_task(self, output_prefix: str = 'uncompressed-') -> None:

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

    def execute_task(self, build: str = '38', output_prefix: str = 'normalized-', reference_file: Path = None) -> None:

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