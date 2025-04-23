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