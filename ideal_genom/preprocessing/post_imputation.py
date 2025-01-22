import os
import subprocess
import psutil
import zipfile
import shutil
from tqdm import tqdm
import time
import threading

from concurrent.futures import ThreadPoolExecutor, as_completed

class PostImputation:

    def __init__(self, input_path:str, output_path:str, output_name:str, dependables:str)->None:
        
        """
        Initializes the AfterImputation class with the given input, output, and dependables directories.

        The goal of this class is to provide a set of functions to process VCF files after imputation and prepare them for downstream analyses. The final output will be binary PLINK files that can be used for association studies.

        Parameters:
        -----------
            input_path (str): The path to the input directory.
            output_path (str): The path to the output directory.
            dependables (str): The path to the dependables directory.

        Raises:
        -------
            ValueError: If input_path is None or not a string.
            ValueError: If output_path is None or not a string.
            ValueError: If dependables is not a string.
            ValueError: If input_path does not exist.
            ValueError: If output_path does not exist.
            ValueError: If dependables does not exist.
        """

        if input_path is None:
            raise ValueError("Input directory is not provided")
        if not isinstance(input_path, str):
            raise ValueError("Input directory should be a string")
        if output_path is None:
            raise ValueError("Output directory is not provided")
        if not isinstance(output_path, str):
            raise ValueError("Output directory should be a string")
        if not isinstance(dependables, str):
            raise ValueError("Dependables directory should be a string")
        if output_name is None:
            raise ValueError("Output name is not provided")
        elif output_name == "":
            raise ValueError("Output name cannot be an empty string")
        if not isinstance(output_name, str):
            raise ValueError("Output name should be a string")
        
        if os.path.isdir(input_path) is False:
            raise ValueError("Input directory does not exist")
        if os.path.isdir(output_path) is False:
            raise ValueError("Output directory does not exist")
        if os.path.isdir(dependables) is False:
            raise ValueError("Dependables directory does not exist")
        
        self.input_path = input_path
        self.output_path= output_path
        self.dependables= dependables
        self.output_name= output_name

        self.results_dir = os.path.join(output_path, 'post_imputation')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        self.ready_dir = os.path.join(self.results_dir, 'analysis_ready')
        if not os.path.exists(self.ready_dir):
            os.mkdir(self.ready_dir)

        pass

    def execute_unzip_chromosome_files(self, password:str, max_workers:int=None)->None:
        """
        Unzips chr*.zip files (chr1.zip to chr22.zip) from the input folder and extracts them into the output folder.

        Parameters:
        -----------
            password (str): The password to unzip the files.
            max_workers (int, optional): The maximum number of worker threads to use for parallel processing. Defaults to None.

        Returns:
        --------
            None

        Raises:
        -------
            ValueError: If password is not provided or is not a string.
        """

        input_folder = self.input_path
        results_dir  = self.results_dir

        if password is None:
            raise ValueError("Password is not provided")
        if not isinstance(password, str):
            raise ValueError("Password should be a string")
        if password == "":
            raise ValueError("Password cannot be an empty string")

        # List of chromosome zip files
        zip_files = [os.path.join(input_folder, f"chr{chr_num}.zip") for chr_num in range(20, 23)]

        self.parallel_unzip_with_password(
            zip_files    =zip_files,
            output_folder=results_dir,
            password     =password,
            max_workers =max_workers
        )
        
        pass

    def execute_filter_variants(self, r2_threshold:float, max_workers:int=None)->None:
        
        """
        Filters genetic variants based on the R2 threshold and processes chromosomes in parallel.

        Parameters:
        -----------
            r2_threshold (float): The R2 threshold for filtering variants. Must be between 0 and 1.
            max_workers (int): The maximum number of worker threads to use for parallel processing. If None, defaults to the number of CPU cores minus 2.

        Raises:
        -------
            ValueError: If r2_threshold is not provided, not a float, or not between 0 and 1. If max_workers is not a positive integer.
        
        Returns:
        --------
            None
        """

        if r2_threshold is None:
            raise ValueError("R2 threshold is not provided")
        if not isinstance(r2_threshold, float):
            raise ValueError("R2 threshold should be a float")
        if r2_threshold < 0 or r2_threshold > 1:
            raise ValueError("R2 threshold should be between 0 and 1")

        results_dir = self.results_dir

        # List of chromosome files
        chr_files = [os.path.join(results_dir, f"chr{chr_num}.dose.vcf.gz") for chr_num in range(20, 23)]

        self.parallel_filter_chromosome(
            chr_files    =chr_files,
            output_folder=results_dir,
            r2_threshold =r2_threshold,
            max_workers  =max_workers
        )

        pass

    def execute_normalize_vcf(self, max_workers:int=None)->None:

        results_dir = self.results_dir

        chr_files = [
            os.path.join(results_dir, f"filtered_chr{chr_num}.dose.vcf.gz") for chr_num in range(20, 23)
        ]

        self.parallel_normalize_chromosome(
            chr_files    =chr_files,
            output_folder=results_dir,
            max_workers  =max_workers
        )

        pass

    def execute_normalize_with_reference(self, reference_genome:str, max_workers:int)->None:

        results_dir = self.results_dir
        dependables = self.dependables

        chr_files = [
            os.path.join(results_dir, f"uncompressed_chr{chr_num}.dose.vcf.gz") for chr_num in range(20, 23)
        ]

        self.parallel_normalize_chromosome_with_reference(
            chr_files    =chr_files,
            output_folder=results_dir,
            reference_file=os.path.join(dependables, reference_genome),
            max_workers  =max_workers
        )

        pass

    def execute_index_vcf(self)->None:

        """
        Executes the indexing of VCF files for chromosomes 1 to 22 using bcftools.
        
        This method checks for the existence of index files (.csi or .tbi) for each chromosome's VCF file. If an index file already exists, it skips the indexing for that chromosome. Otherwise, it builds and executes the bcftools index command to create the index file.
        
        Parameters:
        -----------
            None
        
        Returns:
        --------
            None

        Raises:
        -------
            subprocess.CalledProcessError: If the bcftools index command fails.
        """

        results_dir = self.results_dir

        for chr_num in range(1, 23):
            
            input_vcf = f"normalized_chr{chr_num}.dose.vcf.gz"

            input_file = os.path.join(results_dir, input_vcf)

            if os.path.isfile(input_file+'.csi') or os.path.isfile(input_file+'.tbi'):
                print(f"Index file for {input_vcf} already exists.")
                continue

            # Build the bcftools index command
            command = ["bcftools", "index", input_file]

            # Execute the command
            try:
                subprocess.run(command, check=True)
                print(f"Successfully indexed: {input_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"Error indexing {input_vcf}: {e}")

        pass

    def execute_annotate_vcf(self, annotations_file:str)->None:
        """
        Annotates VCF files for chromosomes 1 to 22 using bcftools.

        This method loops over chromosomes 1 to 22 and annotates the corresponding VCF files
        using the provided annotations file. The annotated VCF files are saved in the results directory.

        Parameters:
        -----------
            annotations_file (str): The path to the annotations file to be used for annotating the VCF files. it is assumed that the annotations file is in the dependables directory and has the format .vcf.gz. If the file it is not indexed, it will be indexed before annotation.

        Returns:
        --------
            None

        Raises:
        -------
            subprocess.CalledProcessError: If the bcftools command fails for any chromosome.
        """


        results_dir = self.results_dir
        dependables = self.dependables

        annotations = os.path.join(dependables, annotations_file)

        if not os.path.exists(annotations+'.tbi') or not os.path.exists(annotations+'.csi'):
            subprocess.run(["bcftools", "index", annotations], check=True)

        for chr_num in range(1, 23):  # Loop over chromosomes 1 to 22

            input_vcf = f"normalized_chr{chr_num}.dose.vcf.gz"
            output_vcf= f"annotated_normalized_chr{chr_num}.dose.vcf.gz"

            input_file = os.path.join(results_dir, input_vcf)
            output_file= os.path.join(results_dir, output_vcf)

            # bcftools command
            command = [
                "bcftools", "annotate",
                "--annotations", annotations,
                "--columns", "ID",
                "--output", output_file,
                input_file
            ]

            # execute bcftools command
            try:
                subprocess.run(command, check=True)
                print(f"Successfully annotated: {output_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"Error annotating {input_vcf}: {e}")
        pass

    
    def execute_concat_vcf(self):
        """
        Concatenates VCF files from chromosome 1 to 22 into a single VCF file.

        This method generates a list of input VCF files for chromosomes 1 through 22,
        then uses the `bcftools concat` command to concatenate them into a single
        output VCF file. The output file is saved in the results directory specified
        by `self.results_dir`.

        The method utilizes multiple threads for the concatenation process, with the
        number of threads being the total CPU count minus four.

        Parameters:
        -----------
            None

        Returns:
        --------
            None

        Raises:
        -------
            subprocess.CalledProcessError: If the `bcftools concat` command fails.
        """

        results_dir = self.results_dir

        # Create a list of input files (chr1 to chr22)
        input_files = [
            os.path.join(results_dir, f"annotated_normalized_chr{chr_num}.dose.vcf.gz") for chr_num in range(1, 23)
        ]

        output_file = os.path.join(results_dir, "annotated_normalized_combined_1_22.vcf.gz")

        max_threads = os.cpu_count() - 4

        # bcftools command
        command = [
            "bcftools", "concat",
            *input_files,  # List of input files
            "--threads", str(max_threads),
            "-Oz",
            "-o", output_file
        ]

        # execute the command
        try:
            subprocess.run(command, check=True)
            print(f"Successfully concatenated and outputted to: {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error concatenating VCF files: {e}")

    def get_plink_files(self, threads:int=None, memory:int=None)->None:
        """
        Generates PLINK binary files from a VCF file using PLINK2.

        This method constructs and executes a PLINK2 command to convert a VCF file into PLINK binary files (.bed, .bim, .fam). The number of threads and the amount of memory to be used by PLINK2 can be specified. If not provided, the method will use default values based on the system's available resources.

        Parameters:
        -----------
        threads (int, optional): The number of threads to use for PLINK2. If not provided,
                     defaults to the number of CPU cores minus two, or 10 if
                     the number of CPU cores cannot be determined.
        memory (int, optional): The amount of memory (in MB) to allocate for PLINK2. If not
                    provided, defaults to two-thirds of the available system memory.
        
        Returns:
        --------
            None
        """
        

        results_dir = self.results_dir
        output_name = self.output_name
        ready_dir   = self.ready_dir

        if threads is None:
            if os.cpu_count() is not None:
                threads = os.cpu_count() - 2
            else:
                threads = 10

        if memory is None:
            # get virtual memory details
            memory_info = psutil.virtual_memory()
            available_memory_mb = memory_info.available / (1024 * 1024)
            memory = round(2*available_memory_mb/3,0)

        input_vcf = os.path.join(results_dir, "annotated_normalized_combined_1_22.vcf.gz")

        # plink2 command
        command = [
            "plink2",
            "--vcf", input_vcf,
            "--snps-only", "just-acgt",
            "--make-bed",
            "--out", os.path.join(ready_dir, output_name),
            "--threads", str(threads),
            "--memory", str(memory)
        ]

        # execute plink2 command
        try:
            subprocess.run(command, check=True)
            print(f"PLINK2 command executed successfully. Output files saved with prefix: {output_name}")
        except subprocess.CalledProcessError as e:
            print(f"Error running PLINK2: {e}")

        pass

    @staticmethod
    def unzip_file_with_password(zip_path: str, output_folder: str, password: str) -> None:
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
        
        password_bytes = bytes(password, 'utf-8')
        
        with zipfile.ZipFile(zip_path, 'r') as zf:

            zf.setpassword(password_bytes)
            
            if zf.testzip() is not None:
                print(f"Corrupted ZIP file: {zip_path}")
                return
            else:
                print(f"Extracting {zip_path}...")
            
            for member in zf.namelist():
                
                if zf.getinfo(member).is_dir():
                    continue
                
                target_path = os.path.join(output_folder, os.path.basename(member))
                
                with zf.open(member) as source, open(target_path, 'wb') as target:
                    target.write(source.read())

        pass

    @staticmethod
    def parallel_unzip_with_password(zip_files: list, output_folder: str, password: str, max_workers: int = None) -> None:
        """
        Unzips multiple password-protected zip files in parallel.

        Parameters:
        -----------
            zip_files (list): List of paths to zip files
            output_folder (str): Folder where files will be extracted
            password (str): Password to decrypt the zip files
            max_workers (int): Maximum number of parallel workers. Defaults to min(8, CPU cores)
        """

        if max_workers is None:
            max_workers = min(8, os.cpu_count())

        start_time = time.time()

        print(f"Active threads before: {threading.active_count()}")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    PostImputation.unzip_file_with_password, 
                    zip_file, 
                    output_folder, 
                    password
                ): zip_file for zip_file in zip_files
            }
            print(f"Active threads after submission: {threading.active_count()}")

            with tqdm(total=len(futures), desc="Extracting files") as pbar:

                for future in as_completed(futures):

                    zip_file = futures[future]
                    
                    pbar.update(1)
                    pbar.set_description(f"🔄 Active jobs: {len(futures) - pbar.n}")
                    
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Error with {zip_file}: {str(e)}")

        total_time = time.time() - start_time
        print(f"\nTotal execution time: {total_time:.2f}s")

        pass

    @staticmethod
    def filter_chromosome(input_chr_file:str, output_chr_file:str, r2_threshold:float):
        
        if not os.path.isfile(input_chr_file):
            raise FileExistsError(f"Input file {input_chr_file} does not exist")

        # bcftools command
        bcf_cmd = [
            "bcftools", "view", "-Oz", "-i", f"R2>{r2_threshold}",
            input_chr_file, "-o", output_chr_file
        ]
        chr_number = os.path.basename(input_chr_file)
        try:
            # execute bcftools command
            subprocess.run(bcf_cmd, check=True)
            print(f"Chromosome {chr_number}: Completed")
        except subprocess.CalledProcessError as e:
            print(f"Chromosome {chr_number}: Failed with error {e}")
        except FileNotFoundError:
            print(f"Chromosome {chr_number}: Input file not found")
        pass

    @staticmethod
    def parallel_filter_chromosome(chr_files:list, output_folder:str, r2_threshold:float, max_workers:int=None)->None:

        if max_workers is None:
            max_workers = min(8, os.cpu_count())

        start_time = time.time()

        print(f"Active threads before: {threading.active_count()}")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    PostImputation.filter_chromosome, 
                    chr_file, 
                    os.path.join(output_folder, f"filtered_{os.path.basename(chr_file)}"),
                    r2_threshold
                ): chr_file for chr_file in chr_files
            }
            print(f"Active threads after submission: {threading.active_count()}")

            with tqdm(total=len(futures), desc="Filtering chromosomes") as pbar:

                for future in as_completed(futures):
                    chr_file = futures[future]
                    pbar.update(1)
                    pbar.set_description(f"🔄 Active jobs: {len(futures) - pbar.n}")
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Error with {chr_file}: {str(e)}")

        total_time = time.time() - start_time
        print(f"\nTotal execution time: {total_time:.2f}s")
        pass

    @staticmethod
    def normalize_chromosome(input_chr_file:str, output_chr_file:str)->None:

        if not os.path.isfile(input_chr_file):
            raise FileExistsError(f"Input file {input_chr_file} does not exist")
        
        # Normalize with `bcftools norm -Ou -m -any`
        bcf_cmd = [
            "bcftools", "norm", "-Ou", "-o", output_chr_file,"-m", "-any", input_chr_file
        ]

        chr_number = os.path.basename(input_chr_file)
        try:
            # execute bcftools command
            subprocess.run(bcf_cmd, check=True)
            print(f"Chromosome {chr_number}: Completed")
        except subprocess.CalledProcessError as e:
            print(f"Chromosome {chr_number}: Failed with error {e}")
        except FileNotFoundError:
            print(f"Chromosome {chr_number}: Input file not found")
        pass

    @staticmethod
    def parallel_normalize_chromosome(chr_files:list, output_folder:str, max_workers:int)->None:

        if max_workers is None:
            max_workers = min(8, os.cpu_count())

        start_time = time.time()

        print(f"Active threads before: {threading.active_count()}")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    PostImputation.normalize_chromosome, 
                    chr_file, 
                    os.path.join(output_folder, f"uncompressed_{os.path.basename(chr_file)}")
                ): chr_file for chr_file in chr_files
            }
            print(f"Active threads after submission: {threading.active_count()}")

            with tqdm(total=len(futures), desc="Filtering chromosomes") as pbar:

                for future in as_completed(futures):
                    chr_file = futures[future]
                    pbar.update(1)
                    pbar.set_description(f"🔄 Active jobs: {len(pbar.n)}")
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Error with {chr_file}: {str(e)}")

        total_time = time.time() - start_time
        print(f"\nTotal execution time: {total_time:.2f}s")
        pass

    @staticmethod
    def normalize_chromosome_with_reference(input_chr_file:str, output_chr_file:str, reference_file:str)->None:
        
        if not os.path.isfile(input_chr_file):
            raise FileExistsError(f"Input file {input_chr_file} does not exist")
        if not os.path.isfile(reference_file):
            raise FileExistsError(f"Reference file {reference_file} does not exist")
        
        # Normalize with reference genome and output to a new file
        bcf_cmd = [
            "bcftools", "norm", "-Oz", "-f", reference_file, "-o", output_chr_file, input_chr_file
        ]

        chr_number = os.path.basename(input_chr_file)
        try:
            # execute bcftools command
            subprocess.run(bcf_cmd, check=True)
            print(f"Chromosome {chr_number}: Completed")
        except subprocess.CalledProcessError as e:
            print(f"Chromosome {chr_number}: Failed with error {e}")
        except FileNotFoundError:
            print(f"Chromosome {chr_number}: Input file not found")
        pass

    @staticmethod
    def parallel_normalize_chromosome_with_reference(chr_files:list, output_folder:str, reference_file:str, max_workers:int=None)->None:

        if max_workers is None:
            max_workers = min(8, os.cpu_count())

        start_time = time.time()

        print(f"Active threads before: {threading.active_count()}")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    PostImputation.normalize_chromosome_with_reference, 
                    chr_file, 
                    os.path.join(output_folder, f"normalized_{os.path.basename(chr_file)}"),
                    reference_file
                ): chr_file for chr_file in chr_files
            }
            print(f"Active threads after submission: {threading.active_count()}")

            with tqdm(total=len(futures), desc="Filtering chromosomes") as pbar:

                for future in as_completed(futures):
                    chr_file = futures[future]
                    pbar.update(1)
                    pbar.set_description(f"🔄 Active jobs: {len(futures) - pbar.n}")
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Error with {chr_file}: {str(e)}")

        total_time = time.time() - start_time
        print(f"\nTotal execution time: {total_time:.2f}s")
        pass