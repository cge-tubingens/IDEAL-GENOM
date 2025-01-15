import os
import subprocess
import psutil
import zipfile

from concurrent.futures import ThreadPoolExecutor

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

        self.ready_dir = os.path.join(output_path, 'analysis_ready')
        if not os.path.exists(self.ready_dir):
            os.mkdir(self.ready_dir)

        pass

    def execute_unzip_chromosome_files(self, password:str)->None:
        """
        Unzips chr*.zip files (chr1.zip to chr22.zip) from the input folder and extracts them into the output folder.

        Returns:
        --------
            None
        """

        input_folder = self.input_path
        results_dir  = self.results_dir

        password_bytes = password.encode('utf-8')

        # Loop through chromosomes 1 to 22
        for chr_num in range(1, 23):
            zip_file = os.path.join(input_folder, f"chr{chr_num}.zip")  # Path to the ZIP file

            # Check if the ZIP file exists
            if not os.path.isfile(zip_file):
                print(f"File not found: {zip_file}")
                continue
            
            # Extract the ZIP file
            try:
                with zipfile.ZipFile(zip_file, 'r') as zf:
                    # Check if the ZIP file is encrypted
                    if zf.pwd is None and zf.testzip() is None:
                        zf.extractall(results_dir, pwd=password_bytes)
                        print(f"Successfully extracted: {zip_file} to {results_dir}")
                    else:
                        print(f"Error: {zip_file} appears to require a different password.")
            except zipfile.BadZipFile:
                print(f"Error: {zip_file} is not a valid ZIP file.")
            except RuntimeError as e:
                print(f"RuntimeError extracting {zip_file}: {e}")
            except Exception as e:
                print(f"Error extracting {zip_file}: {e}")

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
        
        if max_workers is None:
            max_workers = os.cpu_count()-2
        if isinstance(max_workers, float):
            max_workers = max(round(max_workers,0),1)
        elif not isinstance(max_workers, int):
            raise ValueError("Max workers should be a positive integer")

        results_dir = self.results_dir

        def process_chromosome(chr_number:int, input_folder:str, output_folder:str, r2_threshold:float):
            """
            Processes a chromosome VCF file by filtering variants with R2 > r2_threshold using bcftools.

            Parameters:
            -----------
                chr_number (int): The chromosome number to process.
                input_folder (str): The folder containing the input VCF files.
                output_folder (str): The folder where the filtered VCF files will be saved.
                r2_threshold (float): The R2 threshold for filtering variants.

            Raises:
            -------
                subprocess.CalledProcessError: If the bcftools command fails.
                FileNotFoundError: If the input file is not found.

            Prints:
                A message indicating whether the processing of the chromosome was completed or failed.
            """

            input_file = os.path.join(input_folder, f"chr{chr_number}.dose.vcf.gz")
            output_file = os.path.join(output_folder, f"filtered_chr{chr_number}.dose.vcf.gz")

            # bcftools command
            command = [
                "bcftools", "view", "-Oz", "-i", f"R2>{r2_threshold}",
                input_file, "-o", output_file
            ]
            try:
                # execute bcftools command
                subprocess.run(command, check=True)
                print(f"Chromosome {chr_number}: Completed")
            except subprocess.CalledProcessError as e:
                print(f"Chromosome {chr_number}: Failed with error {e}")
            except FileNotFoundError:
                print(f"Chromosome {chr_number}: Input file not found")
            pass

        chromosomes = range(1, 23)

        # Adjust `max_workers` for parallelism
        with ThreadPoolExecutor(max_workers=max_workers) as executor:  
            executor.map(lambda chr: process_chromosome(chr, results_dir, results_dir, r2_threshold), chromosomes)

        print("All processes completed successfully.")
        pass

    def execute_normalize_vcf(self, reference_genome:str)->None:
        
        """
        Normalize VCF files for each chromosome using bcftools.

        This method normalizes VCF files for chromosomes 1 through 22 using the bcftools tool.
        It performs two normalization steps for each chromosome:
        1. Normalize with `bcftools norm -Ou -m -any`.
        2. Normalize with a reference genome and output to a new file.

        Parameters:
        -----------
            reference_genome (str): The filename of the reference genome to use for normalization.

        Raises:
        -------
            Exception: If there is an error during the normalization process for any chromosome.
        
        Returns:
        --------
            None
        """

        results_dir = self.results_dir
        dependables = self.dependables

        for chr_number in range(1, 23):

            input_file = os.path.join(results_dir, f"filtered_chr{chr_number}.dose.vcf.gz")
            temp_file = os.path.join(results_dir, f"uncompressed_chr{chr_number}.dose.vcf.gz")
            output_file = os.path.join(results_dir, f"normalized_chr{chr_number}.dose.vcf.gz")
            reference_file = os.path.join(dependables, reference_genome)

            # Step 1: Normalize with `bcftools norm -Ou -m -any`
            norm_command_1 = [
                "bcftools", "norm", "-Ou", "-o", temp_file,"-m", "-any", input_file
            ]

            # Step 2: Normalize with reference genome and output to a new file
            norm_command_2 = [
                "bcftools", "norm", "-Oz", "-f", reference_file, "-o", output_file, temp_file
            ]

            try:
                # Step 1: Run bcftools norm -Ou -m -any
                norm_proc_1 = subprocess.Popen(norm_command_1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                norm_proc_1.wait()  # Wait for this process to finish

                # Check for errors in Step 1
                if norm_proc_1.returncode != 0:
                    raise Exception(f"Error in normalizing chromosome {chr_number} (Step 1): {norm_proc_1.stderr.read().decode()}")

                # Step 2: Run bcftools norm -Oz with the reference genome
                subprocess.run(norm_command_2, check=True)
                print(f"Chromosome {chr_number}: Normalized successfully")

            except subprocess.CalledProcessError as e:
                print(f"Chromosome {chr_number}: Error during processing: {e}")

        pass

    def execute_annotate_vcf(self, annotations_file:str)->None:
        """
        Annotates VCF files with dbSNP IDs from the given annotation file.

        Parameters:
        ----------
            annotations_file (str): path to the annotations VCF file (e.g., dbSNP file).
            input_vcf_pattern (str): pattern for input VCF files (e.g., 'normalized_chr*.dose.vcf.gz').
            output_vcf_pattern (str): pattern for output VCF files (e.g., 'annotated_normalized_chr*.dose.vcf.gz').
        """

        results_dir = self.results_dir
        dependables = self.dependables

        for chr_num in range(1, 23):  # Loop over chromosomes 1 to 22

            input_vcf = f"normalized_chr{chr_num}.dose.vcf.gz"
            output_vcf= f"annotated_normalized_chr{chr_num}.dose.vcf.gz"

            input_file = os.path.join(results_dir, input_vcf)
            output_file= os.path.join(results_dir, output_vcf)

            annotations = os.path.join(dependables, annotations_file)

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

    def execute_index_vcf(self)->None:

        """
        Index VCF files for chromosomes 1 to 22 using bcftools.

        This method loops over chromosome numbers from 1 to 22, constructs the 
        corresponding VCF file name, and uses the bcftools index command to index each VCF file. If the indexing is successful, a success message is printed. 
        
        If an error occurs during indexing, an error message is printed.
        
        Raises:
        -------
            subprocess.CalledProcessError: If the bcftools index command fails.

        Returns:
        --------
            None
        """

        results_dir = self.results_dir

        for chr_num in range(1, 23):  # Loop over chromosomes 1 to 22
            
            input_vcf = f"annotated_normalized_chr{chr_num}.dose.vcf.gz"

            input_file = os.path.join(results_dir, input_vcf)

            # Build the bcftools index command
            command = ["bcftools", "index", input_file]

            # Execute the command
            try:
                subprocess.run(command, check=True)
                print(f"Successfully indexed: {input_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"Error indexing {input_vcf}: {e}")

        pass
    
    def execute_concat_vcf(self):

        results_dir = self.results_dir

        # Create a list of input files (chr1 to chr22)
        input_files = [f"annotated_normalized_chr{chr_num}.dose.vcf.gz" for chr_num in range(1, 23)]

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
        Runs the PLINK2 command to process a VCF file and generate binary files.

        Parameters:
        -----------
            threads (int, optional): The number of threads to use for PLINK2.
            memory (int, optional): The amount of memory (in MB) to allocate for PLINK2.

        Returns:
        --------
            None

        Raises:
        -------
            subprocess.CalledProcessError: If the PLINK2 command fails.
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
