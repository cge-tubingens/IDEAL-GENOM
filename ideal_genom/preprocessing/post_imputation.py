import os
import subprocess
import psutil
import zipfile
import shutil

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

        # Loop through chromosomes 1 to 22
        os.makedirs(results_dir, exist_ok=True)

        # Convert the password to bytes
        password_bytes = bytes(password, 'utf-8')

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
                    try:
                        zf.setpassword(password_bytes)
                        zf.extractall(results_dir)
                        print(f"Successfully extracted: {zip_file} to {results_dir}")
                    except RuntimeError:
                        print(f"Wrong password for: {zip_file}")

            except zipfile.BadZipFile:
                print(f"Error: {zip_file} is not a valid ZIP file.")
            except RuntimeError as e:
                print(f"RuntimeError extracting {zip_file}: {e}")
            except Exception as e:
                print(f"Error extracting {zip_file}: {e}")

        def move_files_to_parent_and_remove_folders(parent_folder):
            """
            Moves all files from subfolders in the parent folder to the parent folder itself, 
            then removes the now-empty subfolders.

            Parameters:
            - parent_folder: str, the path to the folder containing the subfolders.
            """
            # Loop through all items in the parent folder
            for subfolder in os.listdir(parent_folder):
                subfolder_path = os.path.join(parent_folder, subfolder)

                # Check if it's a directory
                if os.path.isdir(subfolder_path):
                    # Move all files from the subfolder to the parent folder
                    for file_name in os.listdir(subfolder_path):
                        file_path = os.path.join(subfolder_path, file_name)
                        if os.path.isfile(file_path):  # Only move files, skip directories
                            shutil.move(file_path, os.path.join(parent_folder, file_name))

                    # Remove the now-empty subfolder
                    os.rmdir(subfolder_path)
                    print(f"Removed folder: {subfolder_path}")
        
        move_files_to_parent_and_remove_folders(results_dir)

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
