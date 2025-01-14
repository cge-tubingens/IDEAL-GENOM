import os
import subprocess
import psutil
import zipfile

from concurrent.futures import ThreadPoolExecutor

class AfterImputation:

    def __init__(self, input_path:str, output_path:str, dependables:str)->None:
        
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
        
        if os.path.isdir(input_path) is False:
            raise ValueError("Input directory does not exist")
        if os.path.isdir(output_path) is False:
            raise ValueError("Output directory does not exist")
        if os.path.isdir(dependables) is False:
            raise ValueError("Dependables directory does not exist")
        
        self.input_path = input_path
        self.output_path= output_path
        self.dependables = dependables

        self.results_dir = os.path.join(output_path, 'post_imputation')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        pass

    def filter_variants(self, r2_threshold:float, max_workers:int)->None:

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

        input_path = self.input_path
        resuts_dir = self.results_dir

        def process_chromosome(chr_number, input_folder, output_folder):

            input_file = os.path.join(input_folder, f"chr{chr_number}.dose.vcf.gz")
            output_file = os.path.join(output_folder, f"chr{chr_number}_R2_0.3.dose.vcf.gz")

            # Run the bcftools command
            command = [
                "bcftools", "view", "-Oz", "-i", "R2>0.3",
                input_file, "-o", output_file
            ]
            try:
                subprocess.run(command, check=True)
                print(f"Chromosome {chr_number}: Completed")
            except subprocess.CalledProcessError as e:
                print(f"Chromosome {chr_number}: Failed with error {e}")
            except FileNotFoundError:
                print(f"Chromosome {chr_number}: Input file not found")

        chromosomes = range(1, 23)

        with ThreadPoolExecutor(max_workers=os.cpu_count()-2) as executor:  # Adjust `max_workers` for parallelism
            executor.map(lambda chr: process_chromosome(chr, input_path, resuts_dir), chromosomes)

        print("All processes completed successfully.")
        pass

    def normalize_vcf(self, reference_genome):
        """
        This function normalizes VCF files for chromosomes 1 to 22 using `bcftools`.

        Parameters:
        - input_folder (str): Path to the input folder containing VCF files (chr*_R2_0.3.dose.vcf.gz).
        - output_folder (str): Path to the output folder for storing normalized VCF files.
        - reference_genome (str): Path to the reference genome file (e.g., GRCh38).
        """

        results_dir = self.results_dir
        dependables = self.dependables

        # Loop through chromosomes 1 to 22
        for chr_number in range(1, 23):
            # Input and output file paths
            input_file = os.path.join(results_dir, f"chr{chr_number}_R2_0.3.dose.vcf.gz")
            temp_file = os.path.join(results_dir, f"uncompressed_chr{chr_number}_R2_0.3.dose.vcf.gz")
            output_file = os.path.join(results_dir, f"normalized_chr{chr_number}_R2_0.3.dose.vcf.gz")
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

    def annotate_vcf_files(self, annotations_file:str):
        """
        Annotates VCF files with dbSNP IDs from the given annotation file.

        Parameters:
        - annotations_file: str, path to the annotations VCF file (e.g., dbSNP file).
        - input_vcf_pattern: str, pattern for input VCF files (e.g., 'normalized_chr*_R2_0.3.dose.vcf.gz').
        - output_vcf_pattern: str, pattern for output VCF files (e.g., 'annotated_normalized_chr*_R2_0.3.dose.vcf.gz').
        """

        results_dir = self.results_dir
        dependables = self.dependables

        for chr_num in range(1, 23):  # Loop over chromosomes 1 to 22

            input_vcf = f"normalized_chr{chr_num}_R2_0.3.dose.vcf.gz"
            output_vcf = f"annotated_normalized_chr{chr_num}_R2_0.3.dose.vcf.gz"

            input_file = os.path.join(results_dir, input_vcf)
            output_file = os.path.join(results_dir, output_vcf)

            annotations = os.path.join(dependables, annotations_file)

            # Build the bcftools command
            command = [
                "bcftools", "annotate",
                "--annotations", annotations,
                "--columns", "ID",
                "--output", output_file,
                input_file
            ]

            # Execute the command
            try:
                subprocess.run(command, check=True)
                print(f"Successfully annotated: {output_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"Error annotating {input_vcf}: {e}")

    def index_vcf_files(self):
        """
        Indexes VCF files using bcftools index.

        Parameters:
        - input_vcf_pattern: str, pattern for input VCF files (e.g., 'annotated_normalized_chr*_R2_0.3.dose.vcf.gz').
        """

        results_dir = self.results_dir

        for chr_num in range(1, 23):  # Loop over chromosomes 1 to 22
            
            input_vcf = f"annotated_normalized_chr{chr_num}_R2_0.3.dose.vcf.gz"

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
    
    def concat_vcf_files(self):
        """
        Concatenates multiple VCF files using bcftools concat.

        Parameters:
        - input_vcf_pattern: str, pattern for input VCF files (e.g., 'annotated_normalized_chr1_R2_0.3.dose.vcf.gz').
        - output_vcf: str, output file for the concatenated VCF (e.g., 'annotated_normalized_combined_1_22.vcf.gz').
        - threads: int, number of threads to use (default is 28).
        """

        results_dir = self.results_dir

        # Create a list of input files (chr1 to chr22)
        input_files = [f"annotated_normalized_chr{chr_num}_R2_0.3.dose.vcf.gz" for chr_num in range(1, 23)]

        output_file = os.path.join(results_dir, "annotated_normalized_combined_1_22.vcf.gz")

        max_threads = os.cpu_count() - 4

        # Build the bcftools concat command
        command = [
            "bcftools", "concat",
            *input_files,  # List of input files
            "--threads", str(max_threads),
            "-Oz",  # Output in compressed VCF format
            "-o", output_file  # Output file
        ]

        # Execute the command
        try:
            subprocess.run(command, check=True)
            print(f"Successfully concatenated and outputted to: {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error concatenating VCF files: {e}")

    def run_plink2(self, output_prefix:str, threads=28, memory=500000):
        """
        Runs plink2 with the specified parameters to process a VCF file and create binary PLINK files.

        Parameters:
        - vcf_file: str, path to the input VCF file (e.g., 'annotated_normalized_combined_1_22.vcf.gz').
        - output_prefix: str, prefix for the output files (e.g., '/mnt/vdb/imputed_data_luxgiant/unzipped/protocols/step3_finaldata_biallelic_only/final_luxgiant_data').
        - threads: int, number of threads to use (default is 28).
        - memory: int, amount of memory to use in MB (default is 500000).
        """

        results_dir = self.results_dir

        input_vcf = os.path.join(results_dir, "annotated_normalized_combined_1_22.vcf.gz")

        # Build the plink2 command
        command = [
            "plink2",
            "--vcf", input_vcf,
            "--snps-only", "just-acgt",  # Only keep SNPs with ACGT alleles
            "--make-bed",  # Create binary PLINK files
            "--out", os.path.join(results_dir, output_prefix),
            "--threads", str(threads),
            "--memory", str(memory)
        ]

        # Execute the command
        try:
            subprocess.run(command, check=True)
            print(f"PLINK2 command executed successfully. Output files saved with prefix: {output_prefix}")
        except subprocess.CalledProcessError as e:
            print(f"Error running PLINK2: {e}")

        pass