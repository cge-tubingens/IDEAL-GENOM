"""
Module to prepare data for the downstream analysis
"""
import os

from luxgiant_dstream.Helpers import shell_do, delete_temp_files

class PrepDS:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, dependables_path:str) -> None:

        # check if paths are set
        if input_path is None or output_path is None or dependables_path is None:
            raise ValueError("Values for input_path, output_path and dependables_path must be set upon initialization.")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.dependables = dependables_path
        self.config_dict = config_dict
        
        self.files_to_keep = []

        # create results folder
        self.results_dir = os.path.join(output_path, 'preparatory')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        pass

    def exclude_high_ld_hla(self)->dict:

        input_path       = self.input_path
        input_name       = self.input_name
        results_dir      = self.results_dir
        output_name      = self.output_name
        dependables_path = self.dependables

        maf      = self.config_dict['maf']
        geno     = self.config_dict['geno']
        mind     = self.config_dict['mind']
        hwe      = self.config_dict['hwe']
        ind_pair = self.config_dict['indep-pairwise']

        # Check type of maf
        if not isinstance(maf, float):
             raise TypeError("maf should be of type float.")

        # Check type of geno
        if not isinstance(geno, float):
            raise TypeError("geno should be of type float.")

        # Check type of mind
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check if maf is in range
        if maf < 0.05 or maf > 0.1:
            raise ValueError("maf should be between 0.05 and 0.1")
        
        # Check if geno is in range
        if geno < 0.05 or geno > 0.1:
            raise ValueError("geno should be between 0.05 and 0.1")
        
        # Check if mind is in range
        if mind < 0.1 or mind > 0.15:
            raise ValueError("mind should be between 0.1 and 0.15")
        
        # Check if hwe is in range
        if hwe < 0.00000001 or hwe > 0.001:
            raise ValueError("hwe should be between 0.00000001 and 0.001")
        
        # check existence of high LD regions file
        high_ld_regions_file = os.path.join(dependables_path, 'high-LD-regions.txt')
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError(f"File with high LD region was not found: {high_ld_regions_file}")

        step = "ld_prune"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Run plink to exclude high LD regions
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --chr 1-22 --maf {maf} --geno {geno}  --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --make-bed --out {os.path.join(results_dir, output_name+'_prunning')}"

        # LD pruning
        plink_cmd2 = f"plink2 --bfile {os.path.join(results_dir, output_name+'_prunning')} --extract {os.path.join(results_dir, output_name+'_prunning.prune.in')} --make-bed --out {os.path.join(results_dir, output_name+'_LDpruned')} --threads {max_threads}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.files_to_keep.append(output_name+'_LDpruned.bed')
        self.files_to_keep.append(output_name+'_LDpruned.bim')
        self.files_to_keep.append(output_name+'_LDpruned.fam')

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def pca_decomposition(self)->dict:

        results_dir = self.results_dir
        output_name = self.output_name

        pca = self.config_dict['pca']

        step = "pca_decomposition"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Run plink to perform PCA decomposition
        plink_cmd = f"plink --bfile {os.path.join(results_dir, output_name+'_LDpruned')} --pca {pca} --threads {max_threads} --out {os.path.join(results_dir, output_name+'_pca')}"

        shell_do(plink_cmd, log=True)

        self.files_to_keep.append(output_name+'_pca.eigenvec')

        delete_temp_files(self.files_to_keep, results_dir)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict