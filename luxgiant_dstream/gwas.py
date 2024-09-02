import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from adjustText import adjust_text
from luxgiant_dstream.Helpers import shell_do
from luxgiant_dstream.annotate_tools import get_variant_context

class GWAS:

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

        # create results folder
        self.results_dir = os.path.join(output_path, 'gwas_analysis')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create figures folder
        self.plots_dir = os.path.join(output_path, 'plots')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

        pass

    def exclude_high_ld_hla(self)->dict:

        results_dir = self.results_dir
        output_name = self.output_name
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
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, output_name)} --chr 1-22 --maf {maf} --geno {geno}  --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --make-bed --out {os.path.join(results_dir, output_name+'_prunning')}"

        # LD pruning
        plink_cmd2 = f"plink2 --bfile {os.path.join(results_dir, output_name+'_prunning')} --extract {os.path.join(results_dir, output_name+'_prunning.prune.in')} --make-bed --out {os.path.join(results_dir, output_name+'_LDpruned')} --threads {max_threads}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            print(cmd)
            shell_do(cmd, log=True)

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

    def fixed_model_association_analysis(self)->dict:

        output_name= self.output_name
        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir

        maf = self.config_dict['maf']
        mind= self.config_dict['mind']
        hwe = self.config_dict['hwe']
        ci  = self.config_dict['ci']

        step = "association_analysis"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2  # use all available cores
        else:
            max_threads = 10

        # Run plink2 to perform association analysis

        plink2_cmd = f"plink2 --bfile {os.path.join(input_path, input_name)} --adjust --ci {ci} --maf {maf} --mind {mind} --hwe {hwe} --covar {os.path.join(results_dir, output_name+'_pca.eigenvec')} --glm hide-covar omit-ref sex cols=+a1freq,+beta --out {os.path.join(results_dir, output_name+'_glm1')} --threads {max_threads}"

        shell_do(plink2_cmd, log=True)

        df = pd.read_csv(os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), sep="\t")
        rename = {
            '#CHROM'          : 'CHR',	
            'POS'             : 'POS',	
            'ID'              : 'SNP',
            'REF'             : 'A2',	
            'ALT'             : 'ALT',	
            'PROVISIONAL_REF?': 'PROVISIONAL_REF',	
            'A1'              : 'A1',	
            'OMITTED'         : 'OMITTED',	
            'A1_FREQ'         : 'freq',	
            'FIRTH?'          : 'FIRTH',	
            'TEST'            : 'TEST',	
            'OBS_CT'          : 'N',	
            'BETA'            : 'b',	
            'SE'              : 'se',	
            'L95'             : 'L95',	
            'U95'             : 'U95',	
            'Z_STAT'          : 'Z_STAT',	
            'P'               : 'p',	
            'ERRCODE'         : 'ERRCODE'
        }
        df = df.rename(columns=rename)
        df.to_csv(os.path.join(results_dir, output_name+'_glm1.PHENO1.glm.logistic.hybrid'), sep="\t", index=False)

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
    
    def random_model_association_analysis(self)->dict:

        pass
     
    def get_top_hits(self)->dict:

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name

        maf = self.config_dict['maf']

        step = "get_top_hits"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # load results of association analysis
        df = pd.read_csv(os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), sep="\t")

        # prepare .ma file
        df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

        df.to_csv(os.path.join(results_dir, 'cojo_file.ma'), sep="\t", index=False)

        # gcta command
        gcta_cmd = f"gcta64 --bfile {os.path.join(input_path, input_name)} --maf {maf} --cojo-slct --cojo-file {os.path.join(results_dir, 'cojo_file.ma')}   --out {os.path.join(results_dir, 'cojo_file')} --thread-num {max_threads}"

        shell_do(gcta_cmd, log=True)

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
    
    def annotate_hits(self)->dict:

        import time

        results_dir = self.results_dir

        step = "annotate_hits"

        # load the data
        df_hits = pd.read_csv(os.path.join(results_dir, 'cojo_file.jma'), sep="\t")

        df_hits = df_hits[['Chr', 'SNP', 'bp']].copy()

        for k in range(df_hits.shape[0]):
            # get variant context
            chr = df_hits.loc[k, 'Chr']
            pos = df_hits.loc[k, 'bp']

            context = get_variant_context(chr, pos)

            if context is None:
                context = 'NA'
            df_hits.loc[k, 'GENE'] = context[0]

            time.sleep(1.5)

        df_hits = df_hits[['SNP', 'GENE']].copy()

        df_hits.to_csv(os.path.join(results_dir, 'snps_annotated.csv'), sep="\t", index=False)

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

    def manhattan_plot(self)->dict:

        results_dir = self.results_dir
        output_name = self.output_name
        plots_dir   = self.plots_dir
        annotate    = self.config_dict['annotate']

        step = "manhattan_plot"

        # load the data
        df_gwas = pd.read_csv(os.path.join(results_dir, output_name+'_glm.PHENO1.glm.logistic.hybrid'), sep="\t")
        df_gene = pd.read_csv(os.path.join(results_dir, "snps_annotated.csv"), sep="\t")

        # keep columns of interest
        df = df_gwas[['SNP', 'CHR', 'p']].copy()
        df['log10P'] = -np.log10(df['p'])

        # sort values by chromosome
        df = df.sort_values('CHR')

        # to get colors by chromosome
        df['ind'] = range(len(df))
        df_grouped = df.groupby(('CHR'))

        # Subset dataframe to highlight specific SNPs
        snps = df_gene['SNP'].to_list()
        genes = df_gene['GENE'].to_list()

        highlighted_snps = df[df['SNP'].isin(snps)]  # Filter for the SNPs of interest

        # Set the figure size
        fig = plt.figure(figsize=(50, 25))
        ax = fig.add_subplot(111)

        # Colors for alternating chromosomes
        colors = ['grey', 'skyblue']

        # Variables to store labels and positions for the x-axis
        x_labels = []
        x_labels_pos = []

        # Loop over each group by chromosome
        for num, (name, group) in enumerate(df_grouped):
            # Plot each chromosome with alternating colors
            ax.scatter(group['ind'], group['log10P'], color=colors[num % len(colors)], s=10)

            # Add chromosome label and its position
            x_labels.append(name)
            middle_pos = (group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2
            x_labels_pos.append(middle_pos)

            # Set the x-axis labels and positions
            ax.set_xticks(x_labels_pos)
            ax.set_xticklabels(x_labels, fontsize=20)  # Adjust fontsize as needed

            # Plot highlighted SNPs with a different color (red) and larger point size
            ax.scatter(highlighted_snps['ind'], highlighted_snps['log10P'], color='red', s=50, label='Highlighted SNPs')

        if annotate:
            # Add gene names to highlighted SNPs
            texts = []  # A list to store text annotations for adjustment
            for i, row in highlighted_snps.iterrows():
                gene = genes[snps.index(row['SNP'])]  # Get corresponding gene name
                # Add text label to the SNP
                text = ax.text(row['ind'], row['log10P'], gene, fontsize=12, ha='right', va='bottom', color='black',
                           bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))
                texts.append(text)

            # Adjust the text to prevent overlaps using adjustText
            adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'))

        # Set axis limits
        ax.set_xlim([0, len(df)])
        ax.set_ylim([0, max(df['log10P']) + 1])

        # Add horizontal line for genome-wide significance threshold (-log10(5e-8))
        ax.axhline(y=-np.log10(5e-8), color='blue', linestyle='--', linewidth=1)

        # Add labels and title
        ax.set_xlabel('Chromosome', fontsize=24)  # Adjust fontsize as needed
        ax.set_ylabel('-log10(P-value)', fontsize=24)  # Adjust fontsize as needed
        ax.set_title('Manhattan Plot', fontsize=30)  # Add a title if needed

        # Customize the grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        # Display the plot
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'manhattan_plot.png'))

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

    def qq_plot(self)->dict:
    
        import gc

        results_dir      = self.results_dir
        plots_dir        = self.plots_dir

        step = "qq_plot"

        # load the data
        df_gwas = pd.read_csv(os.path.join(results_dir, "gwas_metafile.txt"), sep="\t")

        pvalues = df_gwas['P'].values

        grp = None
        n = len(pvalues)

        log_p = -np.log10(pvalues)
        exp_x = -np.log10((stats.rankdata(pvalues, method='ordinal') - 0.5) / n)

        if grp is not None:  
            # Create a DataFrame with rounded pvalues and exp_x, and grp
            thin = pd.DataFrame({
                'pvalues': np.round(log_p, 3),
                'exp_x': np.round(exp_x, 3),
                'grp': grp
            }).drop_duplicates()
            # Assign the updated group after thinning
            grp = thin['grp'].values
        else:
            # Create a DataFrame with rounded pvalues and exp_x
            thin = pd.DataFrame({
                'pvalues': np.round(log_p, 3),
                'exp_x': np.round(exp_x, 3)
            }).drop_duplicates()
        # Update pvalues and exp_x after thinning
        log_p = thin['pvalues'].values
        exp_x = thin['exp_x'].values

        gc.collect()

        axis_range =  [float(min(log_p.min(), exp_x.min()))-0.5, float(max(log_p.max(), exp_x.max()))+1]

        fig, ax = plt.subplots(figsize=(10,10))

        # Calculate the confidence intervals if draw_conf is True
        def plot_confidence_interval(n, conf_points=1500, conf_col="lightgray", conf_alpha=0.05):
            conf_points = min(conf_points, n - 1)
            mpts = np.zeros((conf_points * 2, 2))

            for i in range(1, conf_points + 1):
                x = -np.log10((i - 0.5) / n)
                y_upper = -np.log10(stats.beta.ppf(1 - conf_alpha / 2, i, n - i))
                y_lower = -np.log10(stats.beta.ppf(conf_alpha / 2, i, n - i))

                mpts[i - 1, 0] = x
                mpts[i - 1, 1] = y_upper
                mpts[conf_points * 2 - i, 0] = x
                mpts[conf_points * 2 - i, 1] = y_lower

            # Plot the confidence interval as a filled polygon
            plt.fill(mpts[:, 0], mpts[:, 1], color=conf_col)

            pass

        plot_confidence_interval(n, conf_points=1500, conf_col="lightgray", conf_alpha=0.05)

        if grp is not None:
            unique_groups = np.unique(grp)
            for group in unique_groups:
                group_mask = (grp == group)
                ax.scatter(exp_x[group_mask], log_p[group_mask], label=f'Group {group}', marker='o')
        else:
            ax.scatter(exp_x, log_p, marker='o')

        # Line y = x
        ax.plot(axis_range, axis_range, color='red', linestyle='--', lw=1)

        ax.set_xlim(axis_range)
        ax.set_ylim(axis_range)

        ax.set_xlabel('Expected (-log10 p-value)')
        ax.set_ylabel('Observed (-log10 p-value)')
        ax.set_aspect('equal')

        ax.grid(True)

        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'qq_plot.png'))

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
