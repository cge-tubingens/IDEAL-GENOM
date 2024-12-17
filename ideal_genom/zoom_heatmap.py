import gzip
import os
import shutil
import warnings
import time
import matplotlib

import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
import pandas as pd
import seaborn as sns
import textalloc as ta

from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrow

from gwaslab.bd_download import download_file
from gwaslab.bd_common_data import gtf_to_all_gene
from gwaslab.g_Log import Log
from gwaslab.util_in_get_sig import annogene

from pyensembl import Genome

from itertools import cycle

from ideal_genom.api_client import VEPEnsemblRestClient

from ideal_genom.Helpers import shell_do

def filter_sumstats(data_df:pd.DataFrame, lead_snp:str, snp_col:str, p_col:str, pos_col:str, chr_col:str, pval_threshold:float=5e-8, radius:int=10e6) -> pd.DataFrame:
    """
    Filter the summary statistics data frame to only include SNPs with p-value less than the threshold
    and the lead SNP
    """

    lead_chr = data_df[data_df[snp_col]==lead_snp][chr_col].values[0]
    lead_pos = data_df[data_df[snp_col]==lead_snp][pos_col].values[0]

    mask_chr = (data_df[chr_col] == lead_chr)
    mask_pval= (data_df[p_col] <= pval_threshold)

    df_filtered = data_df[mask_chr & mask_pval].reset_index(drop=True)

    df_filtered['log10p'] = -np.log10(df_filtered[p_col])
    
    upper_bound = lead_pos + radius
    lower_bound = lead_pos - radius

    mask_upper = (df_filtered[pos_col] <= upper_bound)
    mask_lower = (df_filtered[pos_col] >= lower_bound)

    df_filtered = df_filtered[mask_upper & mask_lower].reset_index(drop=True)   

    return df_filtered

def snp_annotations(data_df:pd.DataFrame, snp_col:str, pos_col:str, chr_col:str, build:str='38', gtf_path:str=None, batch_size:int=100, request_persec:int=15) -> pd.DataFrame:
    
    if gtf_path is None:
        gtf_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
        path_to_gz = os.path.join(os.path.abspath('..'), 'GCF_000001405.40_GRCh38.p14_genomic.gtf.gz')
        path_to_gtf= os.path.join(os.path.abspath('..'), 'GCF_000001405.40_GRCh38.p14_genomic.gtf')
        
        if os.path.exists(path_to_gz) is not True or os.path.exists(path_to_gtf) is not True:
            download_file(gtf_url, path_to_gz)
            
            with gzip.open(path_to_gz, 'rb') as f_in:
                 with open(path_to_gtf, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        gtf_path = path_to_gtf

    variants_toanno = annogene(
            data_df,
            id     =snp_col,
            chrom  =chr_col,
            pos    =pos_col,
            log    =Log(),
            build  =build,
            source ="refseq",
            verbose=True,
            gtf_path=gtf_path
        ).rename(columns={"GENE":"GENENAME"})
    
    if build == '38':
    
        # initialize VEP client
        vep_client = VEPEnsemblRestClient(server='https://rest.ensembl.org', reqs_per_sec=request_persec)

        # list of IDs for the POST request
        snps = variants_toanno[snp_col].to_list()

        # empty DataFrame to store the VEP results
        df_vep = pd.DataFrame()

        # iterate through the list of IDs in batches
        for i in range(0, len(snps), batch_size):

            batch = snps[i:min(i + batch_size, len(snps))]

            response = vep_client.post_vep_request(batch)

            if response:
                batch_df = pd.DataFrame({
                    snp_col: [res['id'] for res in response],
                    'Functional_Consequence': [res['most_severe_consequence'] for res in response]
                })
                df_vep = pd.concat([df_vep, batch_df], ignore_index=True)
            else:
                print("Failed to get response.")

            time.sleep(5)

    elif build == '19' or build == '37':

        # initialize VEP client
        vep_client = VEPEnsemblRestClient(server='https://grch37.rest.ensembl.org', reqs_per_sec=15)

        # list of IDs for the POST request
        snps = variants_toanno[snp_col].to_list()

        # empty DataFrame to store the VEP results
        df_vep = pd.DataFrame()

        # iterate through the list of IDs in batches
        for i in range(0, len(snps), batch_size):

            batch = snps[i:min(i + batch_size, len(snps))]

            response = vep_client.post_vep_request(batch)

            if response:
                batch_df = pd.DataFrame({
                    snp_col: [res['id'] for res in response],
                    'Functional_Consequence': [res['most_severe_consequence'] for res in response]
                })
                df_vep = pd.concat([df_vep, batch_df], ignore_index=True)
            else:
                print("Failed to get response.")

            time.sleep(5)

    variants_toanno = variants_toanno.merge(df_vep, on=snp_col, how='left')
    
    return variants_toanno.drop(columns=['LOCATION'], inplace=False)

def get_gene_information(genes:list, gtf_path:str=None, build:str=38)->pd.DataFrame:

    if gtf_path is None:
        gtf_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
        path_to_gz = os.path.join(os.path.abspath('..'), 'GCF_000001405.40_GRCh38.p14_genomic.gtf.gz')
        path_to_gtf= os.path.join(os.path.abspath('..'), 'GCF_000001405.40_GRCh38.p14_genomic.gtf')
        
        if os.path.exists(path_to_gz) is not True or os.path.exists(path_to_gtf) is not True:
            download_file(gtf_url, path_to_gz)
            
            with gzip.open(path_to_gz, 'rb') as f_in:
                 with open(path_to_gtf, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        gtf_path = path_to_gtf

    gtf_path = gtf_to_all_gene(gtf_path, log=Log())

    if build == '38':
        data = Genome(
            reference_name='GRCh38',
            annotation_name='Refseq',
            gtf_path_or_url=gtf_path
        )
    elif build == '19':
        # data = Genome(
        #     reference_name='GRCh37',
        #     annotation_name='Refseq',
        #     gtf_path_or_url=gtf_path
        # )
        warnings.warn("Build 19 not supported. Using build 38 instead.")

    gene_info = {
        'gene':genes,
        'start':[],
        'end':[],
        'strand':[],
        'length':[]
    }

    for gene in gene_info['gene']:
        try:
            gene_info['start'].append(data.gene_by_id(gene).start)
            gene_info['end'].append(data.gene_by_id(gene).end)
            gene_info['strand'].append(data.gene_by_id(gene).strand)
            gene_info['length'].append(data.gene_by_id(gene).length)
        except:
            gene_info['start'].append(None)
            gene_info['end'].append(None)
            gene_info['strand'].append(None)
            gene_info['length'].append(None)

    return pd.DataFrame(gene_info)

def get_ld_matrix(data_df:pd.DataFrame, snp_col:str, pos_col:str, bfile_folder:str, bfile_name:str, output_path:str)->dict:

    if os.path.exists(output_path) is not True:
        raise FileNotFoundError(f"File {output_path} not found.")
    if os.path.isdir(bfile_folder) is not True:
        raise FileNotFoundError(f"File {bfile_folder} not found.")
    if os.path.exists(os.path.join(bfile_folder, f"{bfile_name}.bim")) is not True:
        raise FileNotFoundError(f"File {bfile_name}.bim not found in {bfile_folder}.")
    if os.path.exists(os.path.join(bfile_folder, f"{bfile_name}.bed")) is not True:
        raise FileNotFoundError(f"File {bfile_name}.bed not found in {bfile_folder}.")
    if os.path.exists(os.path.join(bfile_folder, f"{bfile_name}.fam")) is not True:
        raise FileNotFoundError(f"File {bfile_name}.fam not found in {bfile_folder}.")
    
    if snp_col not in data_df.columns:
        raise ValueError(f"Column {snp_col} not found in the data frame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column {pos_col} not found in the data frame.")
    
    step='get_ld_matrix'

    sorted_data = data_df.sort_values(by=[pos_col], ascending=True).reset_index(drop=True)
    sorted_data[[snp_col]].to_csv(
        os.path.join(bfile_folder, f"{bfile_name}-snplist.txt"),
        header=False,
        index=False,
        sep='\t'
    )

    # plink command
    plink_cmd = f"plink --bfile {os.path.join(bfile_folder, bfile_name)} --extract {os.path.join(bfile_folder, f'{bfile_name}-snplist.txt')} --r2 square  --out {os.path.join(output_path, 'matrix-ld')}"

    # execute PLINK command
    shell_do(plink_cmd, log=True)

    # report
    process_complete = True

    outfiles_dict = {
        'plink_out': output_path
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'output': outfiles_dict
    }
        
    return out_dict

def get_zoomed_data(data_df:pd.DataFrame, lead_snp:str, snp_col:str, p_col:str, pos_col:str, chr_col:str, output_folder:str,pval_threshold:float=5e-6, radius:int=1e6, build='38', batch_size:int=100, request_persec:int=15)->pd.DataFrame:

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("data_df must be a pandas DataFrame.")
    if not isinstance(lead_snp, str):
        raise TypeError("lead_snp must be a string.")
    if not isinstance(snp_col, str):
        raise TypeError("snp_col must be a string.")
    if not isinstance(p_col, str):
        raise TypeError("p_col must be a string.")
    if not isinstance(pos_col, str):
        raise TypeError("pos_col must be a string.")
    if not isinstance(chr_col, str):
        raise TypeError("chr_col must be a string.")
    if not isinstance(pval_threshold, float):
        raise TypeError("pval_threshold must be a float.")
    if not isinstance(radius, int):
        raise TypeError("radius must be an integer.")
    
    if os.path.isdir(output_folder) is not True:
        raise FileNotFoundError(f"Folder {output_folder} not found.")
    
    # filter significant SNPs in the specified region
    filtered_df = filter_sumstats(
        data_df       =data_df, 
        lead_snp      =lead_snp, 
        snp_col       =snp_col, 
        p_col         =p_col, 
        pos_col       =pos_col, 
        chr_col       =chr_col, 
        pval_threshold=pval_threshold, 
        radius        =radius
    )

    if filtered_df.empty:
        raise ValueError("No significant SNPs found in the specified region.")
    
    # annotate the SNPs with gene names and functional consequences
    annotated = snp_annotations(
        data_df=filtered_df, 
        snp_col=snp_col, 
        chr_col=chr_col, 
        pos_col=pos_col,
        build=build,
        batch_size=batch_size,
        request_persec=request_persec
    )

    # scale the position to Mbp
    annotated['Mbp'] = annotated['POS'] / 1e6

    return annotated.drop_duplicates(keep='first').reset_index(drop=True)

def draw_zoomed_heatmap(data_df:pd.DataFrame, lead_snp:str, snp_col:str, p_col:str, pos_col:str, chr_col:str, output_folder:str, pval_threshold:float=5e-6, radius:int=1e6, build='38', gtf_path:str=None, batch_size:int=100, bfile_folder:str=None, bfile_name:str=None, effect_dict:dict={}, extension:str='jpeg', request_persec:int=15)->None:

    annotated = get_zoomed_data(
        data_df       =data_df,
        lead_snp      =lead_snp, 
        snp_col       =snp_col, 
        p_col         =p_col, 
        pos_col       =pos_col, 
        chr_col       =chr_col,
        output_folder =output_folder,
        pval_threshold=pval_threshold, 
        radius        =radius,
        batch_size=batch_size,
        request_persec=request_persec
    )

    annotated['GENENAME'] = annotated['GENENAME'].apply(lambda x: x.split(',')[0])

    effects = annotated['Functional_Consequence'].value_counts(dropna=False).reset_index()

    region = (annotated['Mbp'].min() - 0.05, annotated['Mbp'].max() + 0.05)

    print('\n')
    print(effects)
    print('\n')

    genes = get_gene_information(
        genes=annotated['GENENAME'].unique().tolist(),
        gtf_path=None,
        build='38'
    )
    genes['start_esc'] = genes['start']/1e6
    genes['end_esc']   = genes['end']/1e6

    genes['start_esc'] = genes['start_esc'].apply(lambda x: max(x, region[0]))
    genes['end_esc'] = genes['end_esc'].apply(lambda x: min(x, region[1]))
    
    genes['length_esc']= genes['end_esc'] - genes['start_esc']

    annotated = annotated.merge(genes, left_on='GENENAME', right_on='gene', how='left')
    annotated.to_csv(os.path.join(output_folder, f'zoom_plot_data_for_{lead_snp}.csv'), index=False, sep='\t')

    get_ld_matrix(
        data_df     =annotated,
        snp_col     =snp_col,
        pos_col     =pos_col,
        bfile_folder=bfile_folder,
        bfile_name  =bfile_name,
        output_path =output_folder,
    )

    df_LD = pd.read_csv(
        os.path.join(output_folder, 'matrix-ld.ld'),
        sep=r'\s+',
        header=None,
        index_col=None,
        engine='python'
    )
    ld = df_LD.values

    # plot the heatmap

    N=ld.shape[0]
    ld = np.tril(ld, k=0)
    ldm = np.ma.masked_where(ld==0, ld)

    plt.figure(figsize=(10, 10))

    # Define the overall grid size (9 rows, 1 column)
    ax1 = plt.subplot2grid((9, 1), (0, 0), rowspan=4)  # Top plot (4 rows)
    ax2 = plt.subplot2grid((9, 1), (4, 0), rowspan=1)  # Middle plot (1 row)
    ax3 = plt.subplot2grid((9, 1), (5, 0), rowspan=4)  # Bottom plot (4 rows)

    # Define custom colors
    colors = [
        "#1f77b4",  # Blue
        "#ff7f0e",  # Orange
        "#2ca02c",  # Green
        "#9467bd",  # Purple
        "#8c564b",  # Brown
        "#17becf",  # Cyan
        "#bcbd22",  # Yellow
        "#e377c2"   # Pink
    ]

    # Plot for ax1

    annotated['Hue'] = 'other'

    if len(effect_dict) == 0:

        main_effects = effects['Functional_Consequence'].values[:4]

        for effect in main_effects:

            annotated.loc[annotated['Functional_Consequence'] == effect, 'Hue'] = effect
    else:

        for effect in effect_dict.keys():

            annotated.loc[annotated['Functional_Consequence'] == effect, 'Hue'] = effect_dict[effect]

    # plot other SNPs
    others = annotated[annotated['Hue'] == 'other']
    ax1.scatter(others['Mbp'], others['log10p'], s=15, color='grey', label='', edgecolors='none')

    # plot main effects
    main = annotated[(annotated['Hue'] != 'other') & (annotated['Hue'] != 'Lead Variant')]

    annotated_lead = False
    
    for k, effect in enumerate(main['Hue'].unique()):
        subset = main[main['Hue'] == effect]
        ax1.scatter(subset['Mbp'], subset['log10p'], s=15, color=colors[k], label=effect)
        if lead_snp in subset[snp_col].values: # plot lead SNP
            annotated_lead = True
            lead = subset[subset[snp_col] == lead_snp]
            ax1.scatter(lead['Mbp'], lead['log10p'], s=30, color=colors[k], label='Lead SNP', marker='d')
    
    if annotated_lead is False:
        lead = annotated[annotated[snp_col] == lead_snp]
        ax1.scatter(lead['Mbp'], lead['log10p'], s=30, color='red', label='Lead SNP'+f'{lead["Functional_Consequence"].values[0]}', marker='d')
   
    chr = lead[chr_col].values[0]

    ax1.set_xlim(region)
    ax1.xaxis.set_ticks_position('top')
    ax1.legend(loc='best')
    ax1.set_title(f"Zoom of {lead_snp}", fontsize=12, loc='left')
    ax1.set_ylabel('log10(P)', fontsize=12)
    ax1.xaxis.set_label_position('top')
    ax1.set_xlabel(f'Position on Chr {chr} [Mb]', fontsize=12)

    # Plot for ax2
    ys = cycle([0.1, 0.4, 0.7, 1])

    for i in genes.index:
        symbol, strand = genes.loc[i, 'gene'], genes.loc[i, 'strand']
        start, end, length = genes.loc[i, 'start_esc'], genes.loc[i, 'end_esc'], genes.loc[i, 'length_esc']
        y = next(ys)

        if symbol == lead['GENENAME'].values[0]:
            color = 'red'
        else:
            color = 'black'

        if strand == '+':
            arrow = FancyArrow(start, y, length, 0, width=0.001, head_width=0.03, head_length=0.01, color=color)
            ax2.add_patch(arrow)
            ax2.text(start + 0.5 * length, y + 0.05, symbol, ha='center', size=9)
        elif strand == '-':
            arrow_neg = FancyArrow(end, y, -length, 0, width=0.001, head_width=0.03, head_length=0.01, color=color)
            ax2.add_patch(arrow_neg)
            ax2.text(start + 0.5 * length, y + 0.05, symbol, ha='center', size=9)

    ax2.set_ylim(0, 1.2)
    ax2.set_xlim(region)
    ax2.axis('off')

    base = ax3.transData # to rotate triangle
    rotation = transforms.Affine2D().rotate_deg(180+90+45)
    cmap=matplotlib.cm.Reds
    im=ax3.imshow(ldm, cmap=cmap, transform=rotation+base,aspect='auto')
    ax3.set_xlim([2, 1.41*N])
    ax3.set_ylim([1*N, 2])
    ax3.axis('off')

    # Add a colorbar as the legend
    cbar = plt.colorbar(im, ax=ax3, orientation='horizontal', fraction=0.05, pad=0.2)
    cbar.set_label('LD Value', fontsize=10)  # Adjust label as needed

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f'Zoom for {lead_snp}.{extension}'), dpi=500)
    plt.show()

    return True
