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
    
    vep_client = VEPEnsemblRestClient()

    # Example list of IDs for the POST request
    snps = variants_toanno[snp_col].to_list()

    response = vep_client.post_vep_request(snps)

    if response:
        df_vep = pd.DataFrame({
            snp_col: [ res['id'] for res in response ],
            'consequence': [ res['most_severe_consequence'] for res in response ]
        })

        variants_toanno = variants_toanno.merge(df_vep, on=snp_col, how='left')
    else:
        print("Failed to get response.")
    
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
        data = Genome(
            reference_name='GRCh37',
            annotation_name='Refseq',
            gtf_path_or_url=gtf_path
        )

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

    sorted_data = data_df.sort_values(by=[pos_col], ascending=True).reset_index(drop=True)
    sorted_data[[snp_col]].to_csv(
        os.path.join(bfile_folder, f"{bfile_name}-snplist.txt"),
        header=False,
        index=False,
        sep='\t'
    )

    # plink command
    plink_cmd = f"plink --bfile {os.path.join(bfile_folder, bfile_name)} --extract {os.path.join(bfile_folder, f'{bfile_name}-snplist.txt')} --r2 square  --out {os.path.join(bfile_folder, 'matrix-ld')}"

    # execute PLINK command
    shell_do(plink_cmd, log=True)



    # report
    process_complete = True

    outfiles_dict = {
        'plink_out': bfile_folder
    }

    out_dict = {
        'pass': process_complete,
        # 'step': step,
        'output': outfiles_dict
    }
        
    return out_dict
