import gzip
import os
import shutil
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import textalloc as ta

from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase

from gwaslab.bd_download import download_file
from gwaslab.g_Log import Log
from gwaslab.util_in_get_sig import annogene

def filter_sumstats(data_df:pd.DataFrame, lead_snp:str, snp_col:str, p_col:str, pos_col:str, chr_col:str, pval_threshold:float=5e-8, radius:int=10e6, get_gene_names:bool=True, build:str='38', gtf_path:str=None) -> pd.DataFrame:
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

    if get_gene_names:

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
                df_filtered,
                id     =snp_col,
                chrom  =chr_col,
                pos    =pos_col,
                log    =Log(),
                build  =build,
                source ="refseq",
                verbose=True,
                gtf_path=gtf_path
            ).rename(columns={"GENE":"GENENAME"})   

    return variants_toanno

