"""This module provides functions to create a zoomed heatmap visualization
of SNP associations, gene annotations, and linkage disequilibrium (LD) patterns.

It includes filtering and annotating SNP data, calculating LD matrices,
and generating a three-panel plot with:
1. Association plot with SNPs colored by functional consequences
2. Gene track showing gene locations and orientations
3. LD heatmap showing correlation patterns between SNPs
"""

import os
import time
from matplotlib import cm
import logging

import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
import pandas as pd

from matplotlib.patches import FancyArrow

from ideal_genom.annotations import annotate_snp, gtf_to_all_genes
from ideal_genom.get_references import Ensembl37Fetcher, Ensembl38Fetcher

from pyensembl import Genome
from typing import Optional, Union

from itertools import cycle

from ideal_genom.api_client import VEPEnsemblRestClient

from ideal_genom.Helpers import shell_do

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", force=True)
logger = logging.getLogger(__name__)

def filter_sumstats(data_df: pd.DataFrame, lead_snp: str, snp_col: str, p_col: str, pos_col: str, chr_col: str, pval_threshold: float = 5e-8, radius: Union[float, int] = 10e6) -> pd.DataFrame:
    """Filter GWAS summary statistics based on a lead SNP, p-value threshold and genomic region.
    
    This function filters a DataFrame containing GWAS summary statistics to return variants
    within a specified genomic region around a lead SNP that meet a p-value significance threshold.
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        DataFrame containing GWAS summary statistics
    lead_snp : str
        Identifier of the lead SNP to center the region on
    snp_col : str
        Name of column containing SNP identifiers
    p_col : str
        Name of column containing p-values
    pos_col : str
        Name of column containing genomic positions
    chr_col : str
        Name of column containing chromosome numbers/identifiers
    pval_threshold : float, optional
        P-value significance threshold for filtering variants (default: 5e-8)
    radius : float, optional
        Size of region to include around lead SNP in base pairs (default: 10Mb)
    
    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame containing only variants that:

            - Are on the same chromosome as lead SNP
            - Meet p-value threshold
            - Fall within specified region around lead SNP

        Also includes calculated -log10(p-value) column

    Raises
    ------
    TypeError
        If input parameters are not of the expected types
    ValueError
        If specified columns are not found in the DataFrame
        If lead SNP is not found in the DataFrame
    
    Notes
    -----
    The function adds a 'log10p' column containing -log10 transformed p-values
    to the filtered DataFrame before returning it.
    """

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
    if not isinstance(radius, (float, int)):
        raise TypeError("radius must be a float or an integer.")
    if snp_col not in data_df.columns:
        raise ValueError(f"Column {snp_col} not found in the data frame.")
    if p_col not in data_df.columns:
        raise ValueError(f"Column {p_col} not found in the data frame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column {pos_col} not found in the data frame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column {chr_col} not found in the data frame.")
    
    lead_snp_mask = (data_df[snp_col]==lead_snp)
    if not lead_snp_mask.any():
        raise ValueError(f"Lead SNP {lead_snp} not found in the data frame.")

    lead_chr = data_df[lead_snp_mask][chr_col].values[0]
    lead_pos = data_df[lead_snp_mask][pos_col].values[0]

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

def snp_annotations(data_df: pd.DataFrame, snp_col: str, pos_col: str, chr_col: str, build: str = '38', anno_source: str = 'ensembl', gtf_path: Optional[str] = None, batch_size: int = 100, request_persec: int = 15) -> pd.DataFrame:
    """Annotate SNPs with gene names and functional consequences using Ensembl databases.
    
    This function takes a DataFrame containing SNP information and adds gene name annotations
    and functional consequence annotations using Ensembl VEP (Variant Effect Predictor) API.
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        Input DataFrame containing SNP information
    snp_col : str
        Name of column containing SNP IDs
    pos_col : str 
        Name of column containing genomic positions
    chr_col : str
        Name of column containing chromosome numbers
    build : str, optional
        Genome build version ('38', '37', or '19'), by default '38'
    anno_source : str, optional
        Source for annotations ('ensembl'), by default 'ensembl'
    gtf_path : str, optional
        Path to GTF file for annotations, by default None
    batch_size : int, optional
        Number of SNPs to process in each API request batch, by default 100
    request_persec : int, optional
        Maximum number of API requests per second, by default 15
    
    Returns
    -------
    pd.DataFrame
        Input DataFrame augmented with gene name and functional consequence annotations.
        Added columns:
        - GENENAME: Gene name from Ensembl
        - Functional_Consequence: Most severe consequence from VEP
    
    Raises
    ------
    ValueError
        If specified genome build version is not supported
        If annotation source is not supported
        If the specified columns are not found in the DataFrame
    TypeError
        If input parameters are not of the expected types
    
    Notes
    -----
    Supports genome builds 19/37 and 38 using different Ensembl REST API endpoints.
    Implements rate limiting and batch processing for API requests.
    """

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("data_df must be a pandas DataFrame.")
    if not isinstance(snp_col, str):
        raise TypeError("snp_col must be a string.")
    if not isinstance(pos_col, str):
        raise TypeError("pos_col must be a string.")
    if not isinstance(chr_col, str):
        raise TypeError("chr_col must be a string.")
    if not isinstance(build, str):
        raise TypeError("build must be a string.")
    if not isinstance(anno_source, str):
        raise TypeError("anno_source must be a string.")
    if not isinstance(batch_size, int):
        raise TypeError("batch_size must be an integer.")
    if not isinstance(request_persec, int):
        raise TypeError("request_persec must be an integer.")
    if build not in ['19', '37', '38']:
        raise ValueError(f"Unsupported build version: {build}. Supported versions are '19', '37', and '38'.")
    if anno_source != 'ensembl':
        raise ValueError(f"Unsupported annotation source: {anno_source}. Currently only 'ensembl' is supported.")
    if gtf_path is not None and not isinstance(gtf_path, str):
        raise TypeError("gtf_path must be a string or None.")
    if snp_col not in data_df.columns:
        raise ValueError(f"Column {snp_col} not found in the data frame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column {pos_col} not found in the data frame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column {chr_col} not found in the data frame.")

    variants_toanno = annotate_snp(
            insumstats=data_df,
            chrom  =chr_col,
            pos    =pos_col,
            build  =build,
            source =anno_source,
            gtf_path=gtf_path # type: ignore
        ).rename(columns={"GENE":"GENENAME"})
    
    
    logger.info(f" - Annotated SNPs with gene names.")
    logger.info(f'variants_toanno shape: {variants_toanno.shape}')
    logger.info(f'variants_toanno columns: {variants_toanno.columns.tolist()}')
    
    if build == '38':

        logger.info(" -Assigning functional consequences using Ensembl VEP for protein coding genes for build 38")
    
        # initialize VEP client
        vep_client = VEPEnsemblRestClient(server='https://rest.ensembl.org', reqs_per_sec=request_persec)

        logger.info(" - Using Ensembl VEP API to annotate SNPs with functional consequences.")

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

        logger.info(" -Assigning functional consequences using Ensembl VEP for protein coding genes for build 37")

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
    else:
        raise ValueError(f"Unsupported build version: {build}. Supported versions are '19', '37', and '38'.")

    logger.info(" - Finished annotating SNPs with functional consequences.")
    logger.info(f'df_vep shape: {df_vep.shape}')
    logger.info(f'df_vep columns: {df_vep.columns.tolist()}')

    variants_toanno = variants_toanno.merge(df_vep, on=snp_col, how='left')
    
    return variants_toanno

def get_gene_information(genes:list, gtf_path: Optional[str] = None, build: str = "38", anno_source: str = 'ensembl') -> pd.DataFrame:
    """Retrieves genomic information for a list of genes using Ensembl annotation.
    
    This function fetches start position, end position, strand, and length information
    for each gene in the provided list using either Ensembl GRCh37 or GRCh38 annotations.
    
    Parameters
    ----------
    genes : list
        List of gene IDs (Ensembl format)
    gtf_path : str, optional
        Path to a custom GTF file. If None, will download and use Ensembl GTF.
    build : str, default "38"
        Human genome build version. Supported values: "19", "37", "38"
        Note: "19" and "37" are equivalent.
    anno_source : str, default "ensembl"
        Source of genome annotations. Currently only supports "ensembl"
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing gene information with columns:
            - gene: gene ID
            - start: gene start position
            - end: gene end position
            - strand: gene strand
            - length: gene length
    
    Raises
    ------
    ValueError
        If unsupported build version or annotation source is provided
    FileNotFoundError
        If provided GTF file path does not exist
    TypeError
        If provided GTF path is not a string
    
    Notes
    -----
    When gtf_path is None, the function will automatically download and process
    the appropriate Ensembl GTF file based on the specified build version.
    The function uses the Ensembl Python API to fetch gene information.
    """

    if anno_source == "ensembl":

        if build=="19" or build=="37":

            print(" -Assigning Gene name using ensembl_hg37_gtf for protein coding genes")
  
            if gtf_path is None:

                nsmbl37 = Ensembl37Fetcher()

                nsmbl37.get_latest_release()
                nsmbl37.download_latest()
                nsmbl37.unzip_latest()
                nsmbl37.get_all_genes()

                gtf_path = nsmbl37.gtf_file

            else:
                if not os.path.isfile(gtf_path):
                    raise FileNotFoundError(f"Provided GTF path {gtf_path} does not exist.")
                if not isinstance(gtf_path, str):
                    raise TypeError(f"Provided GTF path {gtf_path} is not a string.")
                print(" -Using user-provided gtf:{}".format(gtf_path))
                
                gtf_path = gtf_to_all_genes(gtf_path)

            gtf_db_path = gtf_path[:-2]+"db" # type: ignore
            
            try:
                data = Genome(
                    reference_name='GRCh38',
                    annotation_name='Ensembl',
                    gtf_path_or_url=gtf_path
                )
            except Exception as e:
                logger.error(f"Error initializing Genome with GTF path {gtf_path}: {e}")
                raise

            if os.path.isfile(gtf_db_path) is False:

                data.index()
        
        elif build=="38":

            print(" -Assigning Gene name using ensembl_hg38_gtf for protein coding genes")

            if gtf_path is None:

                nsmbl38 = Ensembl38Fetcher() 

                nsmbl38.get_latest_release()
                nsmbl38.download_latest()
                nsmbl38.unzip_latest()
                nsmbl38.get_all_genes()

                gtf_path = nsmbl38.all_genes_path

            else:
                print(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_genes(gtf_path)
            
            gtf_db_path = gtf_path[:-2]+"db"

            try:
                data = Genome(
                    reference_name='GRCh38',
                    annotation_name='Ensembl',
                    gtf_path_or_url=gtf_path
                )
            except Exception as e:
                logger.error(f"Error initializing Genome with GTF path {gtf_path}: {e}")
                raise
        else:
            raise ValueError(f"Unsupported build version: {build}. Supported versions are '19', '37', and '38'.")

    else:
        raise ValueError(f"Unsupported annotation source: {anno_source}. Supported sources are 'ensembl'.")
    
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

def get_ld_matrix(data_df: pd.DataFrame, snp_col: str, pos_col: str, bfile_folder: str, bfile_name: str, output_path: str) -> dict:
    """Calculate LD matrix using PLINK for a set of SNPs.
    
    This function takes a DataFrame containing SNP information and calculates the LD (Linkage Disequilibrium)
    matrix using PLINK. The SNPs are first sorted by position, and then PLINK is used to compute
    pairwise r2 values between SNPs.
    
    Parameters
    ----------
    data_df : pandas.DataFrame 
        DataFrame containing SNP information
    snp_col : str 
        Name of the column containing SNP IDs
    pos_col : str 
        Name of the column containing SNP positions
    bfile_folder : str 
        Path to the folder containing PLINK binary files
    bfile_name : str 
        Base name of the PLINK binary files (without extensions)
    output_path : str 
        Path where output files will be saved
    
    Returns
    -------
    dict 
        Dictionary containing:
            - 'pass' : bool 
                True if process completed successfully
            - 'step' : str 
                Name of the processing step ('get_ld_matrix')
            - 'output' : dict 
                Dictionary with output file paths
    
    Raises
    ------
    FileNotFoundError 
        If any required files or directories are not found
    TypeError
        If input parameters are not of the expected types
    ValueError 
        If specified columns are not found in the DataFrame
    """

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
    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("data_df must be a pandas DataFrame.")
    if not isinstance(snp_col, str):
        raise TypeError("snp_col must be a string.")
    if not isinstance(pos_col, str):
        raise TypeError("pos_col must be a string.")
    
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

def get_zoomed_data(data_df: pd.DataFrame, lead_snp: str, snp_col: str, p_col: str, pos_col: str, chr_col: str, output_folder: str, pval_threshold: float = 5e-6, radius: Union[float, int] = 1e6, build: str = '38', anno_source: str = 'ensembl', gtf_path: Optional[str] = None, batch_size: int = 100, request_persec: int = 15) -> pd.DataFrame:
    """Filter and annotate SNP data around a lead SNP within a specified radius.
    
    This function filters significant SNPs in a region around a lead SNP and annotates them with gene
    names and functional consequences. The position values are scaled to Megabase pairs (Mbp).
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        Input DataFrame containing SNP data
    lead_snp : str
        Identifier of the lead SNP to center the region around
    snp_col : str
        Name of the column containing SNP identifiers
    p_col : str
        Name of the column containing p-values
    pos_col : str
        Name of the column containing position information
    chr_col : str
        Name of the column containing chromosome information
    output_folder : str
        Path to the output folder (must exist)
    pval_threshold : float, optional
        P-value threshold for significance filtering (default: 5e-6)
    radius : Union[float, int], optional
        Radius around the lead SNP in base pairs (default: 1e6)
    build : str, optional
        Genome build version ('38' or '37') (default: '38')
    anno_source : str, optional
        Source for annotations ('ensembl' or other supported sources) (default: 'ensembl')
    gtf_path : str, optional
        Path to GTF file for annotations (default: None)
    batch_size : int, optional
        Number of SNPs to process in each batch (default: 100)
    request_persec : int, optional
        Number of API requests per second allowed (default: 15)
    
    Returns
    -------
    pandas.DataFrame
        Filtered and annotated DataFrame with added Mbp column and removed duplicates
    
    Raises
    ------
    TypeError
        If input parameters are not of the expected types
    FileNotFoundError
        If `output_folder` does not exist
    ValueError
        If no significant SNPs are found in the specified region
    
    Notes
    -----
    The function removes duplicate SNPs, keeping the first occurrence only.
    Position values are converted to Megabase pairs in the output DataFrame.
    """

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
    if not isinstance(radius, float) and not isinstance(radius, int):
        raise TypeError("radius must be a float or an integer.")

    if os.path.isdir(output_folder) is not True:
        raise FileNotFoundError(f"Folder {output_folder} not found.")
    
    logger.info(f" - Filtering data for lead SNP: {lead_snp} in region of radius {radius} around the lead SNP.")
    
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

    logger.info(f" - Filtered data shape: {filtered_df.shape}")
    logger.info(f" - Filtered data columns: {filtered_df.columns.tolist()}")

    if filtered_df.empty:
        raise ValueError("No significant SNPs found in the specified region.")
    
    # annotate the SNPs with gene names and functional consequences
    annotated = snp_annotations(
        data_df       =filtered_df, 
        snp_col       =snp_col, 
        chr_col       =chr_col, 
        pos_col       =pos_col,
        build         =build,
        batch_size    =batch_size,
        request_persec=request_persec,
        anno_source   =anno_source,
        gtf_path      =gtf_path
    )

    # scale the position to Mbp
    annotated['Mbp'] = annotated['POS'] / 1e6

    annotated = annotated.drop_duplicates(keep='first').reset_index(drop=True)

    annotated = annotated.drop_duplicates(subset=[snp_col], keep='first').reset_index(drop=True)

    return annotated

def draw_zoomed_heatmap(data_df: pd.DataFrame, lead_snp: str, snp_col: str, p_col: str, pos_col: str, chr_col: str, output_folder: str, bfile_folder: str, bfile_name: str, pval_threshold: float = 5e-6, radius: Union[int, float] = 1e6, build: str ='38', gtf_path: Optional[str] = None, anno_source: str = "ensembl", batch_size: int = 100, effect_dict: dict = dict(), extension: str = 'pdf', request_persec: int = 15) -> bool:
    """Creates a zoomed heatmap visualization around a lead SNP showing LD patterns and gene annotations.
    
    This function generates a three-panel plot:
    1. Association plot with SNPs colored by functional consequences
    2. Gene track showing gene locations and orientations
    3. LD heatmap showing correlation patterns between SNPs
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        Input DataFrame containing GWAS summary statistics
    lead_snp : str
        Identifier of the lead SNP to center the plot around
    snp_col : str
        Column name containing SNP identifiers
    p_col : str 
        Column name containing p-values
    pos_col : str
        Column name containing genomic positions
    chr_col : str
        Column name containing chromosome numbers
    output_folder : str
        Path to save output files
    bfile_folder : str
        Folder containing PLINK binary files
    bfile_name : str
        Base name of PLINK binary files (without extensions)
    pval_threshold : float, optional
        P-value threshold for significance, default 5e-6
    radius : Union[int, float], optional
        Distance in base pairs to plot around lead SNP, default 1e6
    build : str, optional
        Genome build version, default '38'
    gtf_path : str, optional
        Path to custom GTF file, default None
    anno_source : str, optional
        Source for gene annotations ('ensembl' or 'refseq'), default 'ensembl'
    batch_size : int, optional
        Batch size for API requests, default 100
    effect_dict : dict, optional
        Dictionary mapping functional effects to display names, default empty dict
    extension : str, optional
        File extension for output plot, default 'pdf'
    request_persec : int, optional
        Number of API requests per second allowed, default 15
    
    Returns
    -------
    bool
        True if plot was generated successfully

    Raises
    ------
    TypeError
        If input parameters are not of the expected types
    FileNotFoundError
        If specified folders or PLINK binary files do not exist
    ValueError
        If the required columns are not found in the DataFrame
    
    Notes
    -----
    Required input DataFrame must contain columns for SNP IDs, p-values, positions and chromosomes.
    PLINK binary files (.bed, .bim, .fam) must exist in specified folder.
    Generates and saves a zoomed heatmap plot in the specified output folder.
    """

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
    if not isinstance(output_folder, str):
        raise TypeError("output_folder must be a string.")
    if not isinstance(bfile_folder, str):
        raise TypeError("bfile_folder must be a string.")
    if not isinstance(bfile_name, str):
        raise TypeError("bfile_name must be a string.")
    if not isinstance(pval_threshold, float):
        raise TypeError("pval_threshold must be a float.")
    if not isinstance(radius, (float, int)):
        raise TypeError("radius must be a float or an integer.")
    if not isinstance(build, str):
        raise TypeError("build must be a string.")
    if not isinstance(gtf_path, (str, type(None))):
        raise TypeError("gtf_path must be a string or None.")
    if not isinstance(anno_source, str):
        raise TypeError("anno_source must be a string.")
    if not isinstance(batch_size, int):
        raise TypeError("batch_size must be an integer.")
    if not isinstance(effect_dict, dict):
        raise TypeError("effect_dict must be a dictionary.")
    if not isinstance(extension, str):
        raise TypeError("extension must be a string.")
    if not isinstance(request_persec, int):
        raise TypeError("request_persec must be an integer.")
    
    if not os.path.isdir(output_folder):
        raise FileNotFoundError(f"Output folder {output_folder} does not exist.")
    if not os.path.isdir(bfile_folder):
        raise FileNotFoundError(f"PLINK binary folder {bfile_folder} does not exist.")
    if not os.path.exists(os.path.join(bfile_folder, f"{bfile_name}.bim")):
        raise FileNotFoundError(f"PLINK binary file {bfile_name}.bim not found in {bfile_folder}.")
    if not os.path.exists(os.path.join(bfile_folder, f"{bfile_name}.bed")):
        raise FileNotFoundError(f"PLINK binary file {bfile_name}.bed not found in {bfile_folder}.")
    if not os.path.exists(os.path.join(bfile_folder, f"{bfile_name}.fam")):
        raise FileNotFoundError(f"PLINK binary file {bfile_name}.fam not found in {bfile_folder}.")
    
    if snp_col not in data_df.columns:
        raise ValueError(f"Column {snp_col} not found in the data frame.")
    if p_col not in data_df.columns:
        raise ValueError(f"Column {p_col} not found in the data frame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column {pos_col} not found in the data frame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column {chr_col} not found in the data frame.")

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
        batch_size    =batch_size,
        request_persec=request_persec,
        build         =build,
        gtf_path      =gtf_path,
    )

    logger.info(f" - Annotated data shape: {annotated.shape}")

    annotated['GENENAME'] = annotated['GENENAME'].apply(lambda x: x.split(',')[0])

    effects = annotated['Functional_Consequence'].value_counts(dropna=False).reset_index()

    region = (annotated['Mbp'].min() - 0.05, annotated['Mbp'].max() + 0.05)

    logger.info(f" - Region for zoom plot: {region}")
    logger.info(f" - Lead SNP: {lead_snp}")
    logger.info(f" - Number of SNPs in the region: {annotated.shape[0]}")

    genes = get_gene_information(
        genes   =annotated['GENENAME'].unique().tolist(),
        gtf_path=None,
        build   ='38',
        anno_source=anno_source
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
        sep      =r'\s+',
        header   =None,
        index_col=None,
        engine   ='python'
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
    lead = pd.DataFrame()
    
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

    if lead.empty:
        raise ValueError(f"Lead SNP {lead_snp} not found in the data.")
   
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
            arrow = FancyArrow(start, y, length, 0, width=0.001, head_width=0.03, head_length=0.01, color=color) # type: ignore
            ax2.add_patch(arrow)
            ax2.text(start + 0.5 * length, y + 0.05, symbol, ha='center', size=9) # type: ignore
        elif strand == '-':
            arrow_neg = FancyArrow(end, y, -length, 0, width=0.001, head_width=0.03, head_length=0.01, color=color) # type: ignore
            ax2.add_patch(arrow_neg)
            ax2.text(start + 0.5 * length, y + 0.05, symbol, ha='center', size=9) # type: ignore

    ax2.set_ylim(0, 1.2)
    ax2.set_xlim(region)
    ax2.axis('off')

    base = ax3.transData # to rotate triangle
    rotation = transforms.Affine2D().rotate_deg(180+90+45)
    cmap = cm.get_cmap('Reds')
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
