"""This module provides functions to annotate genomic variants with gene information and effects.

Includes functions to:
- Find the closest gene to a given SNP position.
- Map chromosome numbers to identifiers.
- Convert GTF files to a format containing all genes.
- Annotate SNPs with gene names using Ensembl or RefSeq databases.
- Prepare genome data from GTF files.
- Annotate variants with their effects relying on Ensembl VEP.
"""
import os
import logging

import pandas as pd

from pyensembl import Genome

from gtfparse import read_gtf

from ideal_genom.get_references import Ensembl38Fetcher, Ensembl37Fetcher, RefSeqFetcher

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", force=True)
logger = logging.getLogger(__name__)

def get_closest_gene(x, data: Genome, chrom: str = "CHR", pos: str = "POS", max_iter: int = 20000, step: int = 50, source: str = "ensembl", build: str="38") -> tuple:
    
    """Find the closest gene to a given position in the genome.

    This function searches for the closest gene to a specified SNP position in the genome.
    It checks the position in the specified chromosome and returns the distance to the closest gene
    along with the gene name(s). If no gene is found within the specified distance, it returns "intergenic".
    
    Parameters
    ----------
    x : 
        SNP information.
    data : pyensembl.Genome
        An instance of the Genome class containing gene annotations.
    chrom : str, optional 
        The key in the dictionary `x` that corresponds to the chromosome. Default is "CHR".
    pos : str, optional
        The key in the dictionary `x` that corresponds to the position. Default is "POS".
    max_iter : int, optional 
        The maximum number of iterations to search for a gene. Default is 20000.
    step : int, optional
        The step size for each iteration when searching for a gene. Default is 50.
    source : str, optional 
        The source of the gene annotations, either "ensembl" or "refseq". Default is "ensembl".
    build : str, optional 
        The genome build version, used when source is "refseq". Default is "38".
    
    Returns
    -------
    tuple: 
        A tuple containing the distance to the closest gene and the gene name(s). If no gene is found, returns the distance and "intergenic".

    Raises
    ------
    TypeError
        If `data` is not an instance of Genome, or if `chrom` or `pos` are not strings, or if `max_iter` or `step` are not integers.
    ValueError
        If `source` is not "ensembl" or "refseq", or if `build` is not "37" or "38".
    """

    if not isinstance(data, Genome):
        raise TypeError("Data must be an instance of Genome from pyensembl.")
    if not isinstance(chrom, str):
        raise TypeError("Chromosome key must be a string.")
    if not isinstance(pos, str):
        raise TypeError("Position key must be a string.")
    if not isinstance(max_iter, int):
        raise TypeError("max_iter must be an integer.")
    if not isinstance(step, int):
        raise TypeError("step must be an integer.")
    
    if source not in ["ensembl", "refseq"]:
        raise ValueError("Source must be either 'ensembl' or 'refseq'.")
    if build not in ["37", "38"]:
        raise ValueError("Build must be either '37' or '38'.")

    def clean_empty(gene):
        # remove empty elements
        return [x for x in gene if x!=""]    
    
    # from GWASlab
        
    #convert 23,24,25 back to X,Y,MT for EnsemblRelease query
    if source=="ensembl":
            
        contig = get_number_to_chr()[x[chrom]]

    elif source=="refseq":
            
        contig = get_number_to_chr()[x[chrom]]
        contig = get_chr_to_NC(build=build)[contig]
        # for refseq , gene names are stored as gene_id, using gene_ids_at_locus instead
            
        data.gene_names_at_locus = data.gene_ids_at_locus
    
    position = int(x[pos])
        # query
    gene = data.gene_names_at_locus(contig=contig, position=position)

    if len(clean_empty(gene))==0:
        # if not in any gene
        i=0
        while i<=max_iter:
            # using distance to check upstram and downstream region
            distance = i*step
            # upstream
            gene_u = data.gene_names_at_locus(contig=contig, position=position-distance)
            
            # downstream
            gene_d = data.gene_names_at_locus(contig=contig, position=position+distance)
            
            if len(clean_empty(gene_u))>0 and len(clean_empty(gene_d))>0:
                # if found gene uptream and downstream at the same time 
                # go back to last step
                distance = (i-1)*step
                for j in range(0,step,1):
                    # use small step to finemap                        
                    gene_u = data.gene_names_at_locus(contig=contig, position=position-distance-j)
                    gene_d = data.gene_names_at_locus(contig=contig, position=position+distance+j)
                    if len(clean_empty(gene_u))>0:
                        return -distance-j,",".join(gene_u).strip(",")
                    elif len(clean_empty(gene_d))>0:
                        return distance+j,",".join(gene_d).strip(",")
            elif len(clean_empty(gene_u))>0:                    
                # if found gene uptream
                distance = (i-1)*step
                for j in range(0,step,1):
                    gene_u2 = data.gene_names_at_locus(contig=contig, position=position-distance-j)
                    if len(clean_empty(gene_u2))>0:
                        return -distance-j,",".join(gene_u).strip(",")
            elif len(clean_empty(gene_d))>0:
                # if found gene downstream
                distance = (i-1)*step
                for j in range(0,step,1):
                    gene_d2 = data.gene_names_at_locus(contig=contig, position=position+distance+j)
                    if len(clean_empty(gene_d2))>0:
                        return distance+j,",".join(gene_d).strip(",")
            i+=1
            # increase i by 1
        return distance, "intergenic"
    else:
        return 0,",".join(gene).strip(",")
        
def get_number_to_chr(in_chr: bool = False, xymt: list = ["X","Y","MT"], xymt_num: list = [23,24,25], prefix: str = "") -> dict:
    
    """Creates a dictionary mapping chromosome numbers to chromosome identifiers.
    
    This function generates a mapping between chromosome numbers (as keys) and 
    chromosome identifiers (as values), with special handling for sex chromosomes
    and mitochondrial chromosome.
    
    Parameters
    ----------
    in_chr : bool, default=False
        If True, dictionary keys will be strings; if False, keys will be integers.
    xymt : list, default=["X","Y","MT"]
        List of string identifiers for the X, Y, and mitochondrial chromosomes.
    xymt_num : list, default=[23,24,25]
        List of numeric identifiers corresponding to X, Y, and MT chromosomes.
    prefix : str, default=""
        String prefix to add to all chromosome identifiers.
        
    Returns
    -------
    dict
        A dictionary mapping chromosome numbers to chromosome identifiers.
        For autosomal chromosomes (1-199), maps to prefix+number.
        For sex and mitochondrial chromosomes, maps to prefix+X/Y/MT.

    Raises
    ------
    TypeError
        If `in_chr` is not a boolean, `xymt` or `xymt_num` are not lists, or `prefix` is not a string.
    
    Examples
    --------
    >>> get_number_to_chr()
    {1: '1', 2: '2', ..., 23: 'X', 24: 'Y', 25: 'MT', ...}
    
    >>> get_number_to_chr(in_chr=True, prefix="chr")
    {'1': 'chr1', '2': 'chr2', ..., '23': 'chrX', '24': 'chrY', '25': 'chrMT', ...}
    """

    if not isinstance(in_chr, bool):
        raise TypeError("in_chr must be a boolean value.")
    if not isinstance(xymt, list) or not isinstance(xymt_num, list):
        raise TypeError("xymt and xymt_num must be lists.")
    if not isinstance(prefix, str):
        raise TypeError("prefix must be a string.")

    # from GWASlab
    if in_chr is True:
        dic= {str(i):prefix+str(i) for i in range(1,200)}
        dic[str(xymt_num[0])]=prefix+xymt[0]
        dic[str(xymt_num[1])]=prefix+xymt[1]
        dic[str(xymt_num[2])]=prefix+xymt[2]
    else:
        dic= {i:prefix+str(i) for i in range(1,200)}
        dic[xymt_num[0]]=prefix+xymt[0]
        dic[xymt_num[1]]=prefix+xymt[1]
        dic[xymt_num[2]]=prefix+xymt[2]
    return dic

def get_chr_to_NC(build: str = '38', inverse: bool = False) -> dict:
    """Returns a dictionary mapping between chromosome names and NCBI NC identifiers.

    This function provides a mapping between chromosome names (like "1", "X", "MT") 
    and their corresponding NCBI RefSeq accession numbers (like "NC_000001.10") 
    for different human genome builds.

    Parameters
    ----------
    build : str, optional
        The genome build version. Accepted values are "19", "37", or "38".
        Note that builds "19" and "37" return the same mapping.
    inverse : bool, optional 
        If True, returns an inverted dictionary where NC identifiers are keys 
        and chromosome names are values. Defaults to False.

    Returns
    -------
    dict 
        A dictionary mapping chromosome names to NC identifiers (if inverse=False)
        or NC identifiers to chromosome names (if inverse=True).

    Raises
    ------
    TypeError
        If `build` is not a string or `inverse` is not a boolean.
    ValueError
        If `build` is not one of "19", "37", or "38".

    References
    ----------
        https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
    """

    if not isinstance(build, str):
        raise TypeError("Build must be a string.")
    if build not in ["19", "37", "38"]:
        raise ValueError("Build must be one of '19', '37', or '38'.")
    if not isinstance(inverse, bool):
        raise TypeError("Inverse must be a boolean value.")

    # from GWASlab
    if build =="19" or build=="37":
        dic={
        "1":"NC_000001.10",
        "2":"NC_000002.11",
        "3":"NC_000003.11",
        "4":"NC_000004.11",
        "5":"NC_000005.9",
        "6":"NC_000006.11",
        "7":"NC_000007.13",
        "8":"NC_000008.10",
        "9":"NC_000009.11",
        "10":"NC_000010.10",
        "11":"NC_000011.9",
        "12":"NC_000012.11",
        "13":"NC_000013.10",
        "14":"NC_000014.8",
        "15":"NC_000015.9",
        "16":"NC_000016.9",
        "17":"NC_000017.10",
        "18":"NC_000018.9",
        "19":"NC_000019.9",
        "20":"NC_000020.10",
        "21":"NC_000021.8",
        "22":"NC_000022.10",
        "X":"NC_000023.10",
        "Y":"NC_000024.9",
        "MT":"NC_012920.1"}
    elif build=="38":
        dic={
        "1":"NC_000001.11",
        "2":"NC_000002.12",
        "3":"NC_000003.12",
        "4":"NC_000004.12",
        "5":"NC_000005.10",
        "6":"NC_000006.12",
        "7":"NC_000007.14",
        "8":"NC_000008.11",
        "9":"NC_000009.12",
        "10":"NC_000010.11",
        "11":"NC_000011.10",
        "12":"NC_000012.12",
        "13":"NC_000013.11",
        "14":"NC_000014.9",
        "15":"NC_000015.10",
        "16":"NC_000016.10",
        "17":"NC_000017.11",
        "18":"NC_000018.10",
        "19":"NC_000019.10",
        "20":"NC_000020.11",
        "21":"NC_000021.9",
        "22":"NC_000022.11",
        "X":"NC_000023.11",
        "Y":"NC_000024.1",
        "MT":"NC_012920.1"
        }
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic

def gtf_to_all_genes(gtfpath: str) -> str:
    """Extract all gene records from a GTF file and save them to a new file.
    
    This function reads a GTF file, extracts all gene records, and saves them to a new file
    with the suffix '_all_genes.gtf.gz'. If the output file already exists, it will be returned
    without regenerating it.

    Parameters
    ----------
    gtfpath : str 
        Path to the input GTF file.
    
    Returns
    -------
    str 
        Path to the output file containing all gene records.

    Raises
    ------
    TypeError
        If `gtfpath` is not a string.
    
    Notes
    -----
    The function uses the `read_gtf` function for initial parsing and pandas for extraction.
    The function assumes the GTF file has a standard format with gene_id attributes.
    """

    if not isinstance(gtfpath, str):
        raise TypeError("GTF path must be a string.")
    
    all_gene_path = gtfpath[:-6]+"all_genes.gtf.gz"
    
    # if not existing, extract protein coding records and output to a new file
    if not os.path.isfile(all_gene_path):
        
        # get gene list
        logger.info(f" - Extracting genes from {gtfpath}")
        
        gtf = read_gtf(gtfpath,usecols=["feature", "gene_biotype", "gene_id", "gene_name"])

        gene_list = gtf.loc[gtf["feature"]=="gene", "gene_id"].values  # type: ignore
        
        logger.info(f" - Loaded {gene_list} genes.")
        
        # extract entry using csv
        gtf_raw = pd.read_csv(gtfpath,sep="\t",header=None,comment="#",dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[ gtf_raw["_gene_id"].isin(gene_list) ,:]
        gtf_raw = gtf_raw.drop("_gene_id",axis=1)
        
        logger.info(f" - Extracted records are saved to : {all_gene_path} ")

        gtf_raw.to_csv(all_gene_path, header=False, index=False, sep='\t')

    return all_gene_path

def annotate_snp(insumstats: pd.DataFrame, gtf_path: str, chrom: str = "CHR", pos: str = "POS", build: str = "38", source: str = "ensembl") -> pd.DataFrame:
    """Annotate SNPs with nearest gene name(s) using either Ensembl or RefSeq databases.
    
    This function takes a DataFrame containing SNP data and annotates each variant
    with information about the nearest gene(s) based on genomic coordinates.

    Parameters
    ----------
    insumstats : pandas.DataFrame 
        DataFrame containing SNP data with chromosome and position information.
    gtf_path : str 
        Path to the GTF (Gene Transfer Format) file for gene annotations.
    chrom : str, optional
        Column name in the DataFrame that contains chromosome information. Defaults to "CHR".
    pos : str, optional 
        Column name in the DataFrame that contains position information. Defaults to "POS".
    build : str, optional 
        Genome build version. Must be one of "19", "37", or "38". Defaults to "38".
    source : str, optional 
        Source for gene annotation. Must be either "ensembl" or "refseq". Defaults to "ensembl".

    Returns
    -------
    pandas.DataFrame 
        A copy of the input DataFrame with additional gene annotation columns.

    Raises
    ------
    TypeError 
        If input is not a pandas DataFrame or if GTF path is not a string.
    ValueError 
        If required columns are missing in the input DataFrame or if build/source parameters are invalid.
    """

    if not isinstance(insumstats, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    if chrom not in insumstats.columns or pos not in insumstats.columns:
        raise ValueError(f"Input DataFrame must contain columns '{chrom}' and '{pos}'.")
    if build not in ["19", "37", "38"]:
        raise ValueError("Build must be one of '19', '37', or '38'.")
    if source not in ["ensembl", "refseq"]:
        raise ValueError("Source must be either 'ensembl' or 'refseq'.")
    if gtf_path:
        if not isinstance(gtf_path, str):
            raise TypeError("GTF path must be a string or None.")
    else:
        gtf_path = ''
    
    output = insumstats.copy()
    
    is_gtf_path = os.path.isfile(gtf_path)

    logger.info("Starting to annotate variants with nearest gene name(s)...")
    logger.info(f" -Using {build} as genome build")
    logger.info(f"is_gtf_path set to {is_gtf_path}")

    if source == "ensembl":
        logger.info(f" -Using ensembl as source for gene annotation")
        output = annotate_with_ensembl(output, chrom, pos, build, gtf_path, is_gtf_path)
    elif source == "refseq":
        logger.info(f" -Using refseq as source for gene annotation")
        output = annotate_with_refseq(output, chrom, pos, build, gtf_path, is_gtf_path)

    logger.info("Finished annotating variants with nearest gene name(s) successfully!")
    return output

def annotate_with_ensembl(output: pd.DataFrame, chrom: str, pos: str, build: str, gtf_path: str, is_gtf_path: bool) -> pd.DataFrame:
    """Annotate variants with gene information from Ensembl database.
    
    This function adds gene annotations to a DataFrame containing variant information
    by looking up the genomic coordinates in Ensembl data. It adds 'LOCATION' and 'GENE' 
    columns to the input DataFrame.
    
    Parameters
    ----------
    output : pandas.DataFrame
        DataFrame containing variant information with chromosome and position columns.
    chrom : str
        Name of the column in the DataFrame that contains chromosome information.
    pos : str
        Name of the column in the DataFrame that contains position information.
    build : str
        Genome build version to use. Must be one of '19', '37', or '38'.
        Note that '19' and '37' are treated as equivalent (GRCh37).
    gtf_path : str
        Path to GTF file with gene annotations or None to use default paths.
        If None, the appropriate GTF file will be downloaded or used from cache.
    is_gtf_path : bool
        If True, gtf_path is treated as a direct path to a GTF file.
        If False, gtf_path is treated as a directory where the GTF file should be downloaded.
    
    Returns
    -------
    pandas.DataFrame
        The input DataFrame with additional 'LOCATION' and 'GENE' columns containing
        gene annotations from Ensembl.
    
    Raises
    ------
    TypeError
        If output is not a pandas DataFrame or if gtf_path is not a string (when provided).
    ValueError
        If the required columns are not in the DataFrame or if the build is invalid.
    
    Notes
    -----
    The function supports both GRCh37 (build '19' or '37') and GRCh38 (build '38')
    and will download the appropriate annotation files if not already available.
    """

    if not isinstance(output, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    if chrom not in output.columns or pos not in output.columns:
        raise ValueError(f"Input DataFrame must contain columns '{chrom}' and '{pos}'.")
    if build not in ["19", "37", "38"]:
        raise ValueError("Build must be one of '19', '37', or '38'.")
    if gtf_path is not None and not isinstance(gtf_path, str):
        raise TypeError("GTF path must be a string.")
    
    if build in ["19", "37"]:
        logger.info(" -Assigning Gene name using Ensembl GRCh37 for protein coding genes")
        gtf_path = prepare_gtf_path(gtf_path, is_gtf_path, source="ensembl", build="37")
        data = prepare_genome(gtf_path, "GRCh37", "Ensembl")

    elif build == "38":
        logger.info(" -Assigning Gene name using Ensembl GRCh38 for protein coding genes")
        gtf_path = prepare_gtf_path(gtf_path, is_gtf_path, source="ensembl", build="38")
        data = prepare_genome(gtf_path, "GRCh38", "Ensembl")

    output.loc[:, ["LOCATION", "GENE"]] = annotate_variants(output, data, chrom, pos, "ensembl")
    return output

def annotate_with_refseq(output: pd.DataFrame, chrom: str, pos: str, build: str, gtf_path: str, is_gtf_path: bool) -> pd.DataFrame:
    """Annotate genomic variants with RefSeq gene information.
    
    This function adds gene and location annotations to genomic variants using NCBI RefSeq data.
    It processes the input DataFrame and adds two new columns: 'LOCATION' and 'GENE'.
    
    Parameters
    ----------
    output : pandas.DataFrame 
        DataFrame containing variant information to annotate.
    chrom : str 
        Column name in DataFrame that contains chromosome information.
    pos : str 
        Column name in DataFrame that contains position information.
    build : str 
        Genome build version. Must be one of '19', '37', or '38'.
    gtf_path : str 
        Path to the GTF file. If None, a default path will be used.
    is_gtf_path : bool 
        If True, gtf_path is treated as a direct file path.
        If False, gtf_path is treated as a directory.
    
    Returns
    -------
    pandas.DataFrame 
        The input DataFrame with added 'LOCATION' and 'GENE' columns.
    
    Raises
    ------
    TypeError 
        If output is not a pandas DataFrame or if gtf_path is provided but not a string.
    ValueError 
        If required columns are missing from output or if build is invalid.
    
    Notes
    -----
        - For builds '19' and '37', GRCh37 RefSeq annotations are used.
        - For build '38', GRCh38 RefSeq annotations are used.
        - Only protein-coding genes are considered for annotation.
    """

    if not isinstance(output, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    if chrom not in output.columns or pos not in output.columns:
        raise ValueError(f"Input DataFrame must contain columns '{chrom}' and '{pos}'.")
    if build not in ["19", "37", "38"]:
        raise ValueError("Build must be one of '19', '37', or '38'.")
    if gtf_path is not None and not isinstance(gtf_path, str):
        raise TypeError("GTF path must be a string.")
    
    if build in ["19", "37"]:
        logger.info(" -Assigning Gene name using NCBI refseq latest GRCh37 for protein coding genes")
        gtf_path = prepare_gtf_path(gtf_path, is_gtf_path, source="refseq", build="37")
        data = prepare_genome(gtf_path, "GRCh37", "Refseq")
    elif build == "38":
        logger.info(" -Assigning Gene name using NCBI refseq latest GRCh38 for protein coding genes")
        gtf_path = prepare_gtf_path(gtf_path, is_gtf_path, source="refseq", build="38")
        data = prepare_genome(gtf_path, "GRCh38", "Refseq")
    output.loc[:, ["LOCATION", "GENE"]] = annotate_variants(output, data, chrom, pos, "refseq", build)
    return output

def prepare_gtf_path(gtf_path: str, is_gtf_path: bool, source: str, build: str) -> str:
    """Prepares the path to a GTF (Gene Transfer Format) file for annotation purposes.
    
    This function either uses a user-provided GTF file or downloads one from the
    specified source (Ensembl or RefSeq) for the given genome build (GRCh37 or GRCh38).
    If a download is required, it fetches the latest release, unzips it, and processes
    it to extract all genes.
    
    Parameters
    ----------
    gtf_path : str 
        Path to an existing GTF file, or None if one should be downloaded.
    is_gtf_path : bool 
        Flag indicating whether the provided gtf_path is valid.
    source : str 
        Source database for GTF file ('ensembl' or 'refseq').
    build : str 
        Genome build version ('37' or '38').
    
    Returns
    -------
    str 
        Path to the prepared GTF file with all genes.
    
    Notes
    -----
        - If gtf_path is None or is_gtf_path is False, a new GTF file will be downloaded.
        - If a valid gtf_path is provided, it will be processed using gtf_to_all_genes.
        - The function logs the actions being performed.
    """

    if gtf_path is None or not is_gtf_path:
        if source == "ensembl" and build == '37':
            logger.info(" -Downloading Ensembl GRCh37 GTF file")
            fetcher = Ensembl37Fetcher()
        elif source == "ensembl" and build == '38':
            logger.info(" -Downloading Ensembl GRCh38 GTF file")
            fetcher = Ensembl38Fetcher()
        elif source == "refseq" and build == '37':
            logger.info(" -Downloading RefSeq GRCh37 GTF file")
            fetcher = RefSeqFetcher(build="37")
        elif source == "refseq" and build == '38':
            logger.info(" -Downloading RefSeq GRCh38 GTF file")
            fetcher = RefSeqFetcher(build="38")

        fetcher.get_latest_release()
        fetcher.download_latest()
        fetcher.unzip_latest()
        fetcher.get_all_genes()
        gtf_path = fetcher.all_genes_path
        
    else:
        logger.info(f" -Using user-provided gtf:{gtf_path}")
        gtf_path = gtf_to_all_genes(gtf_path)
    return gtf_path

def prepare_genome(gtf_path: str, reference_name: str, annotation_name: str) -> Genome:
    """Prepare a genome annotation by loading or creating a database from a GTF file.

    This function creates a Genome object from a GTF file and ensures that the
    corresponding database is indexed for efficient access.

    Parameters
    ----------
    gtf_path : str
        Path to the GTF (Gene Transfer Format) file
    reference_name : str
        Name of the reference genome
    annotation_name : str
        Name of the annotation

    Returns
    -------
    pyensemble.Genome
        A Genome object initialized with the provided reference and annotation

    Notes
    -----
    If the database file (with extension .db) doesn't exist, this function
    will create it by calling the index() method on the Genome object.
    """
    gtf_db_path = gtf_path[:-2] + "db"
    data = Genome(reference_name=reference_name, annotation_name=annotation_name, gtf_path_or_url=gtf_path)
    if not os.path.isfile(gtf_db_path):
        data.index()
    return data

def annotate_variants(output: pd.DataFrame, data: Genome, chrom: str, pos: str, source: str, build: str='38') -> pd.DataFrame:
    """
    Annotate variants with their closest genes.

    This function processes a DataFrame containing genomic variants and enriches it 
    with gene annotation information by finding the closest gene for each variant.

    Parameters
    ----------
    output : pandas.DataFrame
        DataFrame containing variant information to be annotated.
    data : Genome
        Genome object containing reference data for annotation.
    chrom : str
        Column name in the output DataFrame that contains chromosome information.
    pos : str
        Column name in the output DataFrame that contains position information.
    source : str
        Source of the gene annotation data.
    build : str, default='38'
        Genome build version (default is GRCh38).

    Returns
    -------
    pandas.DataFrame
        DataFrame containing gene annotation information for each variant.

    Notes
    -----
    This function applies the get_closest_gene function to each row in the input DataFrame
    and returns the results as a DataFrame with the same index as the input.
    """
    return pd.DataFrame(
        list(output.apply(lambda x: get_closest_gene(x, data=data, chrom=chrom, pos=pos, source=source, build=build), axis=1)),
        index=output.index
    ).values # type: ignore
