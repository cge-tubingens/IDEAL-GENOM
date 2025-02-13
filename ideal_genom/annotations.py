import os
import logging

import pandas as pd

from pyensembl import Genome

from gtfparse import read_gtf

from ideal_genom.get_references import Ensembl38Fetcher, Ensembl37Fetcher, RefSeqFetcher

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def get_closest_gene(x, data: Genome, chrom: str = "CHR", pos: str = "POS", max_iter: int = 20000, step: int = 50, source: str = "ensembl", build: str="38"):
    
    """
    Find the closest gene to a given position in the genome.
    
    Parameters:
    -----------
    x : 
        SNP information.
    data (pyensembl.Genome): 
        An instance of the Genome class containing gene annotations.
    chrom (str, optional): 
        The key in the dictionary `x` that corresponds to the chromosome. Default is "CHR".
    pos (str, optional): 
        The key in the dictionary `x` that corresponds to the position. Default is "POS".
    max_iter (int, optional): 
        The maximum number of iterations to search for a gene. Default is 20000.
    step (int, optional): 
        The step size for each iteration when searching for a gene. Default is 50.
    source (str, optional): 
        The source of the gene annotations, either "ensembl" or "refseq". Default is "ensembl".
    build (str, optional): 
        The genome build version, used when source is "refseq". Default is "38".
    
    Returns:
    -------
    tuple: 
        A tuple containing the distance to the closest gene and the gene name(s). If no gene is found, returns the distance and "intergenic".
    """

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
        
def get_number_to_chr(in_chr=False,xymt=["X","Y","MT"],xymt_num=[23,24,25],prefix=""):

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

def get_chr_to_NC(build: str, inverse=False):
    #https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13

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

def gtf_to_all_genes(gtfpath: str = None):
    
    all_gene_path = gtfpath[:-6]+"all_genes.gtf.gz"
    
    # if not existing, extract protein coding records and output to a new file
    if not os.path.isfile(all_gene_path):
        
        # get gene list
        logger.info(f" - Extracting genes from {gtfpath}")
        
        gtf = read_gtf(gtfpath,usecols=["feature", "gene_biotype", "gene_id", "gene_name"])

        gene_list = gtf.loc[gtf["feature"]=="gene", "gene_id"].values
        
        logger.info(f" - Loaded {gene_list} genes.")
        
        # extract entry using csv
        gtf_raw = pd.read_csv(gtfpath,sep="\t",header=None,comment="#",dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[ gtf_raw["_gene_id"].isin(gene_list) ,:]
        gtf_raw = gtf_raw.drop("_gene_id",axis=1)
        
        logger.info(f" - Extracted records are saved to : {all_gene_path} ")

        gtf_raw.to_csv(all_gene_path, header=None, index=None, sep="\t")

    return all_gene_path

def annotate_snp(insumstats: pd.DataFrame, chrom: str = "CHR", pos: str = "POS", build: str = "38", source: str = "ensembl", gtf_path: str = None):
    
    output = insumstats.copy()
    
    if source == "ensembl":

        if build=="19" or build=="37":

            logger.info(" -Assigning Gene name using Ensembl GRCh37 for protein coding genes")
  
            if gtf_path is None:

                nsmbl37 = Ensembl37()

                nsmbl37.get_latest_release()
                nsmbl37.download_latest()
                nsmbl37.unzip_latest()
                nsmbl37.get_all_genes()

                gtf_path = nsmbl37.protein_coding_path

            else:
                logger.info(f" -Using user-provided gtf:{gtf_path}")
                
                gtf_path = gtf_to_all_genes(gtf_path)

            gtf_db_path = gtf_path[:-2]+"db"
            
            data = Genome(
                reference_name='GRCh37',
                annotation_name='Ensembl',
                gtf_path_or_url=gtf_path
            )

            if os.path.isfile(gtf_db_path) is False:

                data.index()
            
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x: get_closest_gene(x,data=data,chrom=chrom,pos=pos,source=source), axis=1)), 
                index=output.index).values
        
        elif build=="38":

            logger.info(" -Assigning Gene name using Ensembl GRCh38 for protein coding genes")

            if gtf_path is None:

                nsmbl38 = Ensembl38()

                nsmbl38.get_latest_release()
                nsmbl38.download_latest()
                nsmbl38.unzip_latest()
                nsmbl38.get_all_genes()

                gtf_path = nsmbl38.all_genes_path

            else:
                logger.info(f" -Using user-provided gtf:{gtf_path}")
                gtf_path = gtf_to_all_genes(gtf_path)
            
            gtf_db_path = gtf_path[:-2]+"db"

            data = Genome(
                reference_name='GRCh38',
                annotation_name='Ensembl',
                gtf_path_or_url=gtf_path
            )
            
            if os.path.isfile(gtf_db_path) is False:
                data.index()

            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:get_closest_gene(x,data=data,chrom=chrom,pos=pos,source=source), axis=1)), 
                index=output.index).values
    
    if source == "refseq":

        if build=="19" or build=="37":
        
            logger.info(" -Assigning Gene name using NCBI refseq latest GRCh37 for protein coding genes")
            
            if gtf_path is None or not is_gtf_path:

                refseq37 = RefSeqFetcher(build = "37")

                refseq37.get_latest_release()
                refseq37.download_latest()
                refseq37.unzip_latest()
                refseq37.get_all_genes()

                gtf_path = refseq37.all_genes_path

            else:
                logger.info(f" -Using user-provided gtf:{gtf_path}")
                gtf_path = gtf_to_all_genes(gtf_path)
            
            gtf_db_path = gtf_path[:-2]+"db"
            
            
            data = Genome(
                reference_name='GRCh37',
                annotation_name='Refseq',
                gtf_path_or_url=gtf_path
            )
            

            if os.path.isfile(gtf_db_path) is False:
                data.index()

            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:get_closest_gene(x,data=data,chrom=chrom,pos=pos,source=source,build=build), axis=1)), 
                index=output.index).values
            
        elif build=="38":
            
            logger.info(" -Assigning Gene name using NCBI refseq latest GRCh38 for protein coding genes")
            
            if gtf_path is None or not is_gtf_path:

                refseq38 = RefSeqFetcher(build = "38")

                refseq38.get_latest_release()
                refseq38.download_latest()
                refseq38.unzip_latest()
                refseq38.get_all_genes()

                gtf_path = refseq38.all_genes_path

            else:
                logger.info(f" -Using user-provided gtf:{gtf_path}")
                gtf_path = gtf_to_all_genes(gtf_path)

            gtf_db_path = gtf_path[:-2]+"db"
            
            data = Genome(
                reference_name='GRCh38',
                annotation_name='Refseq',
                gtf_path_or_url=gtf_path
            )
            
            if os.path.isfile(gtf_db_path) is False:
                data.index()
            
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:get_closest_gene(x,data=data,chrom=chrom,pos=pos,source=source,build=build), axis=1)), 
                index=output.index).values
            
    logger.info("Finished annotating variants with nearest gene name(s) successfully!")
    return output