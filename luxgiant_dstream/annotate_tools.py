import requests
import pandas as pd
import numpy as np

def query_ensembl_gene_overlap(chr, start, end, build=38):
    if build == 38:
        server = "http://rest.ensembl.org"
    elif build == 37:
        server = "http://grch37.rest.ensembl.org"
    
    ext = f"/overlap/region/human/{chr}:{start}-{end}?feature=gene"
    headers = {"Content-Type": "application/json"}
    r = requests.get(f"{server}{ext}", headers=headers)
    
    # Check if request was successful
    r.raise_for_status()
    
    restr = r.json()
    
    # Check if the response is empty
    if not isinstance(restr, list) or len(restr) == 0:
        return pd.DataFrame()
    
    # Convert to DataFrame
    return pd.DataFrame(restr)

def process_overlap_restr(restr):
    # Filter for protein-coding genes
    if 'biotype' not in restr.columns:
        return ''
    prot_genes = restr[restr['biotype'] == 'protein_coding']
    
    # Check if 'external_name' column exists and is not null
    if 'external_name' in prot_genes.columns and not prot_genes['external_name'].isnull().all():
        gene_names = prot_genes['external_name']
    else:
        gene_names = prot_genes['id']
    
    # Return the gene names as a comma-separated string
    return ','.join(gene_names)

def process_case2_restr(restr, pos):
    # Filter for protein-coding genes
    restr = restr[restr['biotype'] == 'protein_coding'].copy()
    print("Overlap no, distance")

    # Compute distances from the variant position
    restr['dist1'] = np.abs(restr['start'] - pos)
    restr['dist2'] = np.abs(restr['end'] - pos)
    restr['dist'] = restr[['dist1', 'dist2']].min(axis=1)

    # Order by distance
    reordered_restr = restr.sort_values(by='dist')

    if 'external_name' not in reordered_restr.columns or reordered_restr['external_name'].isnull().all():
        # If all protein-coding genes have null external_name, return the closest gene's Ensembl Stable ID
        gene = reordered_restr.iloc[0]['gene_id']
        dist = reordered_restr.iloc[0]['dist']
        return {'gene': gene, 'dist': dist}

    for _, row in reordered_restr.iterrows():
        gene = row['external_name']
        dist = row['dist']
        if pd.isnull(gene):
            print('[INFO] Skipping gene with no external_name')
            print(row)
            continue
        break

    print(f"Minimum reached for gene {gene} at {dist}")
    
    assert gene is not None and dist is not None, "Gene or distance should not be None"
    print(dist)
    
    return {'gene': gene, 'dist': dist}

def get_variant_context(chr, pos, build=38):
    print(f"Getting context for {chr}:{pos}")

    # Query for the exact position
    restr = query_ensembl_gene_overlap(chr, pos, pos, build)

    # Three possible cases:
    # 1. Variant directly overlaps region of protein-coding gene(s)
    # 2. Variant is intergenic, but protein-coding gene is within a 2Mb window (1Mb either side)
    # 3. Variant is intergenic, and no protein-coding gene present within 2Mb window (1Mb either side)
    
    # Case 1: Check for overlap with a protein-coding gene
    if 'biotype' in restr.columns and 'protein_coding' in restr['biotype'].values:
        # We have at least one protein-coding gene overlapping
        print("Overlap found")
        gene = process_overlap_restr(restr)
        dist = 0
        assert gene is not None  # Ensures gene is not null

    else:
        # Case 2 and 3: Look for genes within 1Mb either side
        st = max(1, pos - 10**6)  # Ensure start is not less than 1
        restr = query_ensembl_gene_overlap(chr, st, pos + 10**6, build)

        if restr.empty or 'protein_coding' not in restr['biotype'].values:
            # Case 3: No protein-coding gene within the window
            return ["none", 0, "intergenic_variant"]
        else:
            # Case 2: Protein-coding gene found within the window
            output = process_case2_restr(restr, pos)
            gene = output['gene']
            dist = output['dist']
    
    return [gene, dist, "coding_variant" if dist == 0 else "near_gene_variant"]