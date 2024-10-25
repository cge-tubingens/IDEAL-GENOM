"""
This module provides functions to process and visualize Manhattan plots for genomic data. It includes functions to compute relative positions of SNPs, find chromosome centers, process data for Manhattan plots, and draw the plots with optional highlighting and annotation.

Functions:
- compute_relative_pos: Compute relative positions and -log10(p-values) for SNPs.
- find_chromosomes_center: Calculate center positions of chromosomes.
- process_manhattan_data: Prepare data for Manhattan plot visualization.
- draw_manhattan: Generate and save a Manhattan plot with optional SNP highlighting and annotation.
"""

import gzip
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import textalloc as ta

from gwaslab.bd_download import download_file
from gwaslab.g_Log import Log
from gwaslab.util_in_get_sig import annogene

def compute_relative_pos(data:pd.DataFrame, chr_col:str='CHR', pos_col:str='POS', p_col:str='p')->pd.DataFrame:
    
    """
    Compute the relative position of probes/SNPs across chromosomes and add a -log10(p-value) column.

    Parameters:
    -----------
    data (pd.DataFrame): 
        Input DataFrame containing genomic data.
    chr_col (str): 
        Column name for chromosome identifiers. Default is 'CHR'.
    pos_col (str): 
        Column name for base pair positions. Default is 'POS'.
    p_col (str): 
        Column name for p-values. Default is 'p'.
    
    Returns:
    --------
    pd.DataFrame: DataFrame with additional columns for relative positions and -log10(p-values).
    """

    # Group by chromosome and compute chromosome size
    chr_grouped = data.groupby(chr_col).agg(chrlength=(pos_col, 'max')).reset_index()

    # Calculate cumulative chromosome length
    chr_grouped['cumulativechrlength'] = chr_grouped['chrlength'].cumsum() - chr_grouped['chrlength']

    # Merge cumulative chromosome length back to the original data
    data = pd.merge(data, chr_grouped[[chr_col, 'cumulativechrlength']], on=chr_col)

    # Sort by chromosome and position
    data = data.sort_values(by=[chr_col, pos_col])

    # Add the relative position of the probe/snp
    data['rel_pos'] = data[pos_col] + data['cumulativechrlength']

    # Drop cumulative chromosome length column
    data = data.drop(columns=['cumulativechrlength'])

    data['log10p']= -np.log10(data[p_col])

    return data

def find_chromosomes_center(data:pd.DataFrame, chr_col:str='CHR', chr_pos_col:str='rel_pos')->pd.DataFrame:
    
    """
    Calculate the center positions of chromosomes in a given DataFrame. This function takes a DataFrame containing chromosome data and calculates the center position for each chromosome based on the specified chromosome column and chromosome position column.
    
    Parameters:
    -----------
    data : pd.DataFrame
        The input DataFrame containing chromosome data.
    chr_col : str, optional
        The name of the column representing chromosome identifiers (default is 'CHR').
    chr_pos_col : str, optional
        The name of the column representing relative positions within chromosomes (default is 'rel_pos').
    
    Returns:
    --------
    pd.DataFrame
        A DataFrame with columns 'CHR' and 'center', where 'CHR' contains chromosome identifiers and 'center' contains the calculated center positions for each chromosome.
    """

    chromosomes = data[chr_col].unique()

    axis_center = pd.DataFrame(columns=['CHR', 'center'])

    for i, chrom in enumerate(chromosomes):

        temp = data[data[chr_col] == chrom].reset_index(drop=True)

        axis_center.loc[i, 'CHR'] = chrom
        axis_center.loc[i, 'center'] = np.round((temp[chr_pos_col].max()+temp[chr_pos_col].min())/2,0)

    return axis_center

def process_manhattan_data(data_df:pd.DataFrame, chr_col:str='CHR', pos_col:str='POS', p_col:str='p')->dict:

    data = compute_relative_pos(
        data_df, 
        chr_col=chr_col, 
        pos_col=pos_col, 
        p_col  =p_col
    )

    axis_center = find_chromosomes_center(data, chr_col=chr_col)

    maxp = np.ceil(data['log10p'].max(skipna=True))

    manhattan_data = {
        'data': data.reset_index(drop=True),
        'axis': axis_center,
        'maxp': maxp
    }

    return manhattan_data

def draw_manhattan(data_df:pd.DataFrame, snp_col:str, chr_col:str, pos_col:str, p_col:str, plot_dir:str, to_highlight:list=None, to_annotate:list=None, build:str='38', gtf_path:str=None)->bool:

    chr_colors           = ['#66c2a5', '#fc8d62']
    ylab                 = "-log10(p)"
    xlab                 = "Chromosome"
    genome_line          = 5e-8
    genome_line_color    = "red"
    suggestive_line      = 1e-5 
    suggestive_line_color= "blue"

    # format data to draw manhattan plot
    plot_data = process_manhattan_data(
        data_df=data_df,
        chr_col=chr_col,
        pos_col=pos_col,
        p_col=p_col
    )

    max_x_axis = plot_data['data']['rel_pos'].max()

    # Create the figure
    fig= plt.figure(figsize=(7.2, 4.7))
    ax = fig.add_subplot(111)

    ax = sns.scatterplot(
        x=plot_data['data']['rel_pos'], 
        y=plot_data['data']['log10p'],
        hue=plot_data['data'][chr_col], 
        palette=chr_colors, 
        ax=ax, 
        s=2, 
        legend=False
    )

    ax.set_ylabel(ylab, fontsize=7)
    ax.set_xlabel(xlab, fontsize=7)
    ax.set_xlim(0, max_x_axis+10)

    x_ticks=plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis']['CHR'].astype(str).tolist()

    ax.set_xticks(ticks=x_ticks)  # Set x-ticks
    ax.set_xticklabels(x_labels)
    ax.tick_params(axis='both', labelsize=7)

    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax.axhline(
            -np.log10(suggestive_line), 
            color    =suggestive_line_color, 
            linestyle='solid', 
            lw       =0.5
        )
    
    if genome_line is not None:
        ax.axhline(
            -np.log10(genome_line), 
            color    =genome_line_color, 
            linestyle='dashed', 
            lw       =0.5
        )

    if to_highlight is not None:

        plot_data['data']["HUE"] = pd.NA
        plot_data['data']["HUE"] = plot_data['data']["HUE"].astype("Int64")
        plot_data['data'].loc[plot_data['data'][snp_col].isin(to_highlight), "HUE"] = 0

        ax = sns.scatterplot(
                data     =plot_data['data'][plot_data['data']["HUE"]==0].reset_index(drop=True), 
                x        ='rel_pos',
                y        ='log10p',
                s        =5,
                ax       =ax,
                edgecolor="black",
                color   ='#CB132D',
            )
        
        plot_data['data'] = plot_data['data'].rename(columns={'log10p':"scaled_P"}, inplace=False)
        
        plot_data['data'], chrom_df = _quick_assign_i_with_rank(
            plot_data['data'], 
            chrpad=0.03, 
            use_rank=False, 
            chrom=chr_col,
            pos=pos_col,
            drop_chr_start=False,
            _posdiccul=None
        )

        highlight_i = plot_data['data'].loc[plot_data['data'][snp_col].isin(to_annotate), "i"].tolist()
        
    if to_annotate is not None:

        variants_toanno = plot_data['data'][plot_data['data'][snp_col].isin(to_annotate)]\
            .reset_index(drop=True)
        
        if (variants_toanno.empty is not True):
            variants_toanno = annogene(
                variants_toanno,
                id     =snp_col,
                chrom  =chr_col,
                pos    =pos_col,
                log    =Log(),
                build  =build,
                source ="ensembl",
                verbose=False
            ).rename(columns={"GENE":"GENENAME"})

        #return plot_data['data']

        ax = annotate_single(
            sumstats=plot_data['data'],
            anno='GENENAME',
            mode="m",
            ax1=ax,
            highlight_i=highlight_i,
            highlight_chrpos=False,
            highlight_anno_args=None,
            to_annotate=variants_toanno,
            anno_d=dict(),
            anno_alias=dict(),
            anno_style="right",
            anno_args=dict(),
            arm_scale=1,
            anno_max_iter=100,
            arm_scale_d=None,
            arm_offset=50,
            anno_adjust=False,
            anno_fixed_arm_length=None,
            maxy=plot_data['maxp'],
            anno_fontsize= 8,
            font_family="DejaVu sans",
            region=None,
            region_anno_bbox_args=None,
            skip=0,
            anno_height=1,
            snpid=snp_col,
            chrom=chr_col,
            pos=pos_col,
            repel_force=0.03,
            verbose=False,
            log=Log(),
            _invert=False
        )    

    plt.tight_layout()
    plt.show()

    return plot_data['data']
