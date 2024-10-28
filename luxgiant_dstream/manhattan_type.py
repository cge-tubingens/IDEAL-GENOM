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
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import textalloc as ta

from matplotlib.axes import Axes

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

    axis_center = pd.DataFrame(columns=[chr_col, 'center'])

    for i, chrom in enumerate(chromosomes):

        temp = data[data[chr_col] == chrom].reset_index(drop=True)

        axis_center.loc[i, chr_col] = chrom
        axis_center.loc[i, 'center'] = np.round((temp[chr_pos_col].max()+temp[chr_pos_col].min())/2,0)

    return axis_center

def manhattan_process_data(data_df:pd.DataFrame, chr_col:str='CHR', pos_col:str='POS', p_col:str='p')->dict:

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

def manhattan_draw(data_df:pd.DataFrame, snp_col:str, chr_col:str, pos_col:str, p_col:str, plot_dir:str, to_highlight:list=None, to_annotate:list=None, build:str='38', gtf_path:str=None)->bool:

    chr_colors           = ['#66c2a5', '#fc8d62']
    ylab                 = "-log10(p)"
    xlab                 = "Chromosome"
    genome_line          = 5e-8
    genome_line_color    = "#636363"
    suggestive_line      = 1e-5 
    suggestive_line_color= "#377eb8"

    # format data to draw manhattan plot
    plot_data = manhattan_process_data(
        data_df=data_df,
        chr_col=chr_col,
        pos_col=pos_col,
        p_col  =p_col
    )

    max_x_axis = plot_data['data']['rel_pos'].max()

    # Create the figure
    fig= plt.figure(figsize=(15, 9.2))
    ax = fig.add_subplot(111)

    # Suppress warnings about the number of chromosomes and just two colors
    warnings.filterwarnings("ignore", category=UserWarning)

    ax = sns.scatterplot(
        x      =plot_data['data']['rel_pos'], 
        y      =plot_data['data']['log10p'],
        hue    =plot_data['data'][chr_col], 
        palette=chr_colors, 
        ax     =ax, 
        s      =3, 
        legend =False
    )

    # set axis labels and font size
    ax.set_ylabel(ylab, fontsize=7)
    ax.set_xlabel(xlab, fontsize=7)

    # set axis limits
    ax.set_xlim(0, max_x_axis+1000)
    ax.set_ylim(0, plot_data['maxp']+10)

    # set x-axis ticks and labels
    x_ticks=plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis']['CHR'].astype(str).tolist()

    ax.set_xticks(ticks=x_ticks)
    ax.set_xticklabels(x_labels)

    # set ticks font size
    ax.tick_params(axis='both', labelsize=7)

    # add suggestive line
    if suggestive_line is not None:
        ax.axhline(
            -np.log10(suggestive_line), 
            color    =suggestive_line_color, 
            linestyle='dashed', 
            lw       =0.5
        )
    
    # add genome-wide line
    if genome_line is not None:
        ax.axhline(
            -np.log10(genome_line), 
            color    =genome_line_color, 
            linestyle='solid', 
            lw       =0.5
        )

    # highlight SNPs
    if to_highlight is not None:            

        plot_data['data']["HUE"] = pd.NA
        plot_data['data']["HUE"] = plot_data['data']["HUE"].astype("Int64")
        plot_data['data'].loc[plot_data['data'][snp_col].isin(to_highlight), "HUE"] = 0

        ax = sns.scatterplot(
                data     =plot_data['data'][plot_data['data']["HUE"]==0].reset_index(drop=True), 
                x        ='rel_pos',
                y        ='log10p',
                s        =10,
                ax       =ax,
                edgecolor="black",
                color   ='#CB132D',
            )

    # annotate SNPs   
    if to_annotate is not None:

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
                source ="refseq",
                verbose=True,
                gtf_path=gtf_path
            ).rename(columns={"GENE":"GENENAME"})

        texts = []  # a list to store text annotations for adjustment
        x = []      # a list to store x-coordinates for adjustment
        y = []      # a list to store y-coordinates for adjustment

        for i, row in variants_toanno.iterrows():

            x.append(row['rel_pos'])
            y.append(row['log10p'])
            texts.append(row['GENENAME'])

        x_lines_coor = np.linspace(0, max_x_axis, 1000).tolist() # list with a gris of x-coordinates for the lines

        ta.allocate(
            ax,              # the axis to which the text will be
            x        =x,     # x-coordinates of the data point to annotate
            y        =y,     # y-coordinates of the data point to annotate
            text_list=texts, # list of text to annotate
            x_scatter=plot_data['data']['rel_pos'], # all scatter points x-coordinates
            y_scatter=plot_data['data']['log10p'],  # all scatter points y-coordinates
            linecolor='black',                      # color of the line connecting the text to the data point
            textsize =7,                            # size of the text (Default to Nature standard)
            bbox     =dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.5),
            x_lines  = [x_lines_coor, x_lines_coor],
            y_lines  = [[suggestive_line]*len(x_lines_coor), [genome_line]*len(x_lines_coor)],
            avoid_label_lines_overlap =True,
            avoid_crossing_label_lines=True,
            min_distance=0.01,
            max_distance=0.4,
            margin      =0.01,
            rotation    =90
        )

    plt.tight_layout()
    plt.savefig(
        os.path.join(plot_dir, "manhattan_plot.jpeg"), dpi=600
    )
    plt.show()

    return True

def miami_process_data(data_top:pd.DataFrame, data_bottom:pd.DataFrame, chr_col:str, pos_col:str, p_col:str)->dict:
    
    """
    Processes Miami plot data by preparing, computing relative positions, and splitting the data.

    Parameters:
    -----------
    data_top (pd.DataFrame): 
        The top part of the data to be processed.
    data_bottom (pd.DataFrame): 
        The bottom part of the data to be processed.
    
    Returns:
    --------
    dict: A dictionary containing the processed data with the following keys:
        - 'upper': DataFrame containing the top part of the processed data.
        - 'lower': DataFrame containing the bottom part of the processed data.
        - 'axis': The center positions of the chromosomes.
        - 'maxp': The maximum -log10(p-value) in the data.
    """

    data_top['split_by']   = 'top'
    data_bottom['split_by']= 'bottom'

    data = pd.concat([data_top, data_bottom], axis=0, ignore_index=True)

    data = compute_relative_pos(data, chr_col=chr_col, pos_col=pos_col, p_col=p_col)

    axis_center = find_chromosomes_center(data, chr_col=chr_col, chr_pos_col='rel_pos')

    maxp = np.ceil(data['log10p'].max(skipna=True))

    miami_data = {
        'upper': data[data['split_by'] == 'top'].reset_index(drop=True),
        'lower': data[data['split_by'] == 'bottom'].reset_index(drop=True),
        'axis': axis_center,
        'maxp': maxp
    }

    return miami_data

def miami_annotate(axes:Axes, gwas_data:pd.DataFrame, annotations:pd.DataFrame, to_annotate:list)->Axes:

    snps = annotations['SNP'].to_list()
    genes= annotations['GENE'].to_list()
    highlighted_snps = pd.merge(annotations, gwas_data[['SNP', 'rel_pos', 'log10p']], on='SNP', how='inner')

    custom_hue_colors = {
        "on_both"      : "#1f77b4",  
        "top_in_bottom": "#2ca02c",
        "bottom_in_top": "#9467bd",
    }

    axes = sns.scatterplot(
        x       =highlighted_snps['rel_pos'], 
        y       =highlighted_snps['log10p'], 
        ax      =axes,
        hue     =highlighted_snps['type'],
        palette =custom_hue_colors,
        size    =10,
        legend  =False,
    )

    snps_to_annotate = highlighted_snps[highlighted_snps['type'].isin(to_annotate)].reset_index(drop=True)

    texts = []  # A list to store text annotations for adjustment
    x = []
    y = []
    for i, row in snps_to_annotate.iterrows():
        gene = genes[snps.index(row['SNP'])]  # Get corresponding gene name

        x.append(row['rel_pos'])
        y.append(row['log10p'])
        texts.append(gene)

    ta.allocate(
        axes,
        x        =x,
        y        =y,
        text_list=texts,
        x_scatter=x,
        y_scatter=y,
        linecolor='black',
        textsize =10,
        bbox     =dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.8),
    )

    return axes

def miami_classify_annotations(df_top_highlts:pd.DataFrame, df_bottom_highlts:pd.DataFrame)->dict:
    
    """
    Splits the annotation data into two parts for the top and bottom plots.

    Parameters:
    -----------
    df_top_highlts (pd.DataFrame): 
        The annotation data for the top plot.
    df_bottom_highlits (pd.DataFrame): 
        The annotation data for the bottom plot.
    
    Returns:
    --------
    tuple: A tuple containing the split annotation data for the top and bottom plots.
    """

    df_both = pd.merge(df_top_highlts, df_bottom_highlts[['SNP']], on='SNP', how='inner')
    df_both['type'] = 'on_both'

    df_top_in_bottom = df_top_highlts[~df_top_highlts['SNP'].isin(df_bottom_highlts['SNP'])].reset_index(drop=True)
    df_top_in_bottom['type'] = 'top_in_bottom'

    df_bottom_in_top = df_bottom_highlts[~df_bottom_highlts['SNP'].isin(df_top_highlts['SNP'])].reset_index(drop=True)
    df_bottom_in_top['type'] = 'bottom_in_top'


    return pd.concat([df_both, df_top_in_bottom, df_bottom_in_top], axis=0).reset_index(drop=True)

def miami_draw(df_top:pd.DataFrame, df_bottom:pd.DataFrame, snp_col:str, chr_col:str, pos_col:str, p_col:str, plots_dir:str, top_highlights:list=[], top_annotations:list=[], bottom_highlights:list=[], bottom_annotations:list=[], gtf_path:str=None)->bool:
    
    """
    Generates a Miami plot from two dataframes and saves the plot to the specified directory.

    Parameters:
    ----------
    df_top (pd.DataFrame): 
       DataFrame containing the data for the upper plot.
    df_bottom (pd.DataFrame): 
       DataFrame containing the data for the lower plot.
    plots_dir (str): 
       Directory where the plot image will be saved.
     
    Returns:
    -------
    bool: 
       True if the plot is successfully created and saved.
     
    The function creates a Miami plot, which is a type of scatter plot used in genomic studies to display 
    p-values from two different datasets. The plot consists of two panels: the upper panel for the first 
    dataset and the lower panel for the second dataset. The x-axis represents the genomic position, and 
    the y-axis represents the -log10(p) values. The function also adds genome-wide and suggestive significance 
    lines to the plot.
    """

    chr_colors           = ['#66c2a5', '#fc8d62']
    upper_ylab           = "-log10(p)" 
    lower_ylab           = "-log10(p)" 
    genome_line          = 5e-8
    genome_line_color    = "red"
    suggestive_line      = 1e-5 
    suggestive_line_color= "blue"

    # Suppress warnings about the number of chromosomes and just two colors
    warnings.filterwarnings("ignore", category=UserWarning)

    # format data to draw miami plot
    plot_data = miami_process_data(df_top, df_bottom, chr_col=chr_col, pos_col=pos_col, p_col=p_col)

    # Set axis labels for upper and lower plot
    def format_ylabel(label):
        return f"{label}\n-log10(p)" if label != "-log10(p)" else r"-log10(p)"
    
    upper_ylab = format_ylabel(upper_ylab)
    lower_ylab = format_ylabel(lower_ylab)

    max_x_axis = max(plot_data['upper']['rel_pos'].max(), plot_data['lower']['rel_pos'].max())+10

    # Create the figure
    plt.figure(figsize=(15, 18.4))

    # Create the upper plot

    ax_upper = plt.subplot(211)
    sns.scatterplot(x=plot_data['upper']['rel_pos'], y=plot_data['upper']['log10p'],
                    hue=plot_data['upper'][chr_col], palette=chr_colors, ax=ax_upper, s=1, legend=False)
    ax_upper.set_ylabel(upper_ylab)
    ax_upper.set_xlim(0, max_x_axis)

    x_ticks=plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis'][chr_col].astype(str).tolist()

    ax_upper.set_xticks(ticks=x_ticks)  # Set x-ticks
    ax_upper.set_xticklabels(x_labels)
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_upper.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.5)
    
    if genome_line is not None:
        ax_upper.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.5)

    # Create the lower plot
    ax_lower = plt.subplot(212)
    sns.scatterplot(x=plot_data['lower']['rel_pos'], y=plot_data['lower']['log10p'],
                    hue=plot_data['lower'][chr_col], palette=chr_colors, ax=ax_lower, s=1, legend=False)
    ax_lower.set_ylabel(lower_ylab)
    ax_lower.set_ylim(plot_data['maxp'], 0)  # Reverse y-axis
    ax_lower.set_xlim(0, max_x_axis)

    ax_lower.set_xticks(ticks=x_ticks) # Set x-ticks
    ax_lower.set_xticklabels([])
    ax_lower.xaxis.set_ticks_position('top')
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_lower.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.5)
    
    if genome_line is not None:
        ax_lower.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.5)
    
    top = set(top_highlights)
    bottom = set(bottom_highlights)

    if len(top_highlights)>0 or len(bottom_highlights)>0:

        both = list(top.intersection(bottom))
        top_only = list(top.difference(bottom))
        bottom_only = list(bottom.difference(top))

        plot_data['upper']['type'] = None
        plot_data['lower']['type'] = None

        plot_data['upper'].loc[plot_data['upper'][snp_col].isin(both), 'type'] = 'on_both'
        plot_data['lower'].loc[plot_data['lower'][snp_col].isin(both), 'type'] = 'on_both'

        plot_data['upper'].loc[plot_data['upper'][snp_col].isin(top_only), 'type'] = 'top_only'
        plot_data['lower'].loc[plot_data['lower'][snp_col].isin(top_only), 'type'] = 'top_only'

        plot_data['upper'].loc[plot_data['upper'][snp_col].isin(bottom_only), 'type'] = 'bottom_only'
        plot_data['lower'].loc[plot_data['lower'][snp_col].isin(bottom_only), 'type'] = 'bottom_only'

        custom_hue_colors = {
        "on_both"    : "#1f77b4",  
        "top_only"   : "#2ca02c",
        "bottom_only": "#9467bd",
        }

        ax_upper = sns.scatterplot(
            data=plot_data['upper'][~plot_data['upper']['type'].isnull()].reset_index(drop=True),
            x       ='rel_pos', 
            y       ='log10p', 
            ax      =ax_upper,
            hue     ='type',
            palette =custom_hue_colors,
            size    =10,
            legend  =False,
        )

        ax_lower = sns.scatterplot(
            data=plot_data['lower'][~plot_data['lower']['type'].isnull()].reset_index(drop=True),
            x       ='rel_pos', 
            y       ='log10p', 
            ax      =ax_lower,
            hue     ='type',
            palette =custom_hue_colors,
            size    =10,
            legend  =False,
        )

    if len(top_annotations)>0 or len(bottom_annotations)>0:

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

        top_variants_toanno = plot_data['upper'][plot_data['upper'][snp_col].isin(top_annotations)]\
            .reset_index(drop=True)
        
        bottom_variants_toanno = plot_data['lower'][plot_data['lower'][snp_col].isin(bottom_annotations)]\
            .reset_index(drop=True)
        
        if (top_variants_toanno.empty is not True):
            top_variants_toanno = annogene(
                top_variants_toanno,
                id     =snp_col,
                chrom  =chr_col,
                pos    =pos_col,
                log    =Log(),
                build  ='38',
                source ="refseq",
                verbose=True,
                gtf_path=gtf_path
            ).rename(columns={"GENE":"GENENAME"})

        if (bottom_variants_toanno.empty is not True):
            bottom_variants_toanno = annogene(
                bottom_variants_toanno,
                id     =snp_col,
                chrom  =chr_col,
                pos    =pos_col,
                log    =Log(),
                build  ='38',
                source ="refseq",
                verbose=True,
                gtf_path=gtf_path
            ).rename(columns={"GENE":"GENENAME"})

        x_lines_coor = np.linspace(0, max_x_axis, 1000).tolist() # list with a gris of x-coordinates for the lines

        texts_upper = []  # a list to store text annotations for adjustment
        x_upper = []      # a list to store x-coordinates for adjustment
        y_upper = []      # a list to store y-coordinates for adjustment

        for i, row in top_variants_toanno.iterrows():

            x_upper.append(row['rel_pos'])
            y_upper.append(row['log10p'])
            texts_upper.append(row['GENENAME'])

        texts_lower = []  # a list to store text annotations for adjustment
        x_lower = []      # a list to store x-coordinates for adjustment
        y_lower = []      # a list to store y-coordinates for adjustment

        for i, row in bottom_variants_toanno.iterrows():

            x_lower.append(row['rel_pos'])
            y_lower.append(row['log10p'])
            texts_lower.append(row['GENENAME'])

        ta.allocate(
            ax_upper,              # the axis to which the text will be
            x        =x_upper,     # x-coordinates of the data point to annotate
            y        =y_upper,     # y-coordinates of the data point to annotate
            text_list=texts_upper, # list of text to annotate
            #x_scatter=plot_data['upper']['rel_pos'], # all scatter points x-coordinates
            #y_scatter=plot_data['upper']['log10p'],  # all scatter points y-coordinates
            linecolor='black',                      # color of the line connecting the text to the data point
            textsize =7,                            # size of the text (Default to Nature standard)
            bbox     =dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.5),
            x_lines  = [x_lines_coor, x_lines_coor],
            y_lines  = [[-np.log10(suggestive_line)]*len(x_lines_coor), [-np.log10(genome_line)]*len(x_lines_coor)],
            avoid_label_lines_overlap =True,
            avoid_crossing_label_lines=True,
            min_distance=0.01,
            max_distance=0.4,
            margin      =0.01,
            rotation    =90
        )

        ta.allocate(
            ax_lower,              # the axis to which the text will be
            x        =x_lower,     # x-coordinates of the data point to annotate
            y        =y_lower,     # y-coordinates of the data point to annotate
            text_list=texts_lower, # list of text to annotate
            #x_scatter=plot_data['lower']['rel_pos'], # all scatter points x-coordinates
            #y_scatter=plot_data['loweer']['log10p'],  # all scatter points y-coordinates
            linecolor='black',                      # color of the line connecting the text to the data point
            textsize =7,                            # size of the text (Default to Nature standard)
            bbox     =dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.5),
            x_lines  = [x_lines_coor, x_lines_coor],
            y_lines  = [[suggestive_line]*len(x_lines_coor), [genome_line]*len(x_lines_coor)],
            avoid_label_lines_overlap =True,
            avoid_crossing_label_lines=True,
            min_distance=0.01,
            max_distance=0.4,
            margin      =0.01,
            rotation    =90
        )


    ax_lower.set_xlabel("Base pair position")
    ax_upper.set_xlabel("")
    
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'miami_plot.png'))
    plt.show()

    return True
