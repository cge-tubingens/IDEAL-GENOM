"""
This module provides functions to process and visualize Manhattan plots for genomic data. It includes functions to compute relative positions of SNPs, find chromosome centers, process data for Manhattan plots, and draw the plots with optional highlighting and annotation.

Functions:
- compute_relative_pos: Compute relative positions and -log10(p-values) for SNPs.
- find_chromosomes_center: Calculate center positions of chromosomes.
- process_manhattan_data: Prepare data for Manhattan plot visualization.
- draw_manhattan: Generate and save a Manhattan plot with optional SNP highlighting and annotation.
"""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import textalloc as ta

from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase

from ideal_genom.annotations import annotate_snp

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

    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")

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

    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if chr_pos_col not in data.columns:
        raise ValueError(f"Column '{chr_pos_col}' not found in the input DataFrame.")

    chromosomes = data[chr_col].unique()

    axis_center = pd.DataFrame(columns=[chr_col, 'center'])

    for i, chrom in enumerate(chromosomes):

        temp = data[data[chr_col] == chrom].reset_index(drop=True)

        axis_center.loc[i, chr_col] = chrom
        axis_center.loc[i, 'center'] = np.round((temp[chr_pos_col].max()+temp[chr_pos_col].min())/2,0)

    return axis_center

def manhattan_process_data(data_df:pd.DataFrame, chr_col:str='CHR', pos_col:str='POS', p_col:str='p')->dict:
    
    """
    Processes the input DataFrame to prepare data for a Manhattan plot.

    Parameters:
    -----------
    data_df : pd.DataFrame
        The input DataFrame containing genomic data.
    chr_col : str (optional)
        The column name for chromosome data. Defaults to 'CHR'.
    pos_col : str (optional)
        The column name for position data. Defaults to 'POS'.
    p_col : str (optional)
        The column name for p-value data. Defaults to 'p'.

    Returns:
    --------
        dict: A dictionary containing processed data for the Manhattan plot with the following keys:
            - 'data' (pd.DataFrame): The processed DataFrame with relative positions and log-transformed p-values.
            - 'axis' (dict): The center positions of each chromosome for plotting.
            - 'maxp' (float): The maximum log-transformed p-value.
    """

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data_df.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")

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

def manhattan_draw(data_df:pd.DataFrame, snp_col:str, chr_col:str, pos_col:str, p_col:str, plot_dir:str, to_highlight:pd.DataFrame=pd.DataFrame(), highlight_hue:str='hue', to_annotate:pd.DataFrame=pd.DataFrame(), gen_col:str=None, build:str='38', anno_source='ensembl', gtf_path:str=None, save_name:str='manhattan_plot.jpeg', genome_line: float = 5e-8, suggestive_line: str = 1e-5, yaxis_margin: float = 10, dpi: int = 500)->bool:

    """
    Draws a Manhattan plot for visualizing GWAS results.

    Parameters:
    -----------
    data_df : pd.DataFrame
        The input DataFrame containing genomic data.
    snp_col : str    
        The column name for SNP identifiers.
    chr_col : str
        The column name for chromosome identifiers.
    pos_col : str
        The column name for base pair positions.
    p_col : str
        The column name for p-values.
    plot_dir : str
        The directory where the plot will be saved.
    to_highlight : list, optional
        A list of SNP identifiers to highlight in the plot. Default is an empty list.
    to_annotate : list, optional
        A list of SNP identifiers to annotate in the plot. Default is an empty list.
    build : str, optional
        The genome build version. Default is '38'.
    gtf_path : str, optional
        The path to the GTF file for gene annotation. If None, the file will be downloaded. Default is None.
    save_name : str, optional
        The name of the file to save the plot as. Default is 'manhattan_plot.jpeg'.

    Returns:
    --------
    bool
        True if the plot is successfully created and saved, False otherwise.
    """

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data_df.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")
    
    if to_highlight is not None:
        if not isinstance(to_highlight, pd.DataFrame):
            raise TypeError("to_highlight must be a list of SNP identifiers.")
    if to_annotate is not None:
        if not isinstance(to_annotate, pd.DataFrame):
            raise TypeError("to_annotate must be a data frame of SNP identifiers.")
    
    if snp_col not in data_df.columns:
        raise ValueError(f"Column '{snp_col}' not found in the input DataFrame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data_df.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")
    
    if not os.path.exists(plot_dir):
        raise FileNotFoundError(f"Directory '{plot_dir}' not found.")
    
    chr_colors           = ['#66c2a5', '#fc8d62']
    ylab                 = "-log10(p)"
    xlab                 = "Chromosome"
    genome_line_color    = "#636363"
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
    fig= plt.figure(figsize=(15, 10))
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
        legend =False,
        edgecolor='none'
    )

    # set axis labels and font size
    ax.set_ylabel(ylab, fontsize=7)
    ax.set_xlabel(xlab, fontsize=7)

    # set axis limits
    ax.set_xlim(0, max_x_axis+1000)
    ax.set_ylim(0, plot_data['maxp'] + yaxis_margin)

    # set x-axis ticks and labels
    x_ticks=plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis'][chr_col].astype(str).tolist()

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

        if to_highlight.empty is not True: 

            df_subset = plot_data['data'].merge(to_highlight, on=snp_col, how='inner')

            if df_subset.empty is not True:
            
                if highlight_hue in df_subset.columns:
                    ax = sns.scatterplot(
                        data     =df_subset, 
                        x        ='rel_pos',
                        y        ='log10p',
                        s        =20,
                        ax       =ax,
                        edgecolor="black",
                        hue      =highlight_hue,
                        palette='Set1',
                    )
                    ax.legend(loc='best', fontsize=9, title=None)
                else:
                    ax = sns.scatterplot(
                            data     =df_subset, 
                            x        ='rel_pos',
                            y        ='log10p',
                            s        =20,
                            ax       =ax,
                            edgecolor="black"
                        )
                
    plt.tight_layout()

    # annotate SNPs   
    if to_annotate is not None and to_annotate.empty is not True:

        if gen_col is not None:

            variants_toanno = plot_data['data'][plot_data['data'][snp_col].isin(to_annotate[snp_col])]\
                .reset_index(drop=True)
            variants_toanno = variants_toanno.merge(to_annotate, on=snp_col, how='left')
            variants_toanno = variants_toanno.rename(columns={gen_col:"GENENAME"})

        else:

            variants_toanno = plot_data['data'][plot_data['data'][snp_col].isin(to_annotate[snp_col])]\
                .reset_index(drop=True)
            
            if (variants_toanno.empty is not True):
            
                variants_toanno = annotate_snp(
                    insumstats=variants_toanno,
                    chrom     =chr_col,
                    pos       =pos_col,
                    build     =build,
                    source    =anno_source,
                    gtf_path  =gtf_path
                ).rename(columns={"GENE":"GENENAME"})

        ax, texts = manhattan_type_annotate(
            axes           =ax, 
            data           =plot_data['data'], 
            variants_toanno=variants_toanno, 
            max_x_axis     =max_x_axis, 
            suggestive_line=suggestive_line, 
            genome_line    =genome_line
        )


    # save the plot

    plt.savefig(
        os.path.join(plot_dir, save_name), dpi=dpi
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

    if not isinstance(data_top, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if not isinstance(data_bottom, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data_top.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_top.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data_top.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")

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

def manhattan_type_annotate(axes:Axes, data:pd.DataFrame, variants_toanno:pd.DataFrame, max_x_axis:float, suggestive_line:float, genome_line:float)->Axes:
    
    """
    Annotates a Manhattan plot with gene names.

    Parameters:
    -----------
    axes : Axes (matplotlib.axes.Axes)
        The matplotlib axes object where the Manhattan plot is drawn.
    data : pd.DataFrame
        DataFrame containing the scatter plot data with columns 'rel_pos' and 'log10p'.
    variants_toanno : pd.DataFrame
        DataFrame containing the variants to annotate with columns 'rel_pos', 'log10p', and 'GENENAME'.
    max_x_axis : float
        The maximum value for the x-axis.
    suggestive_line : float
        The y-value for the suggestive significance line.
    genome_line : float
        The y-value for the genome-wide significance line.

    Returns:
    --------
    Axes
        The matplotlib axes object with annotations.
    list
        A list of text objects for the annotations.
    """

    if not isinstance(axes, Axes):
        raise TypeError("axes must be a matplotlib Axes object.")
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")
    if not isinstance(variants_toanno, pd.DataFrame):
        raise TypeError("variants_toanno must be a pandas DataFrame.")
    if not isinstance(max_x_axis, (int, float, np.int64, np.float64)):
        raise TypeError("max_x_axis must be an integer or float.")
    if not isinstance(suggestive_line, float):
        raise TypeError("suggestive_line must be a float.")
    if not isinstance(genome_line, float):
        raise TypeError("genome_line must be a float.")
    
    x_lines_coor = np.linspace(0, max_x_axis, 1000).tolist() # list with a grid of x-coordinates for the lines
    
    texts = []  # a list to store text annotations for adjustment
    x = []      # a list to store x-coordinates for adjustment
    y = []      # a list to store y-coordinates for adjustment

    for i, row in variants_toanno.iterrows():

        x.append(row['rel_pos'])
        y.append(row['log10p'])
        texts.append(row['GENENAME'])

    allocate = ta.allocate(
            axes,            
            x        =x,     # x-coordinates of the data point to annotate
            y        =y,     # y-coordinates of the data point to annotate
            text_list=texts, # list of text to annotate
            x_scatter=data['rel_pos'], # all scatter points x-coordinates
            y_scatter=data['log10p'],  # all scatter points y-coordinates
            linecolor='black',                      # color of the line connecting the text to the data point
            textsize =8,                            # size of the text (Default to Nature standard)
            bbox     =dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.5),
            x_lines  = [x_lines_coor, x_lines_coor],
            y_lines  = [[-np.log10(suggestive_line)]*len(x_lines_coor), [-np.log10(genome_line)]*len(x_lines_coor)],
            avoid_label_lines_overlap =True,
            avoid_crossing_label_lines=True,
            min_distance=0.01,
            max_distance=0.4,
            margin      =0.01,
            rotation    =90,
            nbr_candidates=300,
            priority_strategy=42,
            plot_kwargs=dict(linestyle=':')
        )

    text_objs = allocate[2]
        
    return axes, text_objs

def miami_draw_anno_lines(renderer:RendererBase, axes:Axes, texts:list, variants_toanno:pd.DataFrame):
    
    """
    Draws annotation lines from text labels to their corresponding data points on a plot.

    Parameters:
    ----------
    renderer : RendererBase 
        The renderer used to draw the plot.
    axes : Axes (matplotlib.axes.Axes) 
        The axes on which the plot is drawn.
    texts : list
        A list of text objects to annotate.
    variants_toanno : pd.DataFrame
        A DataFrame containing the data points to annotate, with columns 'GENENAME', 'rel_pos', and 'log10p'.

    Returns:
    -------
    Axes: The axes with the annotation lines drawn.
    """

    from math import dist

    for k in range(len(texts)):

        text_obj = texts[k]
        bbox = text_obj.get_window_extent(renderer=renderer)
        
        # Transform to data coordinates the corners of the text bbox
        data_coords = axes.transData.inverted().transform(bbox.corners()) # list with coordinates of the text box

        anno = text_obj.get_text()

        data_x = variants_toanno.loc[variants_toanno['GENENAME'] == anno, 'rel_pos'].values[0]
        data_y = variants_toanno.loc[variants_toanno['GENENAME'] == anno, 'log10p'].values[0]

        closest_point = data_coords[0]
        for point in data_coords:
            if dist([data_x, data_y], point) < dist([data_x, data_y], closest_point):
                closest_point = point 

        axes.plot(
            [data_x, closest_point[0]], 
            [data_y, closest_point[1]], 
            color='black', 
            linestyle=':', 
            linewidth=1
        )

    return axes

def miami_draw(df_top: pd.DataFrame, df_bottom: pd.DataFrame, snp_col: str, chr_col: str, pos_col: str, p_col: str, plots_dir: str, top_highlights: list = [], top_annotations: pd.DataFrame=pd.DataFrame(), bottom_highlights: list = [], bottom_annotations: pd.DataFrame = pd.DataFrame(), top_gen_col: str = None, bottom_gen_col: str = None, gtf_path: str = None, source: str = "ensemble", build: str = '38', save_name: str = 'miami_plot.jpeg', legend_top: str = 'top GWAS', legend_bottom: str = 'bottom GWAS', dpi: int = 500)->bool:
    
    """
    Draws a Miami plot (a combination of two Manhattan plots) for visualizing GWAS results.

    Parameters:
    -----------
    df_top : pd.DataFrame
        DataFrame containing the top plot data.
    df_bottom : pd.DataFrame
        DataFrame containing the bottom plot data.
    snp_col : str
        Column name for SNP identifiers.
    chr_col : str
        Column name for chromosome identifiers.
    pos_col : str
        Column name for base pair positions.
    p_col : str
        Column name for p-values.
    plots_dir : str
        Directory where the plot will be saved.
    top_highlights : list, optional
        List of SNPs to highlight in the top plot.
    top_annotations : list, optional
        List of SNPs to annotate in the top plot.
    bottom_highlights : list, optional
        List of SNPs to highlight in the bottom plot.
    bottom_annotations : list, optional
        List of SNPs to annotate in the bottom plot.
    gtf_path : str, optional
        Path to the GTF file for gene annotation. If None, the file will be downloaded.
    save_name : str, optional
        Name of the file to save the plot as. Default is 'miami_plot.jpeg'.
    Returns:
    --------
    bool
        True if the plot is successfully created and saved, False otherwise.
    """
    
    if not isinstance(df_top, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if not isinstance(df_bottom, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in df_top.columns or chr_col not in df_bottom.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in df_top.columns or pos_col not in df_bottom.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in df_top.columns or p_col not in df_bottom.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")
    if not os.path.exists(plots_dir):
        raise FileNotFoundError(f"Directory '{plots_dir}' not found.")
    
    annot_high = [top_highlights,  bottom_highlights]
    for lst in annot_high:
        if not isinstance(lst, list):
            raise TypeError("Annotation lists must be lists of SNP identifiers.")
        for val in lst:
            if not isinstance(val, str):
                raise TypeError("Annotation lists must be lists of SNP identifiers.")
    
    if not isinstance(top_annotations, pd.DataFrame):
        raise TypeError("Annotation data must be a pandas DataFrame.")
    if not isinstance(bottom_annotations, pd.DataFrame):
        raise TypeError("Annotation data must be a pandas DataFrame.")
    
    if not isinstance(save_name, str):
        raise TypeError("save_name must be a string.")
    if not isinstance(legend_top, str):
        raise TypeError("legend_top must be a string.")
    if not isinstance(legend_bottom, str):
        raise TypeError("legend_bottom must be a string.")

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
    fig = plt.figure(figsize=(20, 13))

    # Create the upper plot

    ax_upper = plt.subplot(211)

    sns.scatterplot(
        x      =plot_data['upper']['rel_pos'], 
        y      =plot_data['upper']['log10p'],
        hue    =plot_data['upper'][chr_col], 
        palette=chr_colors, 
        ax     =ax_upper, 
        s      =3, 
        legend =False,
        edgecolor='none'
    )
    ax_upper.set_ylabel(upper_ylab)
    ax_upper.set_xlim(0, max_x_axis)

    x_ticks =plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis'][chr_col].astype(str).tolist()

    ax_upper.set_xticks(ticks=x_ticks)  # Set x-ticks
    ax_upper.set_xticklabels(x_labels)
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_upper.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.7)
    
    if genome_line is not None:
        ax_upper.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.7)

    # Create the lower plot
    ax_lower = plt.subplot(212)
    sns.scatterplot(
        x      =plot_data['lower']['rel_pos'], 
        y      =plot_data['lower']['log10p'],
        hue    =plot_data['lower'][chr_col], 
        palette=chr_colors, 
        ax     =ax_lower, 
        s      =3, 
        legend =False,
        edgecolor='none'
    )
    ax_lower.set_ylabel(lower_ylab)
    ax_lower.set_ylim(plot_data['maxp'], 0)  # Reverse y-axis
    ax_lower.set_xlim(0, max_x_axis)

    ax_lower.set_xticks(ticks=x_ticks) # Set x-ticks
    ax_lower.set_xticklabels([])
    ax_lower.xaxis.set_ticks_position('top')
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_lower.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.7)
    
    if genome_line is not None:
        ax_lower.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.7)
    
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
        "on_both"    : "#7B68EE",  
        "top_only"   : "#008080",
        "bottom_only": "#FF1493",
        }

        ax_upper = sns.scatterplot(
            data    =plot_data['upper'][~plot_data['upper']['type'].isnull()].reset_index(drop=True),
            x       ='rel_pos', 
            y       ='log10p', 
            ax      =ax_upper,
            hue     ='type',
            palette =custom_hue_colors,
            size    =10,
            legend  =False,
            edgecolor='none'
        )

        ax_lower = sns.scatterplot(
            data    =plot_data['lower'][~plot_data['lower']['type'].isnull()].reset_index(drop=True),
            x       ='rel_pos', 
            y       ='log10p', 
            ax      =ax_lower,
            hue     ='type',
            palette =custom_hue_colors,
            size    =10,
            legend  =False,
            edgecolor='none'
        )

    # adjust layout
    plt.tight_layout()

    # render and draw the figure without showing it 
    r = fig.canvas.get_renderer()
    fig.canvas.draw()

    if top_annotations.shape[0]>0:
        
        if top_gen_col is not None:
            
            top_variants_toanno = plot_data['upper'][plot_data['upper'][snp_col].isin(top_annotations[snp_col])]\
                .reset_index(drop=True)
            top_variants_toanno = top_variants_toanno.merge(top_annotations, on=snp_col, how='left')
            top_variants_toanno = top_variants_toanno.rename(columns={top_gen_col:"GENENAME"})

        else:

            top_variants_toanno = plot_data['upper'][plot_data['upper'][snp_col].isin(top_annotations[snp_col])]\
                .reset_index(drop=True)
                
            # get gene names for upper plot
            if (top_variants_toanno.empty is not True):
                top_variants_toanno = annotate_snp(
                    top_variants_toanno,
                    chrom  =chr_col,
                    pos    =pos_col,
                    build  =build,
                    source =source,
                    gtf_path=gtf_path
                ).rename(columns={"GENE":"GENENAME"})

        # annotate upper plot
        ax_upper, texts_upper = manhattan_type_annotate(
            axes           =ax_upper, 
            data           =plot_data['upper'], 
            variants_toanno=top_variants_toanno, 
            max_x_axis     =max_x_axis, 
            suggestive_line=suggestive_line, 
            genome_line    =genome_line
        )

    if bottom_annotations.shape[0]>0:

        if bottom_gen_col is not None:

            bottom_variants_toanno = plot_data['lower'][plot_data['lower'][snp_col].isin(bottom_annotations[snp_col])]\
                .reset_index(drop=True)
            bottom_variants_toanno = bottom_variants_toanno.merge(bottom_annotations, on=snp_col, how='left')
            bottom_variants_toanno = bottom_variants_toanno.rename(columns={bottom_gen_col:"GENENAME"})
        
        else:

            bottom_variants_toanno = plot_data['lower'][plot_data['lower'][snp_col].isin(bottom_annotations[snp_col])]\
            .reset_index(drop=True)

            # get gane names for bottom plot
            if (bottom_variants_toanno.empty is not True):
                bottom_variants_toanno = annotate_snp(
                    bottom_variants_toanno,
                    chrom  =chr_col,
                    pos    =pos_col,
                    build  =build,
                    source =source,
                    gtf_path=gtf_path
                ).rename(columns={"GENE":"GENENAME"})

        
        # annotate lower plot
        ax_lower, texts_lower = manhattan_type_annotate(
            axes           =ax_lower, 
            data           =plot_data['lower'], 
            variants_toanno=bottom_variants_toanno, 
            max_x_axis     =max_x_axis, 
            suggestive_line=suggestive_line, 
            genome_line    =genome_line
        )

    from matplotlib.lines import Line2D

    if len(top_highlights)>0 or len(bottom_highlights)>0:

        if len(both)>0:

            custom_dots = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#7B68EE", markersize=5, label='Top Hits both'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#008080", markersize=5, label=f'Top Hits {legend_top}'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#FF1493", markersize=5, label=f'Top Hits {legend_bottom}'),
            ]

        else:
            custom_dots = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#008080", markersize=5, label=f'Top Hits {legend_top}'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#FF1493", markersize=5, label=f'Top Hits {legend_bottom}'),
            ]

        # Add custom legend
        ax_upper.legend(handles=custom_dots, title='Legend', loc="best", fontsize=10)
    
    ax_lower.set_xlabel("Chromosome")
    ax_upper.set_xlabel("")
    
    # save ad show the plot
    plt.savefig(os.path.join(plots_dir, save_name), dpi=dpi)
    plt.show()

    return True
