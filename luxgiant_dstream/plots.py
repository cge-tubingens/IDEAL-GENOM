"""
Module with function for drawing plots

The module provides functions to draw manhattan and qq plots for GWAS data.

Functions:
----------
manhattan_plot(df_gwas:pd.DataFrame, plots_dir:str, df_annot:pd.DataFrame=None, annotate:bool=False)->bool
    Draw a manhattan plot for GWAS data.
qq_plot(df_gwas:pd.DataFrame, plots_dir:str)->bool
    Draw a qq plot for GWAS data.

The function qq_plot is a Python implementation of the R code provided in the following links: https://gist.github.com/MrFlick/10477946 or here https://github.com/hakyimlab/IntroStatGen/blob/master/gists/qqunif.r

"""

import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from matplotlib.axes import Axes
from pandas.core.groupby.generic import DataFrameGroupBy

from adjustText import adjust_text

def manhattan_plot(df_gwas:pd.DataFrame, plots_dir:str=None, df_annot:pd.DataFrame=None, annotate:bool=False)->bool:

    """
    Function to draw a manhattan plot for GWAS data.
    
    Parameters:
    -----------
    df_gwas : pd.DataFrame
        Dataframe with GWAS summary statistics.
    plots_dir : str
        Path to the directory to save the plot.
    df_annot : pd.DataFrame (default=None)
        Dataframe with annotation information for SNPs of interest.
    annotate : bool (default=False)
        Whether to annotate the plot with gene names.

    Returns:
    --------
    bool
    """
    
    # keep columns of interest
    df = df_gwas.copy()
    df['log10P'] = -np.log10(df['p'])

    # sort values by chromosome
    df = df.sort_values(['CHR', 'bp'])

    # to get colors by chromosome
    df['ind'] = range(len(df))
    df_grouped = df.groupby(('CHR'))

    # Subset dataframe to highlight specific SNPs
    if df_annot is not None:
        snps = df_annot['SNP'].to_list()
        genes = df_annot['GENE'].to_list()
        highlighted_snps = df[df['SNP'].isin(snps)]  # Filter for the SNPs of interest

    # Set the figure size
    fig = plt.figure(figsize=(50, 25))
    ax = fig.add_subplot(111)

    if df_annot is not None:
        ax = draw_chrom_groups(ax, df_grouped, highlighted_snps)
    else:
        ax = draw_chrom_groups(ax, df_grouped, None)

    if annotate:
        # Add gene names to highlighted SNPs
        ax = annotate_plot(ax, highlighted_snps, genes)

    # Set axis limits
    ax.set_xlim([0, len(df)])
    ax.set_ylim([0, max(df['log10P']) + 1])

    # Add horizontal line for genome-wide significance threshold (-log10(5e-8))
    ax.axhline(y=-np.log10(5e-8), color='blue', linestyle='--', linewidth=1)

    # Add labels and title
    ax.set_xlabel('Chromosome', fontsize=24)  # Adjust fontsize as needed
    ax.set_ylabel('-log10(P-value)', fontsize=24)  # Adjust fontsize as needed
    ax.set_title('Manhattan Plot', fontsize=30)  # Add a title if needed

    # Customize the grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    if plots_dir is None:
        return ax
    else:
        # Display the plot
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'manhattan_plot.png'))

        return True
    
def draw_chrom_groups(ax:Axes, df_chr_group:DataFrameGroupBy, highlighted_snps:pd.DataFrame)->Axes:

    """
    Plots chromosome groups on a given Axes object with alternating colors and highlights specific SNPs.
    
    Parameters:
    -----------
    ax (Axes): 
        The matplotlib Axes object where the plot will be drawn.
    df_chr_group (DataFrameGroupBy): 
        Grouped DataFrame by chromosome containing 'ind' and 'log10P' columns.
    highlighted_snps (pd.DataFrame): 
        DataFrame containing SNPs to be highlighted with 'ind' and 'log10P' columns.
    
    Returns:
    --------
    Axes: The Axes object with the plotted chromosome groups and highlighted SNPs.
    
    Notes:
    - Chromosome groups are plotted with alternating colors (grey and skyblue).
    - Highlighted SNPs are plotted in red with a larger point size.
    - X-axis labels are set to chromosome names with positions calculated as the middle of each chromosome group.
    """

    # Colors for alternating chromosomes
    colors = ['grey', 'skyblue']

    # Variables to store labels and positions for the x-axis
    x_labels = []
    x_labels_pos = []

    # Loop over each group by chromosome
    for num, (name, group) in enumerate(df_chr_group):

        # Plot each chromosome with alternating colors
        ax.scatter(group['ind'], group['log10P'], color=colors[num % len(colors)], s=10)

        # Add chromosome label and its position
        x_labels.append(name)
        middle_pos = (group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2
        x_labels_pos.append(middle_pos)

        # Set the x-axis labels and positions
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels, fontsize=20)  # Adjust fontsize as needed

        if highlighted_snps is not None:
            # Plot highlighted SNPs with a different color (red) and larger point size
            ax.scatter(highlighted_snps['ind'], highlighted_snps['log10P'], color='red', s=50, label='Highlighted SNPs')

    return ax

def annotate_plot(ax:Axes, highlighted_snps:pd.DataFrame, snps:list, genes:list)->Axes:
    
    """
    Annotates a plot with gene names corresponding to highlighted SNPs.

    Parameters:
    ----------
    ax (Axes): 
        The matplotlib Axes object to annotate.
    highlighted_snps (pd.DataFrame): 
        DataFrame containing SNPs to be highlighted with columns 'SNP', 'ind', and 'log10P'.
    snps (list): 
        List of SNP identifiers.
    genes (list): 
        List of gene names corresponding to the SNPs.
    
    Returns:
    --------
    Axes: The annotated matplotlib Axes object.
    """



    texts = []  # A list to store text annotations for adjustment
    for i, row in highlighted_snps.iterrows():
        gene = genes[snps.index(row['SNP'])]  # Get corresponding gene name
        # Add text label to the SNP
        text = ax.text(row['ind'], row['log10P'], gene, fontsize=12, ha='right', va='bottom', color='black',
                       bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))
        texts.append(text)
    # Adjust the text to prevent overlaps using adjustText
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'))

    return ax

def qq_plot(df_gwas:pd.DataFrame, plots_dir:str)->bool:

    """
    Function to draw a qq plot for GWAS data.

    Parameters:
    -----------
    df_gwas : pd.DataFrame
        Dataframe with GWAS summary statistics.
    plots_dir : str
        Path to the directory to save the plot.

    Returns:
    --------
    bool
    """
    
    pvalues = df_gwas['p'].values
    grp = None
    n = len(pvalues)

    log_p = -np.log10(pvalues)
    exp_x = -np.log10((stats.rankdata(pvalues, method='ordinal') - 0.5) / n)

    if grp is not None:

        # Create a DataFrame with rounded pvalues and exp_x, and grp
        thin = pd.DataFrame({
            'pvalues': np.round(log_p, 3),
            'exp_x': np.round(exp_x, 3),
            'grp': grp
        }).drop_duplicates()

        # Assign the updated group after thinning
        grp = thin['grp'].values
    else:

        # Create a DataFrame with rounded pvalues and exp_x
        thin = pd.DataFrame({
            'pvalues': np.round(log_p, 3),
            'exp_x': np.round(exp_x, 3)
        }).drop_duplicates()

    # Update pvalues and exp_x after thinning
    log_p = thin['pvalues'].values
    exp_x = thin['exp_x'].values

    axis_range =  [float(min(log_p.min(), exp_x.min()))-0.5, float(max(log_p.max(), exp_x.max()))+1]

    fig, ax = plt.subplots(figsize=(10,10))

    # Calculate the confidence intervals if draw_conf is True
    def plot_confidence_interval(n:int, conf_points:int=1500, conf_col:str="lightgray", conf_alpha:float=0.05)->None:

        """
        Function to plot the confidence interval for the QQ plot.

        Parameters:
        -----------
        n : int
            Number of p-values.
        conf_points : int (default=1500)
            Number of points to plot the confidence interval.
        conf_col : str (default="lightgray")
            Color of the confidence interval.
        conf_alpha : float (default=0.05)
            Alpha value for the confidence interval.

        Returns:
        --------
        None
        """
        
        conf_points = min(conf_points, n - 1)
        mpts = np.zeros((conf_points * 2, 2))

        for i in range(1, conf_points + 1):

            x = -np.log10((i - 0.5) / n)

            y_upper = -np.log10(stats.beta.ppf(1 - conf_alpha / 2, i, n - i))
            y_lower = -np.log10(stats.beta.ppf(conf_alpha / 2, i, n - i))

            mpts[i - 1, 0] = x
            mpts[i - 1, 1] = y_upper
            mpts[conf_points * 2 - i, 0] = x
            mpts[conf_points * 2 - i, 1] = y_lower

        # Plot the confidence interval as a filled polygon
        plt.fill(mpts[:, 0], mpts[:, 1], color=conf_col)
        pass

    plot_confidence_interval(n, conf_points=1500, conf_col="lightgray", conf_alpha=0.05)

    if grp is not None:
        unique_groups = np.unique(grp)
        for group in unique_groups:
            group_mask = (grp == group)
            ax.scatter(exp_x[group_mask], log_p[group_mask], label=f'Group {group}', marker='o')
    else:
        ax.scatter(exp_x, log_p, marker='o')

    # Line y = x
    ax.plot(axis_range, axis_range, color='red', linestyle='--', lw=1)

    ax.set_xlim(axis_range)
    ax.set_ylim(axis_range)

    ax.set_xlabel('Expected (-log10 p-value)')
    ax.set_ylabel('Observed (-log10 p-value)')
    ax.set_aspect('equal')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'qq_plot.png'))

    return True

def miami_plot(df_top:pd.DataFrame, df_bottom:pd.DataFrame, plots_dir:str)->bool:

    chr_colors           = ['grey', 'skyblue']
    upper_ylab           = "-log10(p)" 
    lower_ylab           = "-log10(p)" 
    genome_line          = 5e-8
    genome_line_color    = "red"
    suggestive_line      = 1e-5 
    suggestive_line_color= "blue"

    # format data to draw miami plot
    plot_data = process_miami_data(df_top, df_bottom)

    # Set axis labels for upper and lower plot
    def format_ylabel(label):
        return f"{label}\n-log10(p)" if label != "-log10(p)" else r"-log10(p)"
    
    upper_ylab = format_ylabel(upper_ylab)
    lower_ylab = format_ylabel(lower_ylab)

    # Create the figure
    plt.figure(figsize=(20, 16))

    ax_upper = plt.subplot(211)
    sns.scatterplot(x=plot_data['upper']['rel_pos'], y=plot_data['upper']['log10p'],
                    hue=plot_data['upper']['CHR'], palette=chr_colors, ax=ax_upper, s=1)
    ax_upper.set_ylabel(upper_ylab)
    ax_upper.set_xlim(0, plot_data['axis']['center'].max())

    ax_upper.set_xlabel("")

    x_ticks=plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis']['CHR'].astype(str).tolist()

    ax_upper.set_xticks(ticks=x_ticks)  # Set x-ticks
    ax_upper.set_xticklabels(x_labels)

    # Customize x-axis labels
    #ax_upper.set_xticks(ticks=plot_data['axis']['center'], labels=plot_data['axis']['CHR'])  # Set x-ticks
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_upper.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.5)
    
    if genome_line is not None:
        ax_upper.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.5)
    
    ax_upper.legend([], frameon=False)  # Remove legend

    # Create the lower plot
    ax_lower = plt.subplot(212)
    sns.scatterplot(x=plot_data['lower']['rel_pos'], y=plot_data['lower']['log10p'],
                    hue=plot_data['lower']['CHR'], palette=chr_colors, ax=ax_lower, s=1)
    ax_lower.set_ylabel(lower_ylab)
    ax_lower.set_ylim(plot_data['maxp'], 0)  # Reverse y-axis
    ax_lower.set_xlim(0, plot_data['axis']['center'].max())

    ax_lower.set_xlabel("")

    ax_lower.set_xticks(ticks=x_ticks) # Set x-ticks
    ax_lower.set_xticklabels([])
    ax_lower.xaxis.set_ticks_position('top')
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_lower.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.5)
    
    if genome_line is not None:
        ax_lower.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.5)
    
    ax_lower.legend([], frameon=False)  # Remove legend
    
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'miami_plot.png'))
    plt.show()

    return True

def prepare_data(data_top:pd.DataFrame, data_bottom:pd.DataFrame)->pd.DataFrame:

    data_top['split_by'] = 'top'
    data_bottom['split_by'] = 'bottom'

    joint = pd.concat([data_top, data_bottom], axis=0)

    return joint

def compute_relative_pos(data:pd.DataFrame, chr_col:str='CHR', pos_col:str='bp', p_col:str='p')->pd.DataFrame:

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

    chromosomes = data[chr_col].unique()

    axis_center = pd.DataFrame(columns=['CHR', 'center'])

    for i, chrom in enumerate(chromosomes):

        temp = data[data[chr_col] == chrom].reset_index(drop=True)

        axis_center.loc[i, 'CHR'] = chrom
        axis_center.loc[i, 'center'] = np.round((temp[chr_pos_col].max()+temp[chr_pos_col].min())/2,0)

    return axis_center

def process_miami_data(data_top:pd.DataFrame, data_bottom:pd.DataFrame)->dict:

    data = prepare_data(data_top, data_bottom)

    data = compute_relative_pos(data, chr_col='CHR', pos_col='bp', p_col='p')

    axis_center = find_chromosomes_center(data)

    maxp = np.ceil(data['log10p'].max(skipna=True))

    df_top = data[data['split_by'] == 'top']
    df_bottom = data[data['split_by'] == 'bottom']

    miami_data = {
        'upper': df_top,
        'lower': df_bottom,
        'axis': axis_center,
        'maxp': maxp
    }

    return miami_data
