"""
This module provides functions for generating various plots for GWAS (Genome-Wide Association Studies) data. The plots include QQ plots, beta-beta scatter plots, and trumpet plots for visualizing the effect sizes and power of GWAS results.

Functions:
----------
- qqplot_draw(df_gwas: pd.DataFrame, plots_dir: str, conf_color="lightgray", save_name: str='qq_plot.jpeg') -> bool:
    Draws a QQ plot for GWAS data.
- confidence_interval(n: int, conf_points: int=1500, conf_alpha: float=0.05) -> np.ndarray:
    Computes confidence intervals for the QQ plot.
- beta_beta_draw(gwas_1: pd.DataFrame, gwas_2: pd.DataFrame, p_col: str, beta_col: str, se_col: str, snp_col: str, label_1: str, label_2: str, plot_dir: str, significance: float=5e-8, annotate_coincidents: bool=True, save_name: str='beta_beta.jpeg', draw_error_line: bool=True, draw_reg_line: bool=True) -> bool:
- new_trumpet(df_gwas: pd.DataFrame, df_freq: pd.DataFrame, plot_dir: pd.DataFrame, snp_col: str, chr_col: str, pos_col: str, maf_col: str, beta_col: str, power_ts: list, n_case: int, n_control: int, sample_size: int=None, n_col: str='', sample_size_strategy: str='median', p_col: str=None, prevalence: int=None, mode: str='binary', p_filter: float=5e-8, to_highlight: list=[], to_annotate: list=[], cmap: str="cool", power_sig_level: float=5e-8, build='38', gtf_path: str=None, save_name: str='trumpet_plot.jpeg') -> bool:
    Generates a trumpet plot for visualizing the effect sizes and power of GWAS results.
"""

import os
import gzip
import shutil

import matplotlib
import matplotlib.colors as mc
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import textalloc as ta

from matplotlib.collections import LineCollection

from gwaslab.bd_download import download_file
from gwaslab.g_Log import Log
from gwaslab.util_in_get_sig import annogene
from gwaslab.util_in_calculate_power import get_beta
from gwaslab.util_in_calculate_power import get_beta_binary

def qqplot_draw(df_gwas:pd.DataFrame, plots_dir:str, conf_color="lightgray", save_name:str='qq_plot.jpeg')->bool:

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

    # compute confidence intervals
    mpts = confidence_interval(n, conf_points=1500, conf_alpha=0.05)

    # Plot the confidence interval as a filled polygon
    plt.fill(mpts[:, 0], mpts[:, 1], color=conf_color)

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
    plt.savefig(os.path.join(plots_dir, save_name), dpi=500)
    plt.show()

    return True

def confidence_interval(n:int, conf_points:int=1500, conf_alpha:float=0.05)->np.ndarray:

    """
    Function to confidence intervals for the QQ plot.

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
    ndarray
    """

    if not isinstance(n, int):
        raise ValueError(f"n must be an integer.")
    elif n < 0:
        raise ValueError(f"n must be positive.")
    
    if not isinstance(conf_points, int):
        raise ValueError(f"conf_points must be an integer.")
    elif conf_points < 0:
        raise ValueError(f"conf_points must be positive.")
    
    if not isinstance(conf_alpha, float):
        raise ValueError(f"conf_alpha must be a float.")
    elif conf_alpha < 0 or conf_alpha > 1:
        raise ValueError(f"conf_alpha must be between 0 and 1.")
    
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
    
    return mpts

def beta_beta_draw(gwas_1:pd.DataFrame, gwas_2:pd.DataFrame, p_col:str, beta_col:str, se_col:str, snp_col:str, label_1:str, label_2:str, plot_dir:str, significance:float=5e-8, annotate_coincidents:bool=True, save_name:str='beta_beta.jpeg', draw_error_line:bool=True, draw_reg_line:bool=True)->bool:
    
    """
    Generates a scatter plot comparing the effect sizes (beta values) of two GWAS studies.

    Parameters:
    -----------
    gwas_1 : pd.DataFrame
        DataFrame containing the first GWAS data.
    gwas_2 : pd.DataFrame
        DataFrame containing the second GWAS data.
    p_col : str
        Column name for p-values in the GWAS dataframes.
    beta_col : str
        Column name for beta values in the GWAS dataframes.
    se_col : str
        Column name for standard errors in the GWAS dataframes.
    snp_col : str
        Column name for SNP identifiers in the GWAS dataframes.
    label_1 : str
        Label for the first GWAS study.
    label_2 : str
        Label for the second GWAS study.
    plot_dir : str
        Directory where the plot will be saved.
    significance : float, optional
        Significance threshold for p-values, by default 5e-8.
    annotate_coincidents : bool, optional
        Whether to annotate SNPs that are significant in both GWAS studies, by default True.
    save_name : str, optional
        Name of the saved plot file, by default 'beta_beta.jpeg'.
    draw_error_line : bool, optional
        Whether to draw error bars, by default True.
    draw_reg_line : bool, optional
        Whether to draw a regression line, by default True.

    Returns:
    --------
    bool
        True if the plot is successfully generated and saved.

    Raises:
    -------
    ValueError
        If the input dataframes are not pandas dataframes.
        If the specified columns are not present in the dataframes.
        If the significance level is not a float between 0 and 1.
        If the boolean parameters are not boolean.
    """

    
    
    # check if the dataframes are pandas dataframes
    if not isinstance(gwas_1, pd.DataFrame):
        raise ValueError(f"GWAS 1 dataframe must be a pandas dataframe.")
    if not isinstance(gwas_2, pd.DataFrame):
        raise ValueError(f"GWAS 2 dataframe must be a pandas dataframe.")
    
    # check if the column names are in the dataframe
    if beta_col not in gwas_1.columns or p_col not in gwas_2.columns:
        raise ValueError(f"Column {beta_col} not present in both GWAS dataframes.")
    if p_col not in gwas_1.columns or p_col not in gwas_2.columns:
        raise ValueError(f"Column {p_col} not present in both GWAS dataframes.")
    
    # check if the significance level is a float between 0 and 1
    if not isinstance(significance, float):
        raise ValueError(f"Significance level must be a float.")
    elif significance < 0 or significance > 1:
        raise ValueError(f"Significance level must be between 0 and 1.")
    
    # check if the boolean values are boolean
    if not isinstance(annotate_coincidents, bool):
        raise ValueError(f"annotate_coincidents must be a boolean.")
    if not isinstance(draw_error_line, bool):
        raise ValueError(f"draw_error_line must be a boolean.")
    if not isinstance(draw_reg_line, bool):
        raise ValueError(f"draw_reg_line must be a boolean.")
    
    # check if the se column in dataframes if draw_error_line is True
    if draw_error_line is True:
        if se_col not in gwas_1.columns or se_col not in gwas_2.columns:
            raise ValueError(f"Column {se_col} not present in both GWAS dataframes.")

    # check if the snp column is in the dataframes if annotate_coincidents is True
    if annotate_coincidents is True:
        if snp_col not in gwas_1.columns or snp_col not in gwas_2.columns:
            raise ValueError(f"Column {snp_col} not present in both GWAS dataframes.")
    
    # rename columns to avoid conflicts
    df_gwas1 = gwas_1.copy()
    df_gwas1.columns = [f"{col}_1" for col in df_gwas1.columns if col != snp_col]
    df_gwas2 = gwas_2.copy()
    df_gwas2.columns = [f"{col}_2" for col in df_gwas2.columns if col != snp_col]

    # merge the dataframes
    df = pd.merge(df_gwas1, df_gwas2, on=snp_col, how='inner')

    del df_gwas1, df_gwas2

    # create masks to identify the SNPs that are significant in one or both GWAS
    mask_significance_1 = (df[f'{p_col}_1'] < significance)
    mask_significance_2 = (df[f'{p_col}_2'] < significance)

    on_first = df[mask_significance_1 &  ~mask_significance_2].reset_index(drop=True)[snp_col].to_list()
    on_both  = df[mask_significance_1 &  mask_significance_2].reset_index(drop=True)[snp_col].to_list()
    on_second= df[mask_significance_2 &  ~mask_significance_1].reset_index(drop=True)[snp_col].to_list()

    df[f'P-val<{significance}']= None

    df.loc[df['ID'].isin(on_both), f'P-val<{significance}']  = 'Both'
    df.loc[df['ID'].isin(on_first), f'P-val<{significance}'] = f'{label_1} GWAS'
    df.loc[df['ID'].isin(on_second), f'P-val<{significance}']= f'{label_2} GWAS'

    # set the limits for the plot
    max_beta_x = df[f'{beta_col}_1'].abs().max() + 0.01
    max_beta_y = df[f'{beta_col}_2'].abs().max() + 0.01
    max_coords = max(max_beta_x, max_beta_y)

    x_lim = (-max_coords, max_coords)
    y_lim = (-max_coords, max_coords)

    # determine colors and markers for the plot
    colors = {
        f'{label_1} GWAS': "#1f77b4", 
        'Both'           : "#ff7f0e",
        f'{label_2} GWAS': "#2ca02c"
    }
    markers= {
        f'{label_1} GWAS': 'o', 
        'Both'           : 's',
        f'{label_2} GWAS': '^'
    }

    fig= plt.figure(figsize=(10, 10))
    ax = plt.subplot(111)

    # plot with error bars
    if draw_error_line:

        for category in np.unique(df[f'P-val<{significance}']):

            mask_hue = (df[f'P-val<{significance}'] == category)

            ax.errorbar(
                x         =df[f'{beta_col}_1'][mask_hue], 
                y         =df[f'{beta_col}_2'][mask_hue], 
                xerr      =df[f'{se_col}_1'][mask_hue], 
                yerr      =df[f'{se_col}_2'][mask_hue],
                fmt       =markers[category],
                color     =colors[category], 
                ecolor    ='lightgray', 
                label     =category, 
                elinewidth=0.75, 
                capsize   =0,
                markersize=5
            )
    # plot without error bars
    else:
        ax = sns.scatterplot(
            data=df,
            x=f'{beta_col}_1',
            y=f'{beta_col}_2',
            hue=f'P-val<{significance}',
            palette=colors,
            style=f'P-val<{significance}',
            markers=markers,
            s=50,
            edgecolor='black',
            alpha=0.8
        )

    # draw x and y axis
    ax.axhline(0, color='black', linestyle='solid', lw=0.5)
    ax.axvline(0, color='black', linestyle='solid', lw=0.5)
    
    # draw y=x line
    help_line = np.linspace(-max_beta_x, max_beta_x, 100)
    ax.plot(help_line, help_line, color='black', linestyle='solid', lw=0.5)
    
    # draw regression line
    if draw_reg_line:
    
        result = stats.linregress(df[f'{beta_col}_2'], df[f'{beta_col}_2'])

        ax.plot(
            help_line, 
            result.slope*help_line + result.intercept, 
            color    ='gray', 
            linestyle='dashed', 
            lw       =0.5
        )

    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    
    ax.set_xlabel(f'Per Allele Effect Size {label_1} GWAS', fontsize=7)
    ax.set_ylabel(f'Per Allele Effect Size {label_2} GWAS', fontsize=7)
    
    ax.tick_params(axis='both', which='major', labelsize=7)
    
    plt.legend(loc='best', fontsize=7, title=f'P-value < {significance}')

    plt.tight_layout()

    # render and draw the figure without showing it 
    #r = fig.canvas.get_renderer()
    fig.canvas.draw()

    if annotate_coincidents:
        to_annotate = df[df[f'P-val<{significance}']=='Both'].reset_index(drop=True)
        
        texts=[]
        text_x = []
        text_y = []

        for i, row in to_annotate.iterrows():

            texts.append(row[snp_col])
            text_x.append(row[f'{beta_col}_1'])
            text_y.append(row[f'{beta_col}_2'])

        ta.allocate(
                ax,
                x        =text_x,
                y        =text_y,
                text_list=texts,
                x_scatter=df[f'{beta_col}_1'].to_list(),
                y_scatter=df[f'{beta_col}_2'].to_list(),
                linecolor='black',
                textsize =7,
                bbox     =dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.8),
            )
        
    plt.savefig(os.path.join(plot_dir, save_name), dpi=500)
    plt.show()

    return True
    
def trumpet_draw(df_gwas:pd.DataFrame, df_freq:pd.DataFrame, plot_dir:pd.DataFrame, snp_col:str, chr_col:str, pos_col:str, maf_col:str, beta_col:str, power_ts:list, n_case:int, n_control:int, sample_size:int=None, n_col:str='', sample_size_strategy:str='median', p_col:str=None, prevalence:int=None, mode:str='binary', p_filter:float=5e-8, to_highlight:list=[], to_annotate:list=[], cmap:str="cool", power_sig_level:float=5e-8, build='38', gtf_path:str=None, save_name:str='trumpet_plot.jpeg')->bool:

    if not isinstance(df_gwas, pd.DataFrame):
        raise ValueError(f"GWAS dataframe must be a pandas dataframe.")

    # check if the column names are in the dataframe
    if maf_col not in df_gwas.columns:
        if df_freq is None:
            raise ValueError(f"Column {maf_col} not present in the GWAS dataframe and no frequency dataframe provided.")
        elif not isinstance(df_freq, pd.DataFrame):
            raise ValueError(f"Frequency dataframe must be a pandas dataframe.")
        elif maf_col not in df_freq.columns:
            raise ValueError(f"Column {maf_col} not present in the GWAS dataframe and frequency dataframe.")
        elif df_freq.shape[0]==0:
            raise ValueError(f"Frequency dataframe is empty.")
        elif snp_col not in df_freq.columns or snp_col not in df_gwas.columns:
            raise ValueError(f"Column {snp_col} is not a common column. No possible merge.")
    
    if df_gwas.shape[0]==0:
        raise ValueError(f"GWAS dataframe is empty.")
    if df_gwas is None:
        raise ValueError(f"GWAS dataframe is None.")
    if not isinstance(power_ts, list):
        raise ValueError(f"Power thresholds must be a list.")
    for val in power_ts:
        if not isinstance(val, float):
            raise ValueError(f"Power thresholds must be floats.")
        elif val > 1 or val < 0:
            raise ValueError(f"Power thresholds must be between 0 and 1.")
    if mode not in ['binary', 'quantitative']:
        raise ValueError(f"Mode must be either 'binary' or 'quantitative'.")
    
    if mode == 'binary':
        if n_case is None or n_control is None:
            raise ValueError(f"Number of cases and controls must be provided.")
        elif not isinstance(n_case, int) or not isinstance(n_control, int):
            raise ValueError(f"Number of cases and controls must be integers.")
        elif n_case < 0 or n_control < 0:
            raise ValueError(f"Number of cases and controls must be positive.")
        if prevalence is None:
            prevalence = n_case / (n_case + n_control)

    if mode == 'quantitative':
        if n_col not in df_gwas.columns and sample_size is None:
            raise ValueError(f"Column {n_col} not present in the GWAS dataframe and sample size is unknown.")
        elif sample_size is not None and not isinstance(sample_size, int):
            raise ValueError(f"Sample size must be an integer.")
        elif sample_size < 0:
            raise ValueError(f"Sample size must be positive.")
        
        if n_col in df_gwas.columns:
            if sample_size_strategy not in ['min', 'max', 'median', 'mean']:
                raise ValueError(f"Sample size strategy must be either 'min', 'max', 'median' or 'mean'.")


    if p_filter is not None and p_col is not None:
        if not isinstance(p_filter, float):
            raise ValueError(f"Significance level must be a float.")
        elif p_filter < 0 or p_filter > 1:
            raise ValueError(f"Significance level must be between 0 and 1.")
        elif p_col not in df_gwas.columns:
            raise ValueError(f"Column {p_col} not present in the GWAS dataframe.")
        else:
            # filter the dataframe according to the significance level
            gwas_df = df_gwas[df_gwas[p_col] < p_filter].reset_index(drop=True)
    else:
        gwas_df = df_gwas.copy()

    if not isinstance(to_highlight, list):
        raise ValueError(f"to_highlight must be a list.")
    if not isinstance(to_annotate, list):
        raise ValueError(f"to_annotate must be a list.")
    for val in to_highlight:
        if not isinstance(val, str):
            raise ValueError(f"to_highlight must be a list of strings.")
    for val in to_annotate:
        if not isinstance(val, str):
            raise ValueError(f"to_annotate must be a list of strings.")

            
    if maf_col not in gwas_df.columns:
        df = pd.merge(gwas_df, df_freq, on=snp_col, how='inner')
    
    del gwas_df, df_freq

    # colormap
    cmap_to_use = matplotlib.colormaps[cmap]
    if cmap_to_use.N >100:
        rgba = cmap_to_use(power_ts)
    else:
        rgba = cmap_to_use(range(len(power_ts)))

    output_hex_colors=[]
    for i in range(len(rgba)):
        output_hex_colors.append(mc.to_hex(rgba[i]))

    # generate figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # BETA and MAF range for power calculation
    maf_min_power = np.floor( -np.log10(df[maf_col].min())) + 1
    maf_range=(min(np.power(10.0,-maf_min_power),np.power(10.0,-4)),0.5)
    
    if df[beta_col].max()>3:
        beta_range=(0.0001, df[beta_col].max())
    else:
        beta_range=(0.0001, 3)

    if mode=="quantitative" and sample_size is None:
        if sample_size_strategy == "min":
            sample_size = df[n_col].min() 
        elif sample_size_strategy == "max":
            sample_size = df[n_col].max() 
        elif sample_size_strategy == "mean":
            sample_size = df[n_col].mean() 
        else:
            sample_size = df[n_col].median()

    # generate power lines
    if mode=='binary':
        for i,t in enumerate(power_ts):
            xpower = get_beta_binary(
                eaf_range =maf_range,
                beta_range=beta_range, 
                prevalence=prevalence,
                or_to_rr  =False,
                ncase     =n_case, 
                ncontrol  =n_control, 
                t         =t,
                sig_level =power_sig_level,
                n_matrix  =2000
            )

            xpower2   = xpower.copy()
            xpower2[1]= -xpower2[1] 
            xpower2[1]= xpower2[1] 
            xpower[1] = xpower[1] 
            lines     = LineCollection([xpower2,xpower], label=t,color=output_hex_colors[i])

            ax.add_collection(lines)

    if mode=='quantitative':
        for i,t in enumerate(power_ts):
            xpower = get_beta(
                mode      ="q",          
                eaf_range =maf_range,
                beta_range=beta_range, 
                n         =sample_size,
                t         =t,
                sig_level =power_sig_level,
                n_matrix  =2000
            )
            
            xpower2   = xpower.copy()
            xpower2[1]= -xpower2[1]
            xpower2[1]= xpower2[1]
            xpower[1] = xpower[1]

            lines = LineCollection([xpower2,xpower], label=t,color=output_hex_colors[i],zorder=0)

            ax.add_collection(lines)

    # get absolute value of BETA for scaling
    df["ABS_BETA"] = df[beta_col].abs()
    size_norm = (df["ABS_BETA"].min(), df["ABS_BETA"].max())

    # scatter plot
    ax = sns.scatterplot(data=df,
                    x=maf_col,
                    y=beta_col,
                    size     ="ABS_BETA", 
                    ax       =ax, 
                    sizes    =(20,80),
                    size_norm=size_norm,
                    legend   =True, 
                    edgecolor="black",
                    alpha    =0.8,
                    zorder   =2,
    )

    # add legend for power lines
    h,l = ax.get_legend_handles_labels()
    if len(power_ts)>0:
        l1 = ax.legend(
            h[:int(len(power_ts))],
            l[:int(len(power_ts))], 
            title         ="Power", 
            loc           ="upper right",
            fontsize      =7,
            title_fontsize=7
        )
        for line in l1.get_lines():
            line.set_linewidth(5.0)
        ax.add_artist(l1)

    # add legend for size
    l2 = ax.legend(
        h[int(len(power_ts)):],
        l[int(len(power_ts)):], 
        title="ABS_BETA", 
        loc="lower right",
        fontsize =7, 
        title_fontsize=7
    )

    if to_highlight is not None and len(to_highlight)>0:

        df_high = df[df[snp_col].isin(to_highlight)].reset_index(drop=True)

        ax = sns.scatterplot(
            data=df_high,
            x=maf_col,
            y=beta_col,
            size     ="ABS_BETA", 
            ax       =ax, 
            sizes    =(20,80),
            size_norm=size_norm,
            legend   =False,
            color = 'red',
            edgecolor="black",
            alpha    =0.8,
            zorder   =2,
            )

    ax.tick_params(axis='y', labelsize=7)

    # add X axis
    ax.axhline(y=0,color="black",linestyle="solid", lw=0.5)

    # axis limits
    ylim =(-df[beta_col].abs().max()*1.5, df[beta_col].abs().max()*1.5)
    xlim =(df[maf_col].min()*0.5,0.52)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    plt.tight_layout()

    r = fig.canvas.get_renderer()
    fig.canvas.draw()

    if to_annotate is not None and len(to_annotate)>0:

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

        variants_toanno = df[df[snp_col].isin(to_annotate)].reset_index(drop=True)
        
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

        split = {
            'top':variants_toanno[variants_toanno[beta_col]>0].reset_index(drop=True), 
            'bottom':variants_toanno[variants_toanno[beta_col]<0].reset_index(drop=True)
        }

        for key in split.keys():
            
            texts = []  # a list to store text annotations for adjustment
            x = []      # a list to store x-coordinates for adjustment
            y = []      # a list to store y-coordinates for adjustment

            for i, row in split[key].iterrows():

                x.append(row[maf_col])
                y.append(row[beta_col])
                texts.append(row['GENENAME'])
            
            if key == 'top':
                direction = 'northeast'
            else:
                direction = 'southeast'

            ta.allocate(
                    ax,              # the axis to which the text will be
                    x        =x,     # x-coordinates of the data point to annotate
                    y        =y,     # y-coordinates of the data point to annotate
                    text_list=texts, # list of text to annotate
                    x_scatter=df[maf_col], # all scatter points x-coordinates
                    y_scatter=df[beta_col],  # all scatter points y-coordinates
                    linecolor='black',                      # color of the line connecting the text to the data point
                    textsize =7,                            # size of the text (Default to Nature standard)
                    bbox     =dict(boxstyle='round,pad=0.5', edgecolor='black', facecolor='#f0f0f0', alpha=0.5),
                    avoid_label_lines_overlap =True,
                    avoid_crossing_label_lines=False,
                    min_distance     =0.05,
                    max_distance     =0.5,
                    margin           =0.01,
                    rotation         =90,
                    ha               ='center',
                    direction        =direction,
                    nbr_candidates   =300,
                    priority_strategy=42,
                    plot_kwargs      =dict(linestyle=':')
                )
    
    plt.show()
    plt.savefig(os.path.join(plot_dir, save_name), dpi=500)
    
    return True
