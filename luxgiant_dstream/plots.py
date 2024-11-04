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
import gzip
import shutil

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import textalloc as ta

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.colors as mc
import matplotlib
from adjustText import adjust_text
from gwaslab.bd_download import download_file
from gwaslab.util_in_get_sig import annogene
from gwaslab.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz_plot_mqqplot import _process_highlight
from gwaslab.g_Log import Log
from gwaslab.util_in_calculate_power import get_beta
from gwaslab.util_in_calculate_power import get_beta_binary
from gwaslab.util_in_fill_data import filldata

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

def beta_beta_draw(gwas_1:pd.DataFrame, gwas_2:pd.DataFrame, p_col:str, beta_col:str, se_col:str, snp_col:str, label_1:str, label_2:str, plot_dir:str, significance:float=5e-8, annotate_coincidents:bool=True, save_name:str='beta_beta.jpeg')->bool:

    df_gwas1 = gwas_1.copy()
    df_gwas1.columns = [f"{col}_1" for col in df_gwas1.columns if col != snp_col]
    df_gwas2 = gwas_2.copy()
    df_gwas2.columns = [f"{col}_2" for col in df_gwas2.columns if col != snp_col]

    df = pd.merge(df_gwas1, df_gwas2, on=snp_col, how='inner')

    del df_gwas1, df_gwas2

    mask_significance_1 = (df[f'{p_col}_1'] < significance)
    mask_significance_2 = (df[f'{p_col}_2'] < significance)

    on_first = df[mask_significance_1 &  ~mask_significance_2].reset_index(drop=True)[snp_col].to_list()
    on_both  = df[mask_significance_1 &  mask_significance_2].reset_index(drop=True)[snp_col].to_list()
    on_second= df[mask_significance_2 &  ~mask_significance_1].reset_index(drop=True)[snp_col].to_list()

    df[f'P-val<{significance}']= None

    df.loc[df['ID'].isin(on_both), f'P-val<{significance}']  = 'Both'
    df.loc[df['ID'].isin(on_first), f'P-val<{significance}'] = f'{label_1} GWAS'
    df.loc[df['ID'].isin(on_second), f'P-val<{significance}']= f'{label_2} GWAS'

    max_beta_x = df[f'{beta_col}_1'].abs().max() + 0.01
    max_beta_y = df[f'{beta_col}_2'].abs().max() + 0.01
    max_coords = max(max_beta_x, max_beta_y)

    x_lim = (-max_coords, max_coords)
    y_lim = (-max_coords, max_coords)

    result = stats.linregress(df[f'{beta_col}_2'], df[f'{beta_col}_2'])

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

    # draw x and y axis
    ax.axhline(0, color='black', linestyle='solid', lw=0.5)
    ax.axvline(0, color='black', linestyle='solid', lw=0.5)
    
    # draw y=x line
    help_line = np.linspace(-max_beta_x, max_beta_x, 100)
    ax.plot(help_line, help_line, color='black', linestyle='solid', lw=0.5)
    
    # draw regression line
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

def trumpet_draw(df_gwas:pd.DataFrame,
                      plot_dir:str,
                snpid:str="SNPID",
                mode:str="q",
                chrom:str="CHR",
                pos:str="POS",
                n:str="N",
                p:str="P",
                maf:str="MAF",
                eaf:str="EAF",
                beta:str="BETA",
                ts:list=None,
                anno=None,
                prevalence:float=None,
                ncase:int=None,
                ncontrol:int=None, 
                sig_level:float=5e-8,
                p_level:float=5e-8,             
                maf_range=None,
                beta_range=None, 
                n_matrix:int=1000,
                xscale:str="log",
                yscale_factor=1,
                cmap:str="cool",
                ylim=None,
                xlim=None,
                markercolor="#597FBD",
                hue=None,
                highlight = None,
                highlight_chrpos = False,
                highlight_color="#CB132D",
                highlight_windowkb = 500,
                highlight_anno_args = None,
                scatter_args=None,
                fontsize:str=15,
                font_family:str="DejaVu Sans",
                size:str= "ABS_BETA",
                sizes=None,
                figargs:dict=None,
                build="99",
                anno_set=None,
                anno_alias=None,
                anno_d=None,
                anno_args=None,
                anno_style="expand",
                anno_source = "refseq",
                anno_max_iter=100,
                arm_scale=1,
                repel_force=0.01,
                ylabel="Effect size",
                xlabel="Minor allele frequency",
                xticks = None,
                xticklabels = None,
                yticks = None,
                yticklabels=None,
                sort:str="beta",
                verbose:bool=True,
                log=Log()):
    
    # code adapted from https://github.com/Cloufield/gwaslab/blob/main/src/gwaslab/viz_plot_trumpetplot.py
    
    #Checking columns#################################################################################################################
    matplotlib.rc('font', family=font_family)
    if sizes is None:
        sizes = (20,80)
    if anno_args is None:
        anno_args={"fontsize":12,"fontstyle":"italic"}
    if anno_set is None:
        anno_set=list()
    if anno_alias is None:
        anno_alias=dict()
    if anno_d is None:
        anno_d=dict()
    if ts is None:
        ts = [0.2,0.4,0.6,0.8]
    if xticks is None:
        if xscale== "log":
            xticks = [0.001,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks
        else:
            xticks = [0,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks            
    if figargs is None:
        figargs={"figsize":(10,8)}
    if scatter_args is None:
        scatter_args ={}
    if hue is not None:
        scatter_args["hue"]=hue
    if markercolor is not None:
        scatter_args["color"]=markercolor
    if highlight is None:
        highlight = list()
    #Checking columns#################################################################################################################
    log.write("Start to create trumpet plot...", verbose=verbose)
    
    #parameter check##################################################################################################################
    if (beta not in df_gwas.columns) or (eaf not in df_gwas.columns):
        log.write(" -No EAF or BETA columns. Skipping...", verbose=verbose)
        return None
    if mode=="b":
        if ncase is None or ncontrol is None:
            log.write(" -No scase or scontrol. Skipping...", verbose=verbose)
            return None
        if prevalence is None:
                prevalence= ncase/(ncase + ncontrol)
                log.write(" -Prevalence is not given. Estimating based on scase and scontrol :{}...".format(prevalence), verbose=verbose)
    
    #print settings##################################################################################################################

    log.write(" -Settings:", verbose=verbose)
    log.write("  -Mode: {}".format(mode), verbose=verbose)
    if mode == "q" :
        log.write("  -N: {}".format(n), verbose=verbose)
    if mode == "b" :
        log.write("  -N_CASE: {}".format(ncase), verbose=verbose)
        log.write("  -N_CONTROL: {}".format(ncontrol), verbose=verbose)
        log.write("  -PREVALENCE: {}".format(prevalence), verbose=verbose)
    log.write("  -BETA: {}".format(beta), verbose=verbose)
    log.write("  -Significance level: {}".format(sig_level), verbose=verbose)
    log.write("  -Power thresholds: {}".format(ts), verbose=verbose)
    log.write("  -Power line smoothness: {}".format(n_matrix), verbose=verbose)
    
    #loading columns #################################################################################################################
    cols_to_use = [snpid, beta, eaf, n, p]
    
    if len(highlight)>0: 
        cols_to_use.append(pos)
        cols_to_use.append(chrom)
    
    if anno is not None:
        if anno != "GENENAME":
            if anno!=True:
                log.write(" -Loading column {} for annotation...".format(anno), verbose=verbose)
                if anno not in cols_to_use:
                    cols_to_use.append(anno)
        else:
            cols_to_use.append(pos) if pos not in cols_to_use else cols_to_use
            cols_to_use.append(chrom) if chrom not in cols_to_use else cols_to_use

    if size != "ABS_BETA":
        if size not in cols_to_use:
            cols_to_use.append(size)
    if "hue" in scatter_args.keys():
        cols_to_use.append(scatter_args["hue"]) 
    #filter by p #################################################################################################################
    if p in df_gwas.columns:
        sumstats = df_gwas.loc[df_gwas[p]< p_level,cols_to_use ].copy()
        log.write(" -Excluding variants with P values > {}".format(p_level), verbose=verbose)
    else:
        cols_to_use.remove(p)
        sumstats = df_gwas[[beta,eaf,n]].copy()
    log.write(" -Plotting {} variants...".format(len(sumstats)), verbose=verbose)
    
    #add maf column #################################################################################################################
    if maf not in sumstats.columns:
        sumstats = filldata(sumstats,to_fill=["MAF"],verbose=False)
        is_filpped = (sumstats["MAF"] < sumstats[eaf]) & (sumstats[eaf] > 0.5)& (sumstats["MAF"] < 0.5)
        log.write(" -Flipping {} variants...".format(sum(is_filpped)), verbose=verbose)
        sumstats.loc[is_filpped, beta] = -sumstats.loc[is_filpped, beta]
    
    #configure n #################################################################################################################
    if mode=="q":
        if n == "N":
            n = sumstats["N"].median() 
        elif n == "max":
            n = sumstats["N"].max() 
        elif n == "min":
            n = sumstats["N"].min() 
        elif n == "median":
            n = sumstats["N"].median() 
        elif n == "mean":
            n = sumstats["N"].mean() 
        log.write(" -N for power calculation: {}".format(n), verbose=verbose)

    #configure beta and maf range ###################################################################################################
    if maf_range is None:
        maf_min_power = np.floor( -np.log10(sumstats[maf].min())) + 1
        maf_range=(min(np.power(10.0,-maf_min_power),np.power(10.0,-4)),0.5)
    if beta_range is None:
        if sumstats[beta].max()>3:
            beta_range=(0.0001,sumstats[beta].max())
        else:
            beta_range=(0.0001,3)
    
    #configure power threshold###################################################################################################
    if ts is None:
        ts=[0.3,0.5,0.8]
    
    #configure colormap##########################################################################################################
    cmap_to_use = matplotlib.colormaps[cmap]
    if cmap_to_use.N >100:
        rgba = cmap_to_use(ts)
    else:
        rgba = cmap_to_use(range(len(ts)))
    
    output_hex_colors=[]
    for i in range(len(rgba)):
        output_hex_colors.append(mc.to_hex(rgba[i]))

    if len(highlight)>0:
        sumstats["HUE"] = pd.NA
        sumstats["HUE"] = sumstats["HUE"].astype("Int64")
        sumstats = _process_highlight(
            sumstats          =sumstats, 
            highlight         =highlight, 
            highlight_chrpos  =highlight_chrpos, 
            highlight_windowkb=highlight_windowkb, 
            snpid             =snpid, 
            chrom             =chrom, 
            pos               =pos
        )
    ##################################################################################################
    fig, ax = plt.subplots(**figargs)
    
    ##creating power line############################################################################################
    if mode=="q":
        for i,t in enumerate(ts):
            xpower = get_beta(
                mode      ="q",          
                eaf_range =maf_range,
                beta_range=beta_range, 
                n         =n,
                t         =t,
                sig_level =sig_level,
                n_matrix  =n_matrix
            )
            
            xpower2   = xpower.copy()
            xpower2[1]= -xpower2[1] 
            xpower2[1]= xpower2[1] * yscale_factor
            xpower[1] = xpower[1] * yscale_factor

            lines = LineCollection([xpower2,xpower], label=t,color=output_hex_colors[i],zorder=0)

            ax.add_collection(lines)
    else:
        for i,t in enumerate(ts):
            xpower = get_beta_binary(
                eaf_range =maf_range,
                beta_range=beta_range, 
                prevalence=prevalence,
                or_to_rr  =False,
                ncase     =ncase, 
                ncontrol  =ncontrol, 
                t         =t,
                sig_level =sig_level,
                n_matrix  =n_matrix
            )

            xpower2   = xpower.copy()
            xpower2[1]= -xpower2[1] 
            xpower2[1]= xpower2[1] * yscale_factor
            xpower[1] = xpower[1] * yscale_factor
            lines     = LineCollection([xpower2,xpower], label=t,color=output_hex_colors[i])

            ax.add_collection(lines)


    ###################################################################################################
    # get abs  and convert using scaling factor
    sumstats[beta] = sumstats[beta]*yscale_factor
    sumstats["ABS_BETA"] = sumstats[beta].abs()

    ##################################################################################################
    size_norm = (sumstats["ABS_BETA"].min(), sumstats["ABS_BETA"].max())
    ## if highlight  ##################################################################################################
    dots = sns.scatterplot(data=sumstats,
                    x=maf,
                    y=beta,
                    size=size, 
                    ax=ax, 
                    sizes=sizes,
                    size_norm=size_norm,
                    legend=True, 
                    edgecolor="black",
                    alpha=0.8,
                    zorder=2,
                    **scatter_args)
    
    if len(highlight) >0:
        
        legend   = None
        style    =None
        linewidth=0
        edgecolor="black"

        if pd.api.types.is_list_like(highlight[0]) and highlight_chrpos==False:
            for i, highlight_set in enumerate(highlight):
                scatter_args["color"]=highlight_color[i%len(highlight_color)]
                log.write(" -Highlighting set {} target loci...".format(i+1),verbose=verbose)
                sns.scatterplot(
                    data     =sumstats.loc[sumstats["HUE"]==i], 
                    x        =maf,
                    y        =beta,
                    legend   =legend,
                    style    =style,
                    size     =size,
                    sizes    =sizes,
                    size_norm=size_norm,
                    linewidth=linewidth,
                    zorder   =3+i,
                    ax       =ax,
                    edgecolor=edgecolor,
                    **scatter_args
                )  
        else:
            log.write(" -Highlighting target loci...",verbose=verbose)
            scatter_args["color"]=highlight_color
            sns.scatterplot(
                data     =sumstats.loc[sumstats["HUE"]==0], 
                x        =maf,
                y        =beta,
                legend   =legend,
                size     =size,
                sizes    =sizes,
                size_norm=size_norm,
                zorder   =3,
                ax       =ax,
                edgecolor="black",
                **scatter_args
            )  
    
    ####################################################################################################################
    
    h,l = ax.get_legend_handles_labels()
    if len(ts)>0:
        l1 = ax.legend(h[:int(len(ts))],l[:int(len(ts))], title="Power", loc="upper right",fontsize =fontsize,title_fontsize=fontsize)
        for line in l1.get_lines():
            line.set_linewidth(5.0)
    if hue is None:
        l2 = ax.legend(h[int(len(ts)):],l[int(len(ts)):], title=size, loc="lower right",fontsize =fontsize,title_fontsize=fontsize)
    else:
        l2 = ax.legend(h[int(len(ts)):],l[int(len(ts)):], title=None, loc="lower right",fontsize =fontsize,title_fontsize=fontsize)
    if len(ts)>0:
        ax.add_artist(l1)

    ##################################################################################################

    ax.tick_params(axis='y', labelsize=fontsize)

    ax.axhline(y=0,color="grey",linestyle="dashed")
    
    if xscale== "log":
        ax.set_xscale('log')
        rotation=0
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(min(sumstats[maf].min()/2,0.001/2),0.52)
    else:
        rotation=90    
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(-0.02,0.52)
        
    if xlim is not None:
        ax.set_xlim(xlim)
    
    if ylim is not None:
        ax.set_ylim(ylim)

    if yticks is not None:
        ax.set_yticks(yticks, yticklabels)

    ax.set_ylabel(ylabel,fontsize=fontsize)
    ax.set_xlabel(xlabel,fontsize=fontsize)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)

    ############  Annotation ##################################################################################################
    if len(anno_set)>0:
        variants_toanno = sumstats[sumstats[snpid].isin(anno_set)].reset_index(drop=True)
        if (variants_toanno.empty is not True):
            variants_toanno = annogene(
                variants_toanno,
                id     ="rsID",
                chrom  =chrom,
                pos    =pos,
                log    =log,
                build  =build,
                source =anno_source,
                verbose=verbose,
                gtf_path='/mnt/0A2AAC152AABFBB7/CGE/luxgiant-dstream/GCF_000001405.40_GRCh38.p14_genomic.gtf'
            ).rename(columns={"GENE":"GENENAME"})
            
        texts_u=[]
        texts_d=[]

        if len(variants_toanno)>0:
            maxy = variants_toanno[beta].abs().max()
            variants_toanno["ADJUSTED_i"] = np.nan 
            y_span = 0.5
                
            if sort == "beta" : 
                variants_toanno = variants_toanno.sort_values(by=beta, key=np.abs, ascending= False)
            else:
                variants_toanno = variants_toanno.sort_values(by=maf, key=np.abs, ascending= True)
                
            if anno_style == "expand":

                min_factor=None
                    
                if len(variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"])>1:
                    variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"] =\
                        adjust_text_position(
                            variants_toanno.loc[variants_toanno[beta]>0,maf].values.copy(), 
                            y_span, 
                            repel_force=repel_force,
                            max_iter   =anno_max_iter,
                            log        =log,
                            amode      =xscale,
                            verbose    =verbose,
                            min_factor =min_factor
                        )

                if len(variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"])>1:
                    variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"] = \
                        adjust_text_position(
                            variants_toanno.loc[variants_toanno[beta]<0,maf].values.copy(), 
                            y_span, 
                            repel_force=repel_force,
                            max_iter   =anno_max_iter,
                            log        =log,
                            amode      =xscale,
                            verbose    =verbose,
                            min_factor =min_factor
                        )

            variants_beta_split = [
                variants_toanno.loc[variants_toanno[beta]<0,:], 
                variants_toanno.loc[variants_toanno[beta]>0,:]
            ]
            for variants_toanno_half in variants_beta_split:
                if len(variants_toanno_half)<1:
                    continue
                last_pos = min(variants_toanno_half[maf])/2
                for index, row in variants_toanno_half.iterrows():
                    
                    armB_length_in_point = ax.transData.transform((0,1.1*maxy))[1]-ax.transData.transform((0, abs(row[beta])))[1]
                    armB_length_in_point = armB_length_in_point*arm_scale

                    if anno_style == "right" :
                        #right style
                        if row[maf]>last_pos*(repel_force+1):
                            last_pos=row[maf]
                        else:
                            last_pos*= (repel_force+1)
                    elif anno_style == "expand" :
                        last_pos = row["ADJUSTED_i"]

                    if anno_style == "right"  or anno_style == "expand":
                        if row[beta] >0 :
                            texts_u.append(
                                ax.annotate(
                                    row[anno], 
                                    xy=(row[maf], row[beta]),
                                    xytext=(last_pos , 1.2*maxy),
                                    arrowprops=dict(
                                        relpos=(0,0),
                                        arrowstyle="-|>",
                                        linestyle='--', 
                                        facecolor='black',
                                        connectionstyle="arc,angleA=-90,armA={},angleB=0,armB=0,rad=0".format(armB_length_in_point)
                                    ),
                                    rotation=90,
                                    ha="left",
                                    va="bottom",
                                    **anno_args
                                )
                            )
                        else:
                            texts_d.append(
                                ax.annotate(
                                    row[anno], 
                                    xy=(row[maf], row[beta]),
                                    xytext=(last_pos , -1.2*maxy),
                                    arrowprops=dict(
                                        relpos=(0,1),
                                        arrowstyle="-|>", 
                                        linestyle='--',
                                        facecolor='black',
                                        connectionstyle="arc,angleA=90,armA={},angleB=0,armB=0,rad=0".format(armB_length_in_point)
                                    ),
                                    rotation=90,
                                    ha="left",
                                    va="top",
                                    **anno_args
                                )
                            )
                        
                    if anno_style=="tight":
                        texts_d.append(ax.text(row[maf], row[beta], row[anno]))

            if anno_style=="tight":
                adjust_text(
                    texts_d, 
                    autoalign    =True,
                    precision    =0.001,
                    lim          =1000, 
                    expand_text  =(0.5,0.5), 
                    expand_points=(0.5,0.5),
                    force_objects=(0.5,0.5), 
                    ax           =ax
                )
        
    ############  Annotation ##################################################################################################
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'trumpet_plot.jpeg'), dpi=500)
    plt.show()

    return None
    
def new_trumpet(df_gwas:pd.DataFrame, df_freq:pd.DataFrame, plot_dir:pd.DataFrame, snp_col:str, chr_col:str, pos_col:str, maf_col:str, beta_col:str, power_ts:list, n_case:int, n_control:int, sample_size:int=None, n_col:str='', p_col:str=None, prevalence:int=None, mode:str='binary', p_filter:float=5e-8, to_highlight:list=[], to_annotate:list=[], cmap:str="cool", power_sig_level:float=5e-8, build='38', gtf_path:str=None, save_name:str='trumpet_plot.jpeg')->bool:

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

    else:
        gwas_df = df_gwas.copy()
        del df_gwas
            
    if maf_col not in df_gwas.columns:
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
                n         =n,
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
    ax.axhline(y=0,color="black",linestyle="solid")

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
