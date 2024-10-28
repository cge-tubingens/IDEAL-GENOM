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
import textalloc as ta

from matplotlib.axes import Axes
from pandas.core.groupby.generic import DataFrameGroupBy

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import scipy as sp
from matplotlib.collections import LineCollection
import matplotlib.colors as mc
import matplotlib
from adjustText import adjust_text
from gwaslab.util_in_get_sig import annogene
from gwaslab.viz_aux_annotate_plot import annotate_single
from gwaslab.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz_aux_save_figure import save_figure
from gwaslab.viz_plot_mqqplot import _process_highlight
from gwaslab.g_Log import Log
from gwaslab.util_in_calculate_power import get_beta
from gwaslab.util_in_calculate_power import get_power
from gwaslab.util_in_calculate_power import get_beta_binary
from gwaslab.util_in_fill_data import filldata



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






    

def draw_trumpet_plot(df_gwas:pd.DataFrame,
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
    plt.savefig(os.path.join(plot_dir, 'trumpet_plot.jpeg'))
    plt.show()

    log.write("Finished creating trumpet plot!", verbose=verbose)
    return None
    
