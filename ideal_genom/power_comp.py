# from GWASlab

import pandas as pd
import numpy as np
import scipy.stats as ss

def calculate_power_quantitative(beta: float, eaf: float, sample_size: int, sig_level: float = 5e-8, variance: float =1) -> float:

    """
    Calculate the statistical power for a single beta value for quantitative traits.

    Parameters:
    -----------
    beta (float): 
        The effect size.
    eaf (float): 
        The effect allele frequency.
    n (int): 
        The sample size.
    sig_level (float): 
        The significance level (default is 5e-8).
    variance (float): 
        The variance (default is 1).

    Returns:
    --------
    float: The calculated power.
    """
                
    c = ss.chi2.isf(sig_level,df=1) # critical value for chi-square test

    h2 = 2*eaf*(1-eaf)*(beta**2) # heritability contribution
                
    ncp = sample_size * h2/variance # non-centrality parameter 
                
    power = 1 - ss.ncx2.cdf(c,df=1,nc=ncp) # statistical power
                
    return power

def get_beta_quantitative(
    eaf_range: tuple = (0.0001, 0.5),
    beta_range: tuple = (0.0001, 10),
    t: float = 0,
    n: int = None,
    sig_level: float = 5e-8,
    variance: float = 1,
    n_matrix: int = 500
):
    if t <= 0 or t > 1:
        return pd.DataFrame(columns=["eaf", "beta"])  # Invalid threshold or no computation needed

    if n is None:
        raise ValueError("Sample size 'n' must be specified.")

    # Generate grid of eaf and beta values
    eafs = np.linspace(eaf_range[0], eaf_range[1], n_matrix)
    betas = np.linspace(beta_range[0], beta_range[1], n_matrix)
    eaf_grid, beta_grid = np.meshgrid(eafs, betas, indexing="ij")

    # Compute power for each (eaf, beta) pair
    power_matrix = np.vectorize(calculate_power_quantitative)(
        beta=beta_grid, eaf=eaf_grid, n=n, sig_level=sig_level, variance=variance
    )

    # Extract (eaf, beta) values where power >= threshold `t`
    mask = power_matrix >= t
    eaf_beta = np.column_stack((eaf_grid[mask], beta_grid[mask]))

    return pd.DataFrame(eaf_beta, columns=["eaf", "beta"])

def get_beta_binary(
              prevalence=None,
              ncase=None,
              ncontrol=None,
              eaf_range=(0.0001,0.5),
              beta_range=(0.0001,10),
              t=0,
              sig_level= 5e-8,
              vary=1,
              n_matrix=500,
              or_to_rr=False
             ):
    if t >0:
        def calculate_power_single(
                            beta, 
                            daf, 
                            prevalence,
                            ncase, 
                            ncontrol, 
                            sig_level=5e-8,
                            or_to_rr=False):
                
                aaf = daf**2
                abf = 2 * (daf) * (1 - daf)
                bbf = (1- daf)**2

                if or_to_rr == False:
                    genotype_or = np.exp(beta)
                    genotype_rr = genotype_or
                else:
                    genotype_or = np.exp(beta)
                    genotype_rr = genotype_or/ ((1-prevalence)+(genotype_or*prevalence))
                    # https://jamanetwork.com/journals/jama/fullarticle/188182
                # additive
                x = [ 2*genotype_rr-1, genotype_rr, 1 ] 
                aap= x[0] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
                abp= x[1] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
                bbp= x[2] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
                pcase= (aap * aaf + abp * abf*0.5) / prevalence
                pcontrol=((1-aap )* aaf + (1-abp )* abf*0.5) / (1 - prevalence)
                vcase = pcase *(1-pcase)
                vcontrol =pcontrol *(1-pcontrol)
                num= (pcase - pcontrol)
                den= np.sqrt( (vcase/ncase +  vcontrol/ncontrol)*0.5 )
                u = num / den
                c = ss.norm.isf(sig_level/2)
                power = 1 - ss.norm.cdf(c-u) + ss.norm.cdf(-c-u)
                return power
        
        eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
        betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
        
        print(" -Updating eaf-beta matrix...")
        if or_to_rr ==False:
            print(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....")
        else:
            print(" -OR is converted to GRR using base prevalence: {}".format(prevalence))
        
        eaf_beta_matrix = np.array([
            [calculate_power_single(beta=betas[j], daf=eafs[i], ncase=ncase, ncontrol=ncontrol, prevalence=prevalence, sig_level=sig_level, or_to_rr=or_to_rr) for j in range(n_matrix)]
            for i in range(n_matrix)
        ])
        
        print(" -Extracting eaf-beta combinations with power = {}...".format(t))
        i,j=1,1
        eaf_beta = []
        while i<n_matrix-1 and j<n_matrix-1:
            if eaf_beta_matrix[i,j] < t:
                j+=1
            else:
                i+=1
                eaf_beta.append((eafs[i],betas[j]))
                
    return pd.DataFrame(eaf_beta)