# from GWASlab
import logging

import pandas as pd
import numpy as np
import scipy.stats as ss

from typing import Tuple, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def calculate_power_quantitative(beta: np.ndarray, eaf: np.ndarray, sample_size: int, sig_level: float = 5e-8, variance: float = 1) -> float:

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

def calculate_power_binary(beta, daf, prevalence: float, ncase: int, ncontrol: int, sig_level: float = 5e-8, or_to_rr: bool = False):
    
    """
    Calculate statistical power for genetic association testing under binary trait.
    This function calculates the statistical power for genetic association studies with binary traits
    (e.g., case-control studies) using an additive genetic model.

    Parameters
    ----------
    beta : float
        Effect size (log odds ratio)
    daf : float
        Disease allele frequency (between 0 and 1)
    prevalence : float
        Disease prevalence in the population (between 0 and 1)
    ncase : int
        Number of cases
    ncontrol : int
        Number of controls
    sig_level : float, optional
        Significance level (alpha), default is 5e-8 (genome-wide significance)
    or_to_rr : bool, optional
        If True, converts odds ratio to relative risk using prevalence, default is False

    Returns
    -------
    float
        Statistical power (between 0 and 1)
    
    Notes
    -----
    The function implements power calculation for an additive genetic model using normal
    approximation. When or_to_rr is True, it converts odds ratio to relative risk using
    the formula from Zhang and Yu (JAMA, 1998).
    References
    ----------
    Zhang J, Yu KF. What's the Relative Risk? JAMA. 1998;280(19):1690-1691.
    doi:10.1001/jama.280.19.1690
    """
                
    aaf = daf**2
    abf = 2 * (daf) * (1 - daf)
    bbf = (1- daf)**2

    if not or_to_rr:
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

    pcase= (aap * aaf + abp * abf*0.5) / prevalence
    pcontrol=((1-aap )* aaf + (1-abp )* abf*0.5) / (1 - prevalence)

    vcase = pcase *(1-pcase)
    vcontrol =pcontrol *(1-pcontrol)

    num= (pcase - pcontrol)
    den= np.sqrt( (vcase/ncase +  vcontrol/ncontrol)*0.5 )
    u = num / den

    # Ensure the significance level is divided by 2 for a two-tailed test
    c = ss.norm.isf(sig_level / 2)

    power = 1 - ss.norm.cdf(c-u) + ss.norm.cdf(-c-u)

    return power

def get_beta_quantitative(
    eaf_range: Tuple[float, float] = (0.00001, 0.5),
    beta_range: Tuple[float, float] = (0.00001, 5),
    t: float = 0,
    sample_size: Optional[int] = None,
    sig_level: float = 5e-8,
    variance: float = 1,
    n_matrix: int = 500
) -> pd.DataFrame:
    """
    Calculate beta values for quantitative traits based on power threshold.

    Parameters:
    -----------
    eaf_range : Tuple[float, float]
        Range of effect allele frequencies (min, max)
    beta_range : Tuple[float, float]
        Range of beta values (min, max)
    t : float
        Power threshold (0 to 1)
    sample_size : int
        Sample size
    sig_level : float
        Significance level
    variance : float
        Variance
    n_matrix : int
        Size of the grid for calculations

    Returns:
    --------
    pd.DataFrame
        DataFrame containing eaf and beta values meeting the threshold
    """
    # Input validation
    if not 0 < t <= 1:
        return pd.DataFrame(columns=["eaf", "beta"])
    
    if sample_size is None or sample_size <= 0:
        raise ValueError("Sample size must be a positive integer")
        
    if not all(0 <= x <= 1 for x in eaf_range):
        raise ValueError("EAF range values must be between 0 and 1")
        
    if any(x < 0 for x in beta_range):
        raise ValueError("Beta range values must be non-negative")
        
    if eaf_range[0] > eaf_range[1] or beta_range[0] > beta_range[1]:
        raise ValueError("Range values must be in ascending order (min, max)")

    # Generate grid of eaf and beta values
    eafs = np.linspace(eaf_range[1], eaf_range[0], n_matrix)
    betas = np.linspace(beta_range[0], beta_range[1], n_matrix)
    
    logger.info(f" -Calculating power matrix with parameters: ")
    logger.info(f"  --EAF range: {eaf_range}")
    logger.info(f"  --Beta range: {beta_range}")
    logger.info(f"  --Sample size: {sample_size}")
    logger.info(f"  --Significance level: {sig_level}")
    power_matrix = calculate_power_quantitative(
        beta=betas[np.newaxis, :],  # Make row vector
        eaf=eafs[:, np.newaxis],    # Make column vector
        sample_size=sample_size,
        sig_level=sig_level,
        variance=variance
    )

    # Find threshold boundary more efficiently
    eaf_beta = []
    i, j = 1, 1  # Start from (1,1) to avoid edge effects
    
    while i < n_matrix - 1 and j < n_matrix - 1:
        if power_matrix[i, j] < t:
            j += 1
        else:
            eaf_beta.append({
                'eaf': eafs[i],
                'beta': betas[j],
                'power': power_matrix[i, j]
            })
            i += 1

    # Create DataFrame with results
    result_df = pd.DataFrame(eaf_beta)
    logger.info(f" -Found {result_df.shape[0]} eaf-beta combinations with power >= {t} and columns: {result_df.columns}")
    
    # Add column names and sort if there are results
    if not result_df.empty:
        result_df = result_df.sort_values('eaf', ascending=False)
        
    return result_df[['eaf', 'beta']]  # Return only eaf and beta columns

def get_beta_binary(prevalence: float = None, ncase: int = None, ncontrol: int = None, eaf_range: tuple = (0.0001,0.5), beta_range: tuple = (0.0001,5), t: float = 0, sig_level: float = 5e-8, n_matrix: int = 500, or_to_rr: bool = False):
    
    eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
    betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
        
    print(" -Updating eaf-beta matrix...")
    if or_to_rr ==False:
        logger.info(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....")
    else:
        logger.info(" -OR is converted to GRR using base prevalence: {}".format(prevalence))
        
    power_matrix = calculate_power_binary(
        beta=betas[np.newaxis, :],  # Make row vector
        daf=eafs[:, np.newaxis],    # Make column vector
        ncase=ncase,
        ncontrol=ncontrol,
        prevalence=prevalence,
        sig_level=sig_level,
        or_to_rr=or_to_rr
    )
        
    logger.info(" -Extracting eaf-beta combinations with power = {}...".format(t))
    eaf_beta = []
    i, j = 1, 1  # Start from (1,1) to avoid edge effects
    
    while i < n_matrix - 1 and j < n_matrix - 1:
        if power_matrix[i, j] < t:
            j += 1
        else:
            eaf_beta.append({
                'eaf': eafs[i],
                'beta': betas[j],
                'power': power_matrix[i, j]
            })
            i += 1

     # Create DataFrame with results
    result_df = pd.DataFrame(eaf_beta)
    logger.info(f" -Found {result_df.shape[0]} eaf-beta combinations with power >= {t} and columns: {result_df.columns}")
                
    if not result_df.empty:
        result_df = result_df.sort_values('eaf', ascending=False)
        
    return result_df[['eaf', 'beta']]  # Return only eaf and beta columns