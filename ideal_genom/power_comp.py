# This functions are used to calculate the power of a genetic association study for a given effect size and sample size.
# They were fetched from GWASlab https://github.com/Cloufield/gwaslab and rewritten for better clarity and readability.
# The original code was written by Cloufield and is licensed under the MIT License.

import logging

import pandas as pd
import numpy as np
import scipy.stats as ss

from typing import Tuple, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def calculate_power_quantitative(beta: np.ndarray, eaf: np.ndarray, sample_size: int, sig_level: float = 5e-8, variance: float = 1) -> np.ndarray:
    """
    Calculate statistical power for detecting genetic associations in quantitative traits.
    This function computes the statistical power for quantitative trait association tests
    based on effect sizes (beta), effect allele frequencies (EAF), sample size, 
    significance level, and trait variance.
    
    Parameters
    ----------
    beta : numpy.ndarray
        Array of effect sizes (regression coefficients)
    eaf : numpy.ndarray
        Array of effect allele frequencies (between 0 and 1)
    sample_size : int
        Number of individuals in the study
    sig_level : float, optional
        Significance threshold (default: 5e-8, typical for GWAS)
    variance : float, optional
        Phenotypic variance of the trait (default: 1)
    
    Returns
    -------
    numpy.ndarray
        Array of statistical power values corresponding to each variant
    
    Notes
    -----
    The calculation uses the non-central chi-square distribution to determine 
    the probability of detecting true associations at the specified significance level.
    Power is calculated as 1 minus the cumulative distribution function of the 
    non-central chi-square distribution with 1 degree of freedom.
    """
                
    c = ss.chi2.isf(sig_level,df=1) # critical value for chi-square test

    h2 = 2*eaf*(1-eaf)*(beta**2) # heritability contribution
                
    ncp = sample_size * h2/variance # non-centrality parameter 
                
    power = 1 - ss.ncx2.cdf(c,df=1,nc=ncp) # statistical power
                
    return power

def calculate_power_binary(beta: np.ndarray, daf: np.ndarray, prevalence: float, ncase: int, ncontrol: int, sig_level: float = 5e-8, or_to_rr: bool = False) -> np.ndarray:
    """
    Calculate the statistical power of genetic association tests in a case-control study with binary traits.
    This function computes the power to detect genetic associations in a case-control study design,
    assuming an additive genetic model. It accounts for allele frequencies, effect sizes,
    disease prevalence, and sample sizes to determine the probability of detecting true associations
    at a specified significance level.
    
    Parameters
    ----------
    beta : np.ndarray
        Log odds ratios representing effect sizes of variants.
    daf : np.ndarray
        Disease allele frequencies.
    prevalence : float
        Disease prevalence in the population.
    ncase : int
        Number of cases in the study.
    ncontrol : int
        Number of controls in the study.
    sig_level : float, optional
        Significance level threshold for declaring associations (default: 5e-8, standard for GWAS).
    or_to_rr : bool, optional
        If True, converts odds ratios to risk ratios using the disease prevalence (default: False).
    
    Returns
    -------
    np.ndarray
        Array of statistical power values for each variant, representing the probability
        of detecting a true association at the specified significance level.
    
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

def get_beta_binary(prevalence: float, ncase: int, ncontrol: int, eaf_range: tuple = (0.0001,0.5), beta_range: tuple = (0.0001,5), t: float = 0, sig_level: float = 5e-8, n_matrix: int = 500, or_to_rr: bool = False):
    """
    Find combinations of effect allele frequencies (EAF) and effect sizes (beta) that achieve 
    a specified statistical power for binary traits in genetic association studies.
    
    Parameters
    ----------
    prevalence : float
        Disease prevalence in the population.
    ncase : int
        Number of cases in the study.
    ncontrol : int
        Number of controls in the study.
    eaf_range : tuple, optional
        Range of effect allele frequencies to consider (min, max). Default is (0.0001, 0.5).
    beta_range : tuple, optional
        Range of effect sizes to consider (min, max). Default is (0.0001, 5).
    t : float, optional
        Target power threshold. Combinations with power >= t will be returned. Default is 0.
    sig_level : float, optional
        Significance level (alpha) for power calculation. Default is 5e-8 (genome-wide significance).
    n_matrix : int, optional
        Resolution of the EAF-beta matrix. Default is 500.
    or_to_rr : bool, optional
        If True, converts odds ratios to relative risks using the provided prevalence.
        If False, approximates genetic relative risk (GRR) using odds ratios. Default is False.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing EAF and beta combinations that achieve the target power.
        The DataFrame has columns 'eaf' and 'beta' and is sorted by 'eaf' in descending order.
    
    Notes
    -----
    When prevalence is less than 10%, GRR is very similar to OR if or_to_rr is set to False.
    The function uses a greedy algorithm to find EAF-beta combinations with power >= t.
    """
    
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