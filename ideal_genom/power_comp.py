# from GWASlab

import pandas as pd
import numpy as np
import scipy.stats as ss

def get_beta(mode: str = "q", eaf_range: tuple = (0.0001,0.5), beta_range: tuple =(0.0001,10), t=0, n=None, sig_level: float= 5e-8, vary=1, n_matrix: int = 500):

    if mode=="q":
        if t >0:
            def calculate_power_single(
                                beta, 
                                eaf, 
                                n, 
                                sig_level=5e-8,
                                vary=1):
                
                c = ss.chi2.isf(sig_level,df=1)
                h2 = 2*eaf*(1-eaf)*(beta**2)
                NCP = n * h2/vary
                power = 1 - ss.ncx2.cdf(c,df=1,nc=NCP)
                return power
            
            eaf_beta_matrix = np.zeros((n_matrix,n_matrix),dtype=float)
            eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
            betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
            
            for i in range(n_matrix):
                    eaf_beta_matrix[i,] = calculate_power_single(beta=betas,eaf=eafs[i],n=n,sig_level=sig_level,vary=vary)
            
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
        
        eaf_beta_matrix = np.zeros((n_matrix,n_matrix),dtype=float)
        eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
        betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
        
        print(" -Updating eaf-beta matrix...")
        if or_to_rr ==False:
            print(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....")
        else:
            print(" -OR is converted to GRR using base prevalence: {}".format(prevalence))
        
        for i in range(n_matrix):
                eaf_beta_matrix[i,] = calculate_power_single(beta=betas,
                                                                daf=eafs[i],
                                                                ncase=ncase,
                                                                ncontrol=ncontrol,
                                                                prevalence=prevalence,
                                                                sig_level=sig_level,
                                                                or_to_rr=or_to_rr)
        
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