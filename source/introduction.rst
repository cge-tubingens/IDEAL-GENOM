Introduction
============

IDEAL-GENOM is a a dedicated Python library to run Genome Wide Association Analysis (GWAS) on large-scale genomic data and to ease downstream analysis. It is desgined to be easy to use and to be integrated in a pipeline. This library encompases several years of experience in the field of GWAS and is designed to be used by both beginners and experts in the field. It is designed to be used in a Jupyter notebook or in a Python script.

It is intended to receive the output data of an imputation process and streamline the GWAS analysis. We have incorporated two divergent branches of GWAS, one for fixed effects and another for random effects. The fixed effects branch is based on the `PLINK <https://www.cog-genomics.org/plink2>`_ software and the random effects branch is based on the `GCTA <https://cnsgenomics.com/software/gcta/>`_ software.

The library also contains a set of tools to facilitate the analysis of the results of the GWAS, such as Manhattan plots, QQ plots, Miami plots, effect size comparison plot and trumpet plots. 

Though there are many resources already developed to tackle these problems, we belive that IDEAL-GENOM is a valuable tool for the scientific community, because under the same umbrellam several tasks can be performed and also allowing to use its features independently. 