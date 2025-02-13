# IDEAL-GENOM

This Python package is intended to perform a GWAS pipeline, that starts with the data coming from the imputation process to the files with the summary statistics resulting from a GWAS. We have implemented a 'divergent' pipeline: in one branch we train a fixed effect model and in the other one a random effect model.

# Basic Requirements

The GWAS pipeline rest on three tools: `PLINK 1.9`, `PLINK 2.0`, `GCTA` and `bcftools`. The `IDEAL-GENOM` serves as a wrapper fopr the various pipeline steps. To tun it the above programs must be installed on the system. it is important to remark that `GCTA` is memory intensive, so it is recommended to have a local or virtual machine with a large RAM available.

The pipeline is desgined to seamlessly run with minimal input and produce plots and summary statistics as result. To accomplish this, the following folder structure is expected:

```
projectFolder
    |
    |---inputData
    |
    |---outputData
    |
    |---dependables
    |
    |---configFiles
```

1. The `inputData` folder should contain the files resulting from the imputation process. In the present pipeline we expect that the imputation was done with the Michigan Imputation Server.

2. The `outputData` folder will contain the resultant files of the GWAS pipeline. Below, the pipeline output will be detailed.

3. The `dependables` folder is designed to contain necessary files for the pipeline.

4. The `configFiles` older is essential for the correct functioning of the pipeline. It should contain three configuration files: `parameters.JSON`, `paths.JSON` and `steps.JSON`.

## Configuration Files

These three files contain all the information necessary to run the pipeline.

### GWAS Pipeline parameters

The `parameters.JSON` file contains values for the different parameters used along the pipeline. It expects the following parameters 

```
{
    "post_imputation": {
        "pwd": "<zip_files_password>",
        "r2_thres": 0.3,
        "ref_genome": "<fasta_filename.fa>",
        "annot_vcf": "<annotation_ref.vcf.gz>",
    },
    "preparatory":{
        "maf": 0.01,
        "geno": 0.1,
        "hwe": 5e-6,
        "ind_pair":[50, 5, 0.2],
        "pca": 10 
    },
    "gwas_glm": {
        "maf": 0.01,
        "mind": 0.1,
        "hwe": 5e-8,
        "ci": 0.95,
        "gtf_path": null,
        "build": "38",
        "anno_source": "ensembl"
    },
    "gwas_glmm": {
        "maf": 0.01,
        "gtf_path": null,
        "build": "38",
        "anno_source": "ensembl"
    }
}
```

If you wish to change at least one of the default values, please provide the full information in the configuration file.

### Paths to Project Folders

The `paths.JSON` file contains the addresses to the project folder as well as the prefix of the input and output data. The file must contain the following fields:

```
{
    "input_directory"      : "<path to folder with project input data>",
    "input_prefix"         : "<prefix of the input data>",
    "output_directory"     : "<path to folder where the output data will go>",
    "output_prefix"        : "<prefix for the output data>",
    "dependables_directory": "<path to folder with dependables files>"
}
```

If the CLI is run locally you should provide the full path to the directories.

### Pipeline Steps

The `steps.JSON` file has the following structure:

```
{
    "pos_imputation": true,
    "preparatory"   : true,
    "gwas_glm"      : true,
    "gwas_glmm"     : true
}
```

With the above configuration, all three steps will run seamlessly, which is the recommended initial configuration. If you want to skip some steps, change the value to `false`. For example,

```
{
    "pos_imputation": false,
    "preparatory"   : false,
    "gwas_glm"      : true,
    "gwas_glmm"     : true
}
```

allows you to run only the GWAS with a fixed effect model. Note that an exception will be raised if the preparation step has not been run, as the necessary files for the fixed model to run would not be available.

## Dependable Files

This folder should contain additional files to run the whole pipeline. The structure inside the directory should be as follows:

```
dependables
    |
    |---high-LD-regions.txt
    |---<annotations_file>.vcf.gz
    |---<reference_genome>.fa
    |---<reference_genome>.fa.fai
```

The `high-LD-regions.txt` will be used to prune the raw data in the preparatory steps, while the `<annotations_file>.vcf.gz` (if it is not indexed the pipeline will do it) will be used to get the rsIDs of each SNP and the `<reference_genome>.fa` will be used to normalize the data coming out of imputation.

## Installation and usage

The library can be installed by cloning the GitHub repository:

```
git clone https://github.com/cge-tubingens/IDEAL-GENOM.git
```

### Setting up the environment

The virtual environment can be created using either `Poetry` or `pip`. Since this is a `Poetry`-based project, we recommend using `Poetry`. Once `Poetry` is installed on your system (refer to [`Poetry` documentation](https://python-poetry.org/docs/) for installation details), navigate to the cloned repository folder and run the following command:

```
poetry install
```
It is important to remark that currently the project has been updated to use `Poetry 2.0`.

### Pipeline usage options

#### 1. Inside a virtual environment

After running the `poetry install` activate the virtual environment with 

```
poetry shell
```

 Once the environment is active, you can execute the pipeline with the following command:

```
python3 ideal_genom --path_params <path to parameters.JSON> 
                    --file_folders <path to paths.JSON> 
                    --steps <path to steps.JSON>
```

The first three parameters are the path to the three configuration files.

#### 2. Using `Poetry` directly

One of the benefits of using `Poetry` s that it eliminates the need to activate a virtual environment. Run the pipeline directly with:

```
poetry run python3 ideal_genom --path_params <path to parameters.JSON> 
                               --file_folders <path to paths.JSON> 
                               --steps <path to steps.JSON>
```
#### 3. Jupyter Notebooks

The package includes Jupyter notebooks located in the notebooks folder. Each notebook corresponds to a specific step of the pipeline. Simply provide the required parameters to execute the steps interactively.

Using the notebooks is a great way to gain a deeper understanding of how the pipeline operates.