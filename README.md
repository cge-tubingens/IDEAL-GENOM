# Genome Wide Association Analysis (GWAS) Pipeline

This Python package is designed to execute a GWAS pipeline, encapsulating several years of reseach as well as new experiences at CGE Tübingen.

# Basic Requirements

The GWAS pipeline rest on three tools: `PLINK 1.9`, `PLINK 2.0` and `GCTA`. The `luxgiant-dstream` serves as a wrapper fopr the various pipeline steps. To tun it the above programs must be installed on the system.

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
    "maf" : 0.05,
    "geno": 0.1,
    "mind": 0.1,
    "hwe" : 0.00000005,
    "indep-pairwise": [50, 5, 0.2],
    "pca": 10,
    "ci": 0.95,
    "annotate": true
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
    "prep_ds"   : true,
    "gwas_fixed" : true,
    "gwas_random": true
}
```

With the above configuration, all three steps will run seamlessly, which is the recommended initial configuration. If you want to skip some steps, change the value to `false`. For example,

```
{
    "preps_ds"   : false,
    "gwas_fixed" : true,
    "gwas_random": false
}
```

allows you to run only the GWAS with a fixed effect model. Note that an exception will be raised if the preparation step has not been run, as the necessary files for the fixed model to run would not be available.

## Dependable Files

This folder should contain additional files to run the quality control pipeline. The structure inside the directory should be as follows:

```
dependables
    |
    |---high-LD-regions.txt
```

## Usage

The pipeline is easy to use. Once installed in the system or in a virtual enviroment one needs to run the following command:

```
python3 luxgiant-dstream --path_params <path to parameters.JSON> 
                         --file_folders <path to paths.JSON> 
                         --steps <path to steps.JSON>
```

The first three parameters are the path to the three configuration files.
