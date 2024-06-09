# CLINK #

## Overview ## 
This script performs a genome-wide association study (GWAS) using phenotype and genotype data. It reads phenotype and genotype files, filters SNPs based on minor allele frequency (MAF), handles missing values, scales the data, performs linear regression to find associations, and outputs the results along with Manhattan and QQ plots.

## Requirements ##
To run this script, you need the following Python packages:
- argparse
- pandas
- numpy
- pysam
- matplotlib
- sklearn
- statsmodels
- qqman
You can install them with pip.

## Usage ##
To run the script, use the command line to provide the necessary input files and parameters:
python script_name.py --pheno path_to_pheno_file --geno path_to_geno_file --covar path_to_covar_file --maf minor_allele_frequency_threshold
### Arguments ###
--pheno: Path to the phenotype file (required)
--geno: Path to the genotype file in VCF format (required)
--covar: Path to the covariates file (optional)
--maf: Minor allele frequency threshold (default: 0.05)

## Script Workflow ##
Parsing Arguments: The script uses argparse to handle command-line arguments.
Loading Phenotypes: Reads the phenotype data from a file into a pandas DataFrame.
Reading Genotype Data: Reads the genotype data from a VCF file using pysam.
Calculating MAF: Calculates the minor allele frequency for each SNP.
Filtering SNPs: Filters SNPs based on the provided MAF threshold.
Handling Missing Values: Imputes missing genotype values with the mean of each SNP.
Scaling Genotypes: Scales the genotype data to have a mean of 0 and variance of 1 using StandardScaler.
Performing GWAS: Performs linear regression for each SNP using statsmodels to identify associations with the phenotype.
Outputting Results: Saves the GWAS results to a file and generates Manhattan and QQ plots.

## Output ## 
GWAS Results: A tab-delimited file gwas_results.txt containing the following columns:
CHR: Chromosome number
SNP: SNP identifier
BP: Base pair position
BETA: Regression coefficient
STAT: Test statistic
P: p-value
Plots: Manhattan and QQ plots are displayed to visualize the GWAS results.

## Example Files ##
This contains data from two sources--lab 3 and the Personal Genome Project. lab3_gaws.phen contains the phenotypic data, lab3_gwas.vcf.gz contains the genotypic data, lab3_gwas.vcf.gz.tbi is the index file, and lab3_gwas.assoc.linear is the covariate file for the dataset. Note that the index file is not taken as an argument, but it still needs to be downloaded for the program to run. 
The PGP data is as follows: BMI.csv contains phenotypic data, BMI_merged.vcf.gz contains genotypic data, and BMI_merged.vcf.gz.tbi is the indexed file. 

