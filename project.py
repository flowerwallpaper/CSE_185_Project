import argparse
import pandas as pd
import numpy as np
import pysam
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from qqman import qqman

def parse_args():
    parser = argparse.ArgumentParser(description='Run PLINK analysis.')
    parser.add_argument('--pheno', required=True, help='Path to the phenotype file')
    parser.add_argument('--geno', required=True, help='Path to the genotype file (VCF format)')
    parser.add_argument('--maf', type=float, default=0.05, help='Minor allele frequency threshold (default: 0.05)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load phenotypes
    phenotypes = pd.read_csv(args.pheno, sep='\t', header=None)
    phenotypes.columns = ['Sample', 'ID', 'Phenotype']

    # Read the VCF file
    vcf_file = pysam.VariantFile(args.geno)
    genotype_data = []
    variant_ids = []
    chromosomes = []
    positions = []

    for idx, record in enumerate(vcf_file):
        variant_id = record.id
        chromosomes.append(record.chrom)
        positions.append(record.pos)
        variant_ids.append(variant_id)
        
        genotypes = []
        for sample in record.samples.values():
            genotype = sample['GT']
            if None in genotype:
                genotypes.append(np.nan)
            else:
                allele_count = sum(genotype)
                genotypes.append(allele_count)
        genotype_data.append(genotypes)

    print("Genotype file parsed")

    # Convert genotype data to a NumPy array
    genotype_array = np.array(genotype_data).T  # Transpose to have samples as rows and variants as columns

    print("Genotype is numpy array")

    # Calculate Minor Allele Frequency (MAF)
    maf = np.nanmean(genotype_array / 2, axis=0)  # Assuming diploid organisms, hence divide by 2

    print("MAF calculated")

    # Filter SNPs with MAF > threshold
    maf_threshold = args.maf
    indices = np.where(maf > maf_threshold)[0]
    genotype_array_filtered = genotype_array[:, indices]
    variant_ids_filtered = [variant_ids[i] for i in indices]
    chromosomes_filtered = [chromosomes[i] for i in indices]
    positions_filtered = [positions[i] for i in indices]

    print("Genotype array filtered based on MAF threshold")

    # Handle missing values by imputing with the mean of the column
    nan_mask = np.isnan(genotype_array_filtered)
    means = np.nanmean(genotype_array_filtered, axis=0)
    genotype_array_filtered[nan_mask] = np.take(means, np.where(nan_mask)[1])

    print("Handled missing values")

    # Scale the genotypes to have mean 0 and variance 1
    scaler = StandardScaler()
    scaled_genotype_array = scaler.fit_transform(genotype_array_filtered)

    print("Scaled genotypes")

    # Create the np.array of phenotypes
    pts = phenotypes['Phenotype'].to_numpy()

    print("Shape of pts:", pts.shape)
    print("Shape of X:", genotype_array_filtered.shape)

    # Perform GWAS
    pvals = []
    betas = []
    stats = []
    for i in range(scaled_genotype_array.shape[1]):
        X = sm.add_constant(scaled_genotype_array[:, i])
        model = sm.OLS(pts, X)
        results = model.fit()
        betas.append(results.params[1])
        stats.append(results.tvalues[1])
        pvals.append(results.pvalues[1])

    # Prepare data for output
    output_df = pd.DataFrame({
        'CHR': chromosomes_filtered,
        'SNP': variant_ids_filtered,
        'BP': positions_filtered,
        'BETA': betas,
        'STAT': stats,
        'P': pvals
    })

    # Save to file
    output_df.to_csv('gwas_results.txt', sep='\t', index=False)

    # Manhattan and QQ plot
    data = pd.DataFrame({
        'CHR': chromosomes_filtered,
        'BP': positions_filtered,
        'P': pvals,
        'SNP': variant_ids_filtered
    })

    fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
    fig.set_size_inches((15, 5))

    qqman.manhattan(data, ax=ax0)
    qqman.qqplot(data, ax=ax1)

    plt.show()

if __name__ == "__main__":
    main()
