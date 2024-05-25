#pip install pandas vcfpy statsmodels 
import pandas as pd
import vcfpy
from statsmodels.formula.api import ols
import numpy as np
import matplotlib.pyplot as plt
import gzip

# # Get phenotypes
# phenotypes = pd.read_csv('lab3_gwas.phen')

# # Decompress and read the VCF file
# with gzip.open('lab3_gwas.vcf.gz', 'rt') as f:
#     vcf_reader = vcfpy.Reader(f)
#     samples = vcf_reader.header.samples.names

#     genotype_data = []
#     variant_ids = []
#     for record in vcf_reader:
#         variant_ids.append(f"{record.CHROM}_{record.POS}")
#         genotypes = []
#         for call in record.calls:
#             if call.data['GT'] == '0/0':
#                 genotypes.append(0)
#             elif call.data['GT'] == '0/1':
#                 genotypes.append(1)
#             elif call.data['GT'] == '1/1':
#                 genotypes.append(2)
#             else:
#                 genotypes.append(None)
#         genotype_data.append(genotypes)


# Sample phenotypes data
phenotypes_data = {
    'Sample': ['Sample1', 'Sample2', 'Sample3', 'Sample4'],
    'Phenotype': [0.75, 0.60, 0.90, 0.80]
}

# Sample genotype data
genotype_data = {
    'Sample': ['Sample1', 'Sample2', 'Sample3', 'Sample4'],
    'Variant1': [0, 1, 0, 2],
    'Variant2': [1, 2, 0, 1],
    'Variant3': [2, 0, 1, 1]
}
phenotypes = pd.DataFrame(phenotypes_data)
genotype_df = pd.DataFrame(genotype_data)
variant_ids = ['Variant1', 'Variant2', 'Variant3']

# #make genomic data DataFrame
# genotype_df = pd.DataFrame(genotype_data, index=variant_ids, columns=samples).T
# genotype_df.reset_index(inplace=True)
# genotype_df.rename(columns={'index': 'Sample'}, inplace=True)

data = pd.merge(phenotypes, genotype_df, on='Sample')

#get p-values
p_values = []
for variant in variant_ids:
    model = ols(f'Phenotype ~ {variant}', data=data).fit()
    p_value = model.pvalues[1]
    p_values.append(p_value)

#results DataFrame
results = pd.DataFrame({'Variant': variant_ids, 'p_value': p_values})
results['-log10(p_value)'] = -1 * np.log10(results['p_value'])

#plot
plt.scatter(range(len(results)), results['-log10(p_value)'])
plt.xlabel('Variants')
plt.ylabel('-log10(p_value)')
plt.title('Manhattan Plot')
plt.show()
