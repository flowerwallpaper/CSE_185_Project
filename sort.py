import pandas as pd

# Load the GWAS results file
gwas_results = pd.read_csv('lab3_gwas_covar.assoc.linear', sep='\s+')

# Sort the DataFrame by the p-values in ascending order
sorted_gwas_results = gwas_results.sort_values(by='P')

# Save the sorted DataFrame to a new file
sorted_gwas_results.to_csv('sorted_gwas_results.txt', sep='\t', index=False)

print("GWAS results sorted by p-value and saved to 'sorted_gwas_results.txt'")
