import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

file_path = "C:/Users/salab/Desktop/UCI_project/PROJECT/Transcriptome/Age/Combined_Expression_Age.csv"
combined_data = pd.read_csv(file_path)

disease_data = combined_data[combined_data['Group'] == 'Disease']

gene_identifiers = combined_data.columns[4:]

correlation_results = []

for gene_id in gene_identifiers:
    gene_expression = disease_data[gene_id]
    age = disease_data['Age']
    correlation, p_value = pearsonr(age, gene_expression)
    correlation_results.append({
        'Gene_ID': gene_id,
        'Correlation': correlation,
        'P-value': p_value
    })

correlation_df = pd.DataFrame(correlation_results)

top_positive_genes = correlation_df.nlargest(5, 'Correlation')
top_negative_genes = correlation_df.nsmallest(5, 'Correlation')

top_genes = pd.concat([top_positive_genes, top_negative_genes])

print(top_genes)

plt.figure(figsize=(10, 6))
plt.barh(top_genes['Gene_ID'], top_genes['Correlation'], color=['red' if x > 0 else 'steelblue' for x in top_genes['Correlation']])
plt.xlabel('Correlation with Age')
plt.ylabel('Gene')

plt.gca().set_yticklabels(top_genes['Gene_ID'], fontsize=10, fontweight='bold')
plt.gca().invert_yaxis()


plt.show()
