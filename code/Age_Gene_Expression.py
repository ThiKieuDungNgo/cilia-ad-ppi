import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import seaborn as sns

# Load the combined data
file_path = "Combined_Expression_Age.csv"
combined_data = pd.read_csv(file_path)


gene_identifiers = combined_data.columns[4:]  # Assuming first four columns are metadata (Sample, Gender, Group, Age)
combined_data[gene_identifiers] = combined_data[gene_identifiers].apply(pd.to_numeric, errors='coerce')


control_data = combined_data[combined_data['Group'] == 'Control']
disease_data = combined_data[combined_data['Group'] == 'Disease']


correlation_results = []
trend_results = []
impact_results = []

for gene in gene_identifiers:
    gene_expression_disease = disease_data[gene]
    gene_expression_control = control_data[gene]
    

    correlation, p_value_corr = pearsonr(disease_data['Age'], gene_expression_disease)
    correlation_results.append({
        'Gene_ID': gene,
        'Correlation': correlation,
        'P-value': p_value_corr
    })
    
   
    disease_X = sm.add_constant(disease_data['Age'])
    disease_model = sm.OLS(gene_expression_disease, disease_X).fit()
    disease_slope = disease_model.params['Age']
    disease_pvalue = disease_model.pvalues['Age']
    

    control_X = sm.add_constant(control_data['Age'])
    control_model = sm.OLS(gene_expression_control, control_X).fit()
    control_slope = control_model.params['Age']
    control_pvalue = control_model.pvalues['Age']
    

    trend_results.append({
        'Gene': gene,
        'Control_Slope': control_slope,
        'Control_P-value': control_pvalue,
        'Disease_Slope': disease_slope,
        'Disease_P-value': disease_pvalue,
        'R-squared': disease_model.rsquared
    })


correlation_df = pd.DataFrame(correlation_results)
top_positive_genes = correlation_df.nlargest(5, 'Correlation')
top_negative_genes = correlation_df.nsmallest(5, 'Correlation')
top_genes = pd.concat([top_positive_genes, top_negative_genes])


plt.figure(figsize=(10, 6))
plt.barh(top_genes['Gene_ID'], top_genes['Correlation'], color=['red' if x > 0 else 'steelblue' for x in top_genes['Correlation']])
plt.xlabel('Correlation with Age')
plt.ylabel('Gene')
plt.gca().invert_yaxis()
plt.tight_layout()
save_path = 'C:/Users/salab/Desktop/UCI_project/PROJECT/Transcriptome/Age/Top_Genes_Correlation_with_Age.png'
plt.savefig(save_path)
plt.show()


trend_df = pd.DataFrame(trend_results)


trend_df['Control_Adj_P-value'] = multipletests(trend_df['Control_P-value'], method='fdr_bh')[1]
trend_df['Disease_Adj_P-value'] = multipletests(trend_df['Disease_P-value'], method='fdr_bh')[1]


significant_control_trends = trend_df[trend_df['Control_Adj_P-value'] < 0.05]
significant_disease_trends = trend_df[trend_df['Disease_Adj_P-value'] < 0.05]


print("Significant Control Trends:\n", significant_control_trends)
print("Significant Disease Trends:\n", significant_disease_trends)


for gene in significant_disease_trends['Gene']:
    # Fit linear model for disease group
    disease_X = sm.add_constant(disease_data['Age'])
    disease_model = sm.OLS(disease_data[gene], disease_X).fit()
    disease_slope = disease_model.params['Age']
    disease_pvalue = disease_model.pvalues['Age']
    
    # Store the impact results
    impact_results.append({
        'Gene': gene,
        'Disease_Slope': disease_slope,
        'Disease_P-value': disease_pvalue,
        'R-squared': disease_model.rsquared
    })


impact_df = pd.DataFrame(impact_results)
print("Impact of Gene Expression on Alzheimer's Progression:\n", impact_df)


high_r_squared_genes = impact_df[impact_df['R-squared'] > 0.7]


if not high_r_squared_genes.empty:
    for gene in high_r_squared_genes['Gene']:
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='Age', y=gene, hue='Group', data=combined_data)
        sns.regplot(x='Age', y=gene, data=disease_data, scatter=False, color='red', label='Disease Trend')
        plt.xlabel('Age')
        plt.ylabel('Expression Level')
        plt.legend(title='Group')
        filename = f"C:/Users/salab/Desktop/UCI_project/PROJECT/Transcriptome/{gene}_expression_plot.png"
        plt.savefig(filename)
        plt.show()
else:
    print("No genes with high R-squared values found to visualize.")
