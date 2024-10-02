import pandas as pd
from scipy.stats import ttest_ind
import pandas as pd


datasets = [
    {'expression': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE53987\GSE53987_gene_Final.csv',
     'metadata': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE53987\GSE53987_Metadata.csv',
     'disease_label': 'bipolar disorder',
     'control_label': 'control'},
    
    {'expression': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE48350\GSE48350_gene_Final.csv',
     'metadata': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE48350\GSE48350_Metadata.csv',
     'disease_label': 'disease',
     'control_label': 'control'},
    
    {'expression': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE29378\GSE29378_gene_Final.csv',
     'metadata': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE29378\GSE29378_Metadata1.csv',
     'disease_label': 'Disease',
     'control_label': 'Control'},
    
    {'expression': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE28146\GSE28146_gene_Final.csv',
     'metadata': r'C:\Users\salab\Desktop\UCI_project\PROJECT\GSE28146\GSE28146_Metadata.csv',
     'disease_label': 'Disease',
     'control_label': 'Control'}
]


def find_group_column(metadata_df):
    print("Columns in metadata:", metadata_df.columns)  # Debugging: print the columns
    possible_columns = ['Group', 'disease state', 'condition', 'phenotype']
    for col in possible_columns:
        if col in metadata_df.columns:
            print(f"Using column: {col}")  # Debugging: indicate which column is being used
            print(f"Unique values in '{col}':", metadata_df[col].unique())  # Debugging: print unique values in the column
            return col
    raise KeyError("No suitable 'Group' or disease/control identifying column found in the metadata file.")


def perform_differential_expression(expression_path, metadata_path, disease_label, control_label):
    # Load gene expression and metadata
    expression_df = pd.read_csv(expression_path)
    metadata_df = pd.read_csv(metadata_path)


    print("First few rows of metadata:")
    print(metadata_df.head())  # Debugging: print the metadata for inspection
    

    group_col = find_group_column(metadata_df)
    print(f"Using '{group_col}' column for disease/control identification.")

    metadata_df[group_col] = metadata_df[group_col].str.strip()


    disease_samples = metadata_df[metadata_df[group_col].str.lower() == disease_label.lower()]['Sample']
    control_samples = metadata_df[metadata_df[group_col].str.lower() == control_label.lower()]['Sample']

    print(f"Disease Samples: {disease_samples.tolist()}")
    print(f"Control Samples: {control_samples.tolist()}")

    if control_samples.empty or disease_samples.empty:
        raise ValueError(f"No samples found for {disease_label} or {control_label}. Check the metadata labels.")

    disease_expression = expression_df[['Gene_name'] + disease_samples.tolist()]
    control_expression = expression_df[['Gene_name'] + control_samples.tolist()]

    disease_expression.set_index('Gene_name', inplace=True)
    control_expression.set_index('Gene_name', inplace=True)

    fold_change = disease_expression.mean(axis=1) / control_expression.mean(axis=1)
    p_values = ttest_ind(disease_expression, control_expression, axis=1, equal_var=False).pvalue

    return pd.DataFrame({'fold_change': fold_change, 'p_value': p_values})

combined_degs = pd.DataFrame()

for dataset in datasets:
    print(f"Processing dataset: {dataset['expression']}")
    degs = perform_differential_expression(
        dataset['expression'],
        dataset['metadata'],
        dataset['disease_label'],
        dataset['control_label']
    )
    combined_degs = pd.concat([combined_degs, degs])

def perform_differential_expression(expression_path, metadata_path, disease_label, control_label, metadata_col):

    expression_df = pd.read_csv(expression_path)
    metadata_df = pd.read_csv(metadata_path)

    disease_samples = metadata_df[metadata_df[metadata_col].str.lower().str.strip() == disease_label.lower()]['Sample']
    control_samples = metadata_df[metadata_df[metadata_col].str.lower().str.strip() == control_label.lower()]['Sample']

    if control_samples.empty or disease_samples.empty:
        raise ValueError(f"No samples found for {disease_label} or {control_label}. Check the metadata labels.")

    disease_expression = expression_df[['Gene_name'] + disease_samples.tolist()]
    control_expression = expression_df[['Gene_name'] + control_samples.tolist()]

    disease_expression.set_index('Gene_name', inplace=True)
    control_expression.set_index('Gene_name', inplace=True)

    fold_change = disease_expression.mean(axis=1) / control_expression.mean(axis=1)
    p_values = ttest_ind(disease_expression, control_expression, axis=1, equal_var=False).pvalue

    expression_df.set_index('Gene_name', inplace=True)
    expression_df['fold_change'] = fold_change
    expression_df['p_value'] = p_values

    significant_degs = expression_df[((expression_df['fold_change'] > 2.5) | (expression_df['fold_change'] < 0.2)) & (expression_df['p_value'] < 0.05)]
    return significant_degs.reset_index()

datasets = [
    {'expression': r'GSE53987_gene_Final.csv',
     'metadata': r'CGSE53987_Metadata.csv',
     'disease_label': 'bipolar disorder',
     'control_label': 'control',
     'metadata_col': 'disease state'},

    {'expression': r'GSE48350_gene_Final.csv',
     'metadata': r'GSE48350_Metadata.csv',
     'disease_label': 'disease',
     'control_label': 'control',
     'metadata_col': 'Group'},

    {'expression': r'GSE29378_gene_Final.csv',
     'metadata': r'GSE29378_Metadata1.csv',
     'disease_label': 'disease',
     'control_label': 'control',
     'metadata_col': 'Group'},

    {'expression': r'GSE28146_gene_Final.csv',
     'metadata': r'GSE28146_Metadata.csv',
     'disease_label': 'disease',
     'control_label': 'control',
     'metadata_col': 'Group'}
]

combined_degs = pd.DataFrame()

for dataset in datasets:
    print(f"Processing dataset: {dataset['expression']}")
    degs = perform_differential_expression(
        dataset['expression'],
        dataset['metadata'],
        dataset['disease_label'],
        dataset['control_label'],
        dataset['metadata_col']
    )
    combined_degs = pd.concat([combined_degs, degs])

print(combined_degs.head())

combined_degs.to_csv(r'combined_degs.csv', index=False)

cilia_interactome_file_path = r'Cilia_Interactome.csv'
cilia_interactome_df = pd.read_csv(cilia_interactome_file_path)

cilia_genes = set(cilia_interactome_df['Official Symbol Interactor A']).union(set(cilia_interactome_df['Official Symbol Interactor B']))


overlap_genes = combined_degs[combined_degs['Gene_name'].isin(cilia_genes)]

print("Overlapping genes between DEGs and cilia interactome:")
print(overlap_genes.head())

pathway_analysis_path = r'C:\Users\salab\Desktop\UCI_project\PROJECT\Drug_interaction\Pathway_GSE.csv'
pathways_df = pd.read_csv(pathway_analysis_path)


degs_path = r'Significant_Pathways_GeneDrug_Interactions_GSE.csv'
degs_df = pd.read_csv(degs_path)

relevant_pathways = pathways_df['Description'].tolist()

relevant_genes = set()
for genes in pathways_df['geneID']:
    relevant_genes.update(genes.split('/'))

filtered_degs_df = degs_df[degs_df['Gene_name'].isin(relevant_genes)]

gene_drug_path = r'GeneDrug_interactions.csv'
gene_drug_df = pd.read_csv(gene_drug_path)

filtered_gene_drug_df = gene_drug_df[gene_drug_df['Gene_name'].isin(filtered_degs_df['Gene_name'])]


filtered_gene_drug_path = r'Filtered_Significant_Pathways_GeneDrug_Interactions_GSE.csv'
filtered_gene_drug_df = pd.read_csv(filtered_gene_drug_path)

pathway_analysis_path = r'Pathway_GSE.csv'
pathways_df = pd.read_csv(pathway_analysis_path)

relevant_gene_ids = pathways_df['geneID'].str.split('/').sum()

filtered_gene_drug_df = filtered_gene_drug_df[filtered_gene_drug_df['gene_claim_name'].isin(relevant_gene_ids)]

key_genes_path = r'Filtered_Nervous_System_Pathways.csv'
key_genes_df = pd.read_csv(key_genes_path)

key_genes = set(gene for gene_list in key_genes_df['geneID'] for gene in gene_list.split('/'))

gene_drug_interaction_path = r'GeneDrug_interactions.csv'
gene_drug_df = pd.read_csv(gene_drug_interaction_path)

nervous_system_drug_interactions = gene_drug_df[gene_drug_df['Gene_name'].isin(key_genes)]



key_genes_path = r'Significant_Nervous_System_Pathways.csv'
key_genes_df = pd.read_csv(key_genes_path)

key_genes = set(gene for gene_list in key_genes_df['geneID'] for gene in gene_list.split('/'))

gene_drug_interaction_path = r'GeneDrug_interactions.csv'
gene_drug_df = pd.read_csv(gene_drug_interaction_path)

nervous_system_drug_interactions = gene_drug_df[gene_drug_df['Gene_name'].isin(key_genes)]

approved_drug_interactions = nervous_system_drug_interactions[nervous_system_drug_interactions['approved'] == True]









