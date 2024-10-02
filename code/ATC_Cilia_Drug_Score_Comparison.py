import pandas as pd
import matplotlib.pyplot as plt

# Load the key genes in nervous system pathways
key_genes_path = r'Significant_Nervous_System_Pathways.csv'
key_genes_df = pd.read_csv(key_genes_path)

# Extract unique genes from the geneID column
key_genes = set(gene for gene_list in key_genes_df['geneID'] for gene in gene_list.split('/'))

# Load the Gene-Drug interaction data
gene_drug_interaction_path = r'C:\Users\salab\Desktop\UCI_project\PROJECT\GeneDrug_interactions.csv'
gene_drug_df = pd.read_csv(gene_drug_interaction_path)

# Filter for key genes in the nervous system pathways
nervous_system_drug_interactions = gene_drug_df[gene_drug_df['Gene_name'].isin(key_genes)]

# Further filter to include only approved drugs
approved_drug_interactions = nervous_system_drug_interactions[nervous_system_drug_interactions['approved'] == True]

# Display the filtered gene-drug interactions with approved drugs
print("Filtered Gene-Drug interaction data with approved drugs:")
print(approved_drug_interactions.head())

non_cilia_df = non_cilia_df[non_cilia_df['approved'] == True]
categorized_significant_pathways_df = categorized_significant_pathways_df[categorized_significant_pathways_df['approved'] == True]

comparison_drugs = [
    'NORTRIPTYLINE', 'FLUSPIRILENE', 'HEROIN', 
    'METHYLPHENIDATE HYDROCHLORIDE', 'VALPROIC ACID', 
    'BRAIN-DERIVED NEUROTROPHIC FACTOR', 
    'LITHIUM CITRATE, LITHIUM CARBONATE', 'TRIFLUOPERAZINE', 
    'CITALOPRAM', 'PROPOFOL'
]

focus_atc_categories = [
    'Nervous system', 'Opioid Analgesic', 'Centrally acting sympathomimetics', 
    'Antiepileptics', 'Neurotrophic Factor', 'Antipsychotics', 
    'Psychoanaleptics, Antidepressants, SSRIs', 'General anesthetics'
]

non_cilia_filtered = non_cilia_df[
    non_cilia_df['drug_name'].isin(comparison_drugs) & 
    non_cilia_df['ATC_Category'].isin(focus_atc_categories)
].copy()

cilia_filtered = categorized_significant_pathways_df[
    categorized_significant_pathways_df['drug_name'].isin(comparison_drugs) & 
    categorized_significant_pathways_df['ATC_Category'].isin(focus_atc_categories)
].copy()

non_cilia_filtered.loc[:, 'group'] = 'Non-Cilia'
cilia_filtered.loc[:, 'group'] = 'Cilia'

combined_df = pd.concat([non_cilia_filtered, cilia_filtered])

avg_interaction_scores = combined_df.groupby(['drug_name', 'group', 'ATC_Category'])['interaction_score'].mean().reset_index()

avg_interaction_scores_pivot = avg_interaction_scores.pivot(
    index=['drug_name', 'ATC_Category'], 
    columns='group', 
    values='interaction_score'
).fillna(0).reset_index()

avg_interaction_scores_pivot['label'] = avg_interaction_scores_pivot['drug_name'] + '\n(' + avg_interaction_scores_pivot['ATC_Category'] + ')'

avg_interaction_scores_pivot['Cilia'] = avg_interaction_scores_pivot['Cilia'].fillna(0)
avg_interaction_scores_pivot['Non-Cilia'] = avg_interaction_scores_pivot['Non-Cilia'].fillna(0)

avg_interaction_scores_pivot.set_index('label', inplace=True)

ax = avg_interaction_scores_pivot[['Cilia', 'Non-Cilia']].plot(
    kind='bar', 
    stacked=True, 
    color=['blue', 'red'], 
    figsize=(12, 8)
)

plt.xlabel('Drug Name (ATC Category)')
plt.ylabel('Average Interaction Score')
plt.xticks(rotation=45, ha='right', fontsize=15)
plt.tight_layout()

plt.show()

