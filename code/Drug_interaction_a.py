import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import math


drugbank_path = r"DrugBank.csv"
comprehensive_list_path = r"comprehensive_list.csv"
filtered_output_path = r"Pathways_GeneDrug_Interactions.csv"
nervous_system_pathways_file_path = r'Nervous_System_Pathways.csv'
nervous_system_drug_interactions_output_path = r'Nervous_System_Drug_Interactions.csv'


drugbank_df = pd.read_csv(drugbank_path)
comprehensive_list_df = pd.read_csv(comprehensive_list_path)

filtered_df = drugbank_df[drugbank_df['Gene_name'].isin(comprehensive_list_df['Gene'])]


filtered_df.to_csv(filtered_output_path, index=False)
print(f"Filtered data saved to: {filtered_output_path}")


filtered_gene_drug_df = pd.read_csv(filtered_output_path)


nervous_system_pathways_df = pd.read_csv(nervous_system_pathways_file_path)


significant_genes_in_nervous_system_pathways = set()
for genes in nervous_system_pathways_df['geneID']:
    significant_genes_in_nervous_system_pathways.update(genes.split('/'))


nervous_system_drug_interactions_df = filtered_gene_drug_df[
    filtered_gene_drug_df['Gene_name'].isin(significant_genes_in_nervous_system_pathways)
]


nervous_system_drug_interactions_df.to_csv(nervous_system_drug_interactions_output_path, index=False)
print(f"Nervous system drug interactions saved to: {nervous_system_drug_interactions_output_path}")


df = pd.read_csv(nervous_system_drug_interactions_output_path)

# Filter out genes that don't have any drug interactions and where "approved" is TRUE
df_filtered = df[(df['drug_name'].notna()) & (df['approved'] == True)]

G = nx.Graph()


for _, row in df_filtered.iterrows():
    gene = row['Gene_name']
    drug = row['drug_name']
    G.add_node(gene, type='gene', color='blue')
    G.add_node(drug, type='drug', color='red')
    G.add_edge(gene, drug, color='lightgray')


to_remove = [node for node, degree in dict(G.degree()).items() if degree == 0]
G.remove_nodes_from(to_remove)


node_colors = [G.nodes[node]['color'] for node in G]
edge_colors = [G.edges[edge]['color'] for edge in G.edges()]


pos = nx.spring_layout(G, k=0.15)

def adjust_positions(pos, adjustment_nodes, adjustment_factor=1.1):
    adjusted_pos = {}
    for node, (x, y) in pos.items():
        if node in adjustment_nodes:
            angle = math.atan2(y, x)
            adjusted_pos[node] = (x * adjustment_factor, y * adjustment_factor)
        else:
            adjusted_pos[node] = (x, y)
    return adjusted_pos

adjustment_nodes = ['4-PHENYLBUTYRIC ACID', 'BUTYRYLHYDROXAMIC ACID']
pos = adjust_positions(pos, adjustment_nodes, adjustment_factor=1.2)

label_pos = {node: (pos[node][0], pos[node][1] + 0.1) for node in G.nodes()}

plt.figure(figsize=(8, 8))  # Increase figure size for better readability
nx.draw(G, pos, with_labels=False, node_color=node_colors, edge_color=edge_colors, node_size=1000, font_size=12, font_weight='bold')
nx.draw_networkx_labels(G, label_pos, font_size=12, font_weight='bold')


plt.title('Gene-Drug Interaction Network')
plt.savefig(graph_output_path, dpi=300)
plt.show()


