import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import math


comprehensive_list_file_path = r'comprehensive_list.csv'
comprehensive_df = pd.read_csv(comprehensive_list_file_path)


gene_drug_file_path = r'GeneDrug_interactions.csv'
gene_drug_df = pd.read_csv(gene_drug_file_path)


comprehensive_genes = comprehensive_df['Gene'].tolist()
filtered_gene_drug_df = gene_drug_df[gene_drug_df['Gene_name'].isin(comprehensive_genes)]

nervous_system_pathways_file_path = r'Nervous_System_Pathways.csv'
nervous_system_pathways_df = pd.read_csv(nervous_system_pathways_file_path)

significant_genes_in_nervous_system_pathways = set()
for genes in nervous_system_pathways_df['geneID']:
    significant_genes_in_nervous_system_pathways.update(genes.split('/'))

nervous_system_drug_interactions_df = filtered_gene_drug_df[filtered_gene_drug_df['Gene_name'].isin(significant_genes_in_nervous_system_pathways)]


nervous_system_drug_interactions_df.to_csv(r'Nervous_System_Drug_Interactions.csv', index=False)

file_path = r'\Nervous_System_Drug_Interactions.csv'
df = pd.read_csv(file_path)

# Filter out genes that don't have any drug interactions and where "approved" is TRUE
df_filtered = df[(df['drug_name'].notna()) & (df['approved'] == True)]

# Create the graph
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


pos = nx.spring_layout(G, k=0.15)  # Decrease 'k' to shorten distances between nodes


def adjust_positions(pos, adjustment_nodes, adjustment_factor=1.1):
    adjusted_pos = {}
    for node, (x, y) in pos.items():
        if node in adjustment_nodes:
            adjusted_pos[node] = (x * adjustment_factor, y * adjustment_factor)
        else:
            adjusted_pos[node] = (x, y)
    return adjusted_pos

adjustment_nodes = ['4-PHENYLBUTYRIC ACID', 'BUTYRYLHYDROXAMIC ACID']
pos = adjust_positions(pos, adjustment_nodes, adjustment_factor=1.2)

label_pos = {node: (pos[node][0], pos[node][1] + 0.1) for node in G.nodes()}  # Adjust the label position to be above the node


plt.figure(figsize=(8, 8))  # Increase figure size for better readability
nx.draw(G, pos, with_labels=False, node_color=node_colors, edge_color=edge_colors, node_size=1000, font_size=12, font_weight='bold')  # Increase node_size and font_size
nx.draw_networkx_labels(G, label_pos, font_size=12, font_weight='bold')
plt.title('Gene-Drug Interaction Network')


output_path = r'C:\Users\salab\Desktop\UCI_project\PROJECT\Drug_interaction\APRROVE_gene_drug_network.png'
plt.savefig(output_path, dpi=300)

plt.show()
