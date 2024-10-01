
cilia_genes_df = pd.read_csv(cilia_genes_path)
ppis_df = pd.read_csv(ppis_path, dtype={30: str, 31: str, 32: str})
alzheimers_genes_df = pd.read_csv(alzheimers_genes_path)
geneage_human_df = pd.read_csv(geneage_human_path)
hedgehog_pathway_df = pd.read_csv(hedgehog_pathway_path)
neurogenesis_df = pd.read_csv(neurogenesis_path)

cilia_genes_set = set(cilia_genes_df['Gene_name'])

alzheimers_genes_set = set(alzheimers_genes_df['mappedGenes1']).union(set(alzheimers_genes_df['mappedGenes2']))

# Extract gene names from other datasets
geneage_human_genes_set = set(geneage_human_df['Gene_name'])
hedgehog_genes_set = set(hedgehog_pathway_df['Gene_name'])
neurogenesis_genes_set = set(neurogenesis_df['Gene_name'])
filtered_cilia_genes_set = cilia_genes_set.intersection(
    alzheimers_genes_set.union(
        geneage_human_genes_set,
        hedgehog_genes_set,
        neurogenesis_genes_set
    )
)

filtered_cilia_genes_set = set(filtered_cilia_genes_set)
print(len(filtered_cilia_genes_set))
filtered_ppis = ppis_df[
    (ppis_df['Official Symbol Interactor A'].isin(filtered_cilia_genes_set)) &
    (ppis_df['Official Symbol Interactor B'].isin(filtered_cilia_genes_set))
]

positive_examples = filtered_ppis.copy()
positive_examples['label'] = 1

# Extract all unique genes from the PPIs dataset
all_genes = set(ppis_df['Official Symbol Interactor A']).union(set(ppis_df['Official Symbol Interactor B']))
all_genes_list = list(all_genes)

negative_pairs = []
while len(negative_pairs) < len(positive_examples):
    gene_a, gene_b = random.sample(all_genes_list, 2)
    if gene_a != gene_b and not (
        ((ppis_df['Official Symbol Interactor A'] == gene_a) & (ppis_df['Official Symbol Interactor B'] == gene_b)).any() or
        ((ppis_df['Official Symbol Interactor A'] == gene_b) & (ppis_df['Official Symbol Interactor B'] == gene_a)).any()
    ):
        negative_pairs.append([gene_a, gene_b])

negative_examples = pd.DataFrame(negative_pairs, columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])
negative_examples['label'] = 0

data = pd.concat([positive_examples[['Official Symbol Interactor A', 'Official Symbol Interactor B', 'label']], negative_examples])

data = data.sample(frac=1).reset_index(drop=True)

all_genes = list(set(data['Official Symbol Interactor A']).union(set(data['Official Symbol Interactor B'])).union(cilia_genes_set))
le = LabelEncoder()
le.fit(all_genes)
data['Interactor A Encoded'] = le.transform(data['Official Symbol Interactor A'])
data['Interactor B Encoded'] = le.transform(data['Official Symbol Interactor B'])

X = data[['Interactor A Encoded', 'Interactor B Encoded']]
y = data['label']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

print(f"Shape of X_train: {X_train.shape}")
print(f"Shape of X_test: {X_test.shape}")
print(f"Shape of y_train: {y_train.shape}")
print(f"Shape of y_test: {y_test.shape}")

cv_scores = cross_val_score(RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced'), X_train, y_train, cv=5)

model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

from sklearn.metrics import precision_recall_curve, auc

y_pred_proba = model.predict_proba(X_test)[:, 1]
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
pr_auc = auc(recall, precision)


def predict_in_batches(model, le, filtered_cilia_genes_set, all_genes, batch_size=100000, threshold=1):
    all_possible_pairs = []
    novel_ppis_list = []
    batch_count = 0

    for cilia_gene in filtered_cilia_genes_set:
        for gene in all_genes:
            if cilia_gene != gene:
                all_possible_pairs.append([cilia_gene, gene])

                # Process in batches
                if len(all_possible_pairs) == batch_size:
                    batch_count += 1
                    print(f"Processing batch {batch_count}")
                    possible_pairs_df = pd.DataFrame(all_possible_pairs, columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])
                    possible_pairs_df['Interactor A Encoded'] = le.transform(possible_pairs_df['Official Symbol Interactor A'])
                    possible_pairs_df['Interactor B Encoded'] = le.transform(possible_pairs_df['Official Symbol Interactor B'])

                    # Predict novel PPIs
                    possible_X = possible_pairs_df[['Interactor A Encoded', 'Interactor B Encoded']]
                    possible_pairs_df['Prediction_Prob'] = model.predict_proba(possible_X)[:, 1]
                    possible_pairs_df['Prediction'] = (possible_pairs_df['Prediction_Prob'] >= threshold).astype(int)

                    # Collect novel PPIs
                    novel_ppis_batch = possible_pairs_df[possible_pairs_df['Prediction'] == 1]
                    novel_ppis_list.append(novel_ppis_batch)

                    # Clear the list for the next batch
                    all_possible_pairs.clear()

    if all_possible_pairs:
        batch_count += 1
        print(f"Processing batch {batch_count}")
        possible_pairs_df = pd.DataFrame(all_possible_pairs, columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])
        possible_pairs_df['Interactor A Encoded'] = le.transform(possible_pairs_df['Official Symbol Interactor A'])
        possible_pairs_df['Interactor B Encoded'] = le.transform(possible_pairs_df['Official Symbol Interactor B'])

        possible_X = possible_pairs_df[['Interactor A Encoded', 'Interactor B Encoded']]
        possible_pairs_df['Prediction_Prob'] = model.predict_proba(possible_X)[:, 1]
        possible_pairs_df['Prediction'] = (possible_pairs_df['Prediction_Prob'] >= threshold).astype(int)

        novel_ppis_batch = possible_pairs_df[possible_pairs_df['Prediction'] == 1]
        novel_ppis_list.append(novel_ppis_batch)

    novel_ppis = pd.concat(novel_ppis_list, ignore_index=True)
    return novel_ppis

novel_ppis = predict_in_batches(model, le, filtered_cilia_genes_set, all_genes, batch_size=100000, threshold=1)

known_ppis = filtered_ppis[['Official Symbol Interactor A', 'Official Symbol Interactor B']].copy()
novel_ppis = novel_ppis[['Official Symbol Interactor A', 'Official Symbol Interactor B']].copy()
known_ppis.loc[:, 'Type'] = 'Known'
novel_ppis.loc[:, 'Type'] = 'Novel'
interactome = pd.concat([known_ppis, novel_ppis]).drop_duplicates()

nodes = []

for gene in filtered_cilia_genes_set:
    if gene in set(known_ppis['Official Symbol Interactor A']).union(set(known_ppis['Official Symbol Interactor B'])):
        nodes.append({'Gene': gene, 'Type': 'Ciliary Known'})
    elif gene in set(novel_ppis['Official Symbol Interactor A']).union(set(novel_ppis['Official Symbol Interactor B'])):
        nodes.append({'Gene': gene, 'Type': 'Ciliary Novel'})
    else:
        nodes.append({'Gene': gene, 'Type': 'Ciliary'})

for gene in set(known_ppis['Official Symbol Interactor A']).union(set(known_ppis['Official Symbol Interactor B'])):
    if gene not in filtered_cilia_genes_set:
        nodes.append({'Gene': gene, 'Type': 'Known'})

for gene in set(novel_ppis['Official Symbol Interactor A']).union(set(novel_ppis['Official Symbol Interactor B'])):
    if gene not in filtered_cilia_genes_set:
        nodes.append({'Gene': gene, 'Type': 'Novel'})

nodes_df = pd.DataFrame(nodes)


X = known_ppis.shape[0]
Y = novel_ppis.shape[0]
Z = len(set(known_ppis['Official Symbol Interactor A']).union(set(known_ppis['Official Symbol Interactor B'])))
A = len(set(novel_ppis['Official Symbol Interactor A']).union(set(novel_ppis['Official Symbol Interactor B'])))
B = novel_ppis.shape[0]
C = len(set(novel_ppis['Official Symbol Interactor A']).union(set(novel_ppis['Official Symbol Interactor B'])).intersection(filtered_cilia_genes_set))
D = len(filtered_cilia_genes_set.difference(set(known_ppis['Official Symbol Interactor A']).union(set(known_ppis['Official Symbol Interactor B']))))

results = []

for gene in filtered_cilia_genes_set:
    known_interactors = set()
    novel_interactors = set()

    known_interactors.update(known_ppis[known_ppis['Official Symbol Interactor A'] == gene]['Official Symbol Interactor B'])
    known_interactors.update(known_ppis[known_ppis['Official Symbol Interactor B'] == gene]['Official Symbol Interactor A'])

    novel_interactors.update(novel_ppis[novel_ppis['Official Symbol Interactor A'] == gene]['Official Symbol Interactor B'])
    novel_interactors.update(novel_ppis[novel_ppis['Official Symbol Interactor B'] == gene]['Official Symbol Interactor A'])

    known_count = len(known_interactors)
    novel_count = len(novel_interactors)
    
    known_interactors = sorted(known_interactors)
    novel_interactors = sorted(novel_interactors)

    novel_interactors_str = ', '.join(novel_interactors)
    
    results.append([gene, known_count, novel_count, novel_interactors_str])

results_df = pd.DataFrame(results, columns=['Gene', 'Known PPIs', 'Novel PPIs', 'Novel Interactors'])
