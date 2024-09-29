import pandas as pd

#import proact scores for rocesses
proact_scores = pd.read_csv(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\ProAct scores\gtex_score_matrix.csv")
processes_and_genes = pd.read_csv(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\ProAct scores\processes_with_genes_and_gos_v8.csv")
#melt the processes and genes dataframe for genes in the 'genes' column that are seperated by ';'
df_genes_expanded = processes_and_genes.assign(Gene=processes_and_genes['genes'].str.split(';')).explode('Gene')
df_genes_expanded['Gene'] = df_genes_expanded['Gene'].str.strip()

#Load mitocarta ENSGs
mitocarta = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Human.MitoCarta3.0.xls", sheet_name='A Human MitoCarta3.0')
mitocarta['EnsemblGeneID_mapping_version_20200130'] = mitocarta['EnsemblGeneID_mapping_version_20200130'].str.split('|').str[0]
mitocarta_ENSGs = list(mitocarta['EnsemblGeneID_mapping_version_20200130'].unique())
intersection_genes = set(df_genes_expanded['Gene']).intersection(set(mitocarta_ENSGs))

df_genes_filtered = df_genes_expanded[df_genes_expanded['Gene'].isin(mitocarta_ENSGs)]
nested_dict = {}
for gene in df_genes_filtered['Gene'].unique():
    gene_dict = {}
    gene_processes = df_genes_filtered[df_genes_filtered['Gene'] == gene]

    for _, row in gene_processes.iterrows():
        process = row['process']
        tissue_scores = proact_scores[proact_scores['Process'] == process].drop(columns=['Process']).to_dict(orient='index')

        # Get the first item from the dictionary
        if tissue_scores:
            gene_dict[process] = list(tissue_scores.values())[0]
        else:
            gene_dict[process] = {}

    nested_dict[gene] = gene_dict

tissue_columns_with_mean = [tissue + '-mean_proact' for tissue in proact_scores.columns.drop('Process')]
mean_scores_df = pd.DataFrame(index=mitocarta_ENSGs, columns=tissue_columns_with_mean)

for gene, processes_dict in nested_dict.items():
    for tissue in proact_scores.columns.drop('Process'):
        tissue_scores = []
        for process, tissue_dict in processes_dict.items():
            if tissue in tissue_dict:
                tissue_scores.append(tissue_dict[tissue])

        # Calculate the mean score for this tissue
        if tissue_scores:
            mean_scores_df.loc[gene, tissue + '-mean_proact'] = sum(tissue_scores) / len(tissue_scores)
        else:
            mean_scores_df.loc[gene, tissue + '-mean_proact'] = None  # Or 0 if you prefer to fill missing values with zero
#import to csv
mean_scores_df.to_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\ProAct_mean_scores_genes causality.csv")
#create proact mean scores file for mitomodel tissues
mitoModel_tissues_associations = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx", sheet_name='Sheet1')
tissues_dict = {}
for _, row in mitoModel_tissues_associations.iterrows():
    if pd.notna(row['mean_proact_scores']):
        tissues_dict[row['mitoModel']] = row['mean_proact_scores'].replace(' ','').split('|')

tissues_mitomodel_with_mean = [tissue+ '_mean proact'for tissue in tissues_dict.keys()]
#match the tissues from the mitomodel to the mean_scores_df dataframe, if a mitomodel tissue has more than one assiciated proact tissue, take the mean of the proact tissues

mean_scores_df_mitomodel = pd.DataFrame(index=mitocarta_ENSGs, columns=tissues_mitomodel_with_mean)
for mitomodel_tissue, proact_tissues in tissues_dict.items():
    proact_tissues = [tissue for tissue in proact_tissues]
    mean_scores_df_mitomodel[mitomodel_tissue + '_mean proact'] = mean_scores_df[proact_tissues].mean(axis=1, skipna=True)

mean_scores_df_mitomodel.to_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\ProAct_mean_scores_causality in tissues.csv")





print('meow')

