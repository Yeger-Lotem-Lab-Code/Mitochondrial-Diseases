import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
#bio_pro
json_file_bio_pro = r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Panther enrichment analysis - mito vs prot cod\bio_pro\analysis.json"
# Open and load the JSON file
with open(json_file_bio_pro, 'r') as file:
    data = json.load(file)
enrichment_bio_pro = data['overrepresentation']['result']
sorted_enrichment_bio_pro = sorted(enrichment_bio_pro, key=lambda x: x['input_list']['pValue'])
top_10_lowest_pValue = sorted_enrichment_bio_pro[1:11]
enrichment_data = []
#create a df with columns term, pvalue and mapped IDs
for dicti in top_10_lowest_pValue:
    pval = dicti['input_list']['pValue']
    term = dicti['term']['label']
    mapped_ids = dicti['input_list']['mapped_id_list']['mapped_id']
    #append to the df
    enrichment_data.append({'Term': term, 'pValue': pval, 'Genes': mapped_ids})
enrichment_bio_pro = pd.DataFrame(enrichment_data)

#cell_comp
json_file_cell_comp = r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Panther enrichment analysis - mito vs prot cod\cellu_comp\analysis.json"
# Open and load the JSON file
with open(json_file_cell_comp, 'r') as file:
    data = json.load(file)
enrichment_cell_comp = data['overrepresentation']['result']
sorted_enrichment_cell_comp = sorted(enrichment_cell_comp, key=lambda x: x['input_list']['pValue'])
#remove the third element in the list
sorted_enrichment_cell_comp = sorted_enrichment_cell_comp[:2] + sorted_enrichment_cell_comp[3:]
top_10_lowest_pValue = sorted_enrichment_cell_comp[:10]
enrichment_data = []
#create a df with columns term, pvalue and mapped IDs
for dicti in top_10_lowest_pValue:
    pval = dicti['input_list']['pValue']
    term = dicti['term']['label']
    mapped_ids = dicti['input_list']['mapped_id_list']['mapped_id']
    #append to the df
    enrichment_data.append({'Term': term, 'pValue': pval, 'Genes': mapped_ids})
enrichmennt_cell_comp = pd.DataFrame(enrichment_data)

#mol_fun
json_file_mol_fun = r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Panther enrichment analysis - mito vs prot cod\molec_func\analysis.json"
# Open and load the JSON file
with open(json_file_mol_fun, 'r') as file:
    data = json.load(file)
enrichment_mol_fun = data['overrepresentation']['result']
sorted_enrichment_mol_fun = sorted(enrichment_mol_fun, key=lambda x: x['input_list']['pValue'])
top_10_lowest_pValue = sorted_enrichment_mol_fun[1:11]
enrichment_data = []
#create a df with columns term, pvalue and mapped IDs
for dicti in top_10_lowest_pValue:
    pval = dicti['input_list']['pValue']
    term = dicti['term']['label']
    mapped_ids = dicti['input_list']['mapped_id_list']['mapped_id']
    #append to the df
    enrichment_data.append({'Term': term, 'pValue': pval, 'Genes': mapped_ids})
enrichment_mol_fun = pd.DataFrame(enrichment_data)


#for each enrichment table, create a new df that has the number of genes associate to each term
enrichment_bio_pro['num_genes'] = enrichment_bio_pro['Genes'].apply(len)
enrichmennt_cell_comp['num_genes'] = enrichmennt_cell_comp['Genes'].apply(len)
enrichment_mol_fun['num_genes'] = enrichment_mol_fun['Genes'].apply(len)



#plot histogram for the top 10 enriched terms in each category
plt.figure(figsize=(20, 20))
sns.barplot(data=enrichment_bio_pro.head(10), x='Term', y='num_genes', palette='viridis')
#y and x titles fontsize
plt.xlabel('Term', fontsize=20)
plt.ylabel('Number of genes', fontsize=20)
plt.title('Top 10 enriched GO terms in Biological Process', fontsize=20)
plt.xticks(rotation=45, fontsize=16, ha='right')
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig('bio_pro_enrichment.png')
plt.show()

plt.figure(figsize=(20, 20))
sns.barplot(data=enrichmennt_cell_comp.head(10), x='Term', y='num_genes', palette='viridis')
plt.xlabel('Term', fontsize=20)
plt.ylabel('Number of genes', fontsize=20)
plt.title('Top 10 enriched GO terms in Cellular Component', fontsize=20)
plt.xticks(rotation=45, fontsize=16, ha='right')
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig('cell_comp_enrichment.png')
plt.show()

plt.figure(figsize=(20, 20))
sns.barplot(data=enrichment_mol_fun.head(10), x='Term', y='num_genes', palette='viridis')
plt.xlabel('Term', fontsize=20)
plt.ylabel('Number of genes', fontsize=20)
plt.title('Top 10 enriched GO terms in Molecular Function',fontsize=20)
plt.xticks(rotation=45, fontsize=16, ha='right')
plt.yticks(fontsize=16)
plt.tight_layout(pad=3.0)
plt.savefig('mol_fun_enrichment.png')
plt.show()

#create a df with mitocarta genes as rows and top 10 enriched terms as columns
mitocarta_genes = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Human.MitoCarta3.0.xls", sheet_name='A Human MitoCarta3.0')[['Symbol', 'EnsemblGeneID_mapping_version_20200130']]
mitocarta_genes.columns = ['Gene', 'ENSG']
top10_bio_pro = enrichment_bio_pro['Term'].apply(lambda x: 'biological_process_' + x).tolist()
top10_cell_comp = enrichmennt_cell_comp['Term'].apply(lambda x: 'cellular_component_' + x).tolist()
top10_mol_fun = enrichment_mol_fun['Term'].apply(lambda x: 'molecular_function_' + x).tolist()
top_terms = pd.concat([pd.Series(top10_bio_pro), pd.Series(top10_cell_comp), pd.Series(top10_mol_fun)]).unique()
GO_terms_enriched_df = pd.DataFrame(0, index=range(len(mitocarta_genes)), columns=['Gene', 'ENSG'] + list(top_terms))
GO_terms_enriched_df['Gene'] = mitocarta_genes['Gene']
GO_terms_enriched_df['ENSG'] = mitocarta_genes['ENSG']

def update_enriched_df(enrichment_df, go_terms_df, aspect_prefix):
    for term in enrichment_df['Term']:
        genes = enrichment_df[enrichment_df['Term'] == term]['Genes'].values[0]
        column_name = aspect_prefix + term
        for gene in genes:
            if gene in go_terms_df['Gene'].values:
                go_terms_df.loc[go_terms_df['Gene'] == gene, column_name] = 1
    return go_terms_df

GO_terms_enriched_df = update_enriched_df(enrichment_bio_pro, GO_terms_enriched_df, 'biological_process_')
GO_terms_enriched_df = update_enriched_df(enrichmennt_cell_comp, GO_terms_enriched_df, 'cellular_component_')
GO_terms_enriched_df = update_enriched_df(enrichment_mol_fun, GO_terms_enriched_df, 'molecular_function_')


#GO_terms_enriched_df.to_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Panther enrichment analysis - mito vs prot cod\GO_terms_enriched_mitocarta_genes vs prot cod.csv", index=False)

# Load causal genes data
df = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_pref_updated.csv")
columns_to_keep = ['ENSG ID', 'Label'] + [col for col in df.columns if any(term in col for term in ['biological_process', 'cellular_component', 'molecular_function'])]
GO_terms_df = df[columns_to_keep]

melted_df = GO_terms_df.melt(id_vars=['ENSG ID', 'Label'], var_name='Term', value_name='Connected')
melted_df = melted_df[melted_df['Connected'] == 1]

melted_df['Aspect'] = melted_df['Term'].apply(lambda x: 'Biological Process' if 'biological_process' in x else
                                                       'Cellular Component' if 'cellular_component' in x else
                                                       'Molecular Function' if 'molecular_function' in x else 'Other')

# Separate dataframes for each aspect
biological_process_df = melted_df[melted_df['Aspect'] == 'Biological Process']
cellular_component_df = melted_df[melted_df['Aspect'] == 'Cellular Component']
molecular_function_df = melted_df[melted_df['Aspect'] == 'Molecular Function']

# Plotting
plt.figure(figsize=(12, 8))
sns.barplot(data=biological_process_df, x='Term', y='Connected', hue='Label', estimator=len)
plt.xlabel('GO Terms')
plt.ylabel('Number of Genes')
plt.title('Number of Genes Connected to Each GO Term (Biological Process)')
plt.legend(title='Gene Type', labels=['Non-Causal', 'Causal'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 8))
sns.barplot(data=cellular_component_df, x='Term', y='Connected', hue='Label', estimator=len)
plt.xlabel('GO Terms')
plt.ylabel('Number of Genes')
plt.title('Number of Genes Connected to Each GO Term (Cellular Component)')
plt.legend(title='Gene Type', labels=['Non-Causal', 'Causal'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 8))
sns.barplot(data=molecular_function_df, x='Term', y='Connected', hue='Label', estimator=len)
plt.xlabel('GO Terms')
plt.ylabel('Number of Genes')
plt.title('Number of Genes Connected to Each GO Term (Molecular Function)')
plt.legend(title='Gene Type', labels=['Non-Causal', 'Causal'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()




print('meow')
