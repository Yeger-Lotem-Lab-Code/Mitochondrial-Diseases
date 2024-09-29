import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import obonet
import numpy as np

def plot_histogram(df, aspect):
    terms_counts = df['GO term name'].value_counts()[1:]  # Remove the first term which is the most common term
    top20_terms = terms_counts.head(10)  # Top 20 terms

    plt.figure(figsize=(16, 14))  # Increase figure size
    ax = top20_terms.plot.bar()
    ax.set_xlabel('Term', fontsize=20)  # Increase font size of x-axis label
    ax.set_ylabel('Number of Associated genes', fontsize=20, weight='bold')
    ax.set_title('Top 10 ' +aspect+' GO terms (Descending Order)', fontsize=20, weight='bold')

    for label in ax.get_xticklabels():
        label.set_fontsize(20)
        label.set_fontweight('bold')
        label.set_rotation(45)
        label.set_ha('right')

    for label in ax.get_yticklabels():
        label.set_fontsize(20)
        label.set_fontweight('bold')
    plt.tight_layout()
    plt.show()

    return terms_counts.head(10)

def create_feature_dataset(aspect_df, top10_terms):
    feature_dataset = pd.DataFrame(index=aspect_df['Gene stable ID'].unique(), columns=top10_terms.index.tolist())
    feature_dataset = feature_dataset.fillna(0)
    for index, row in aspect_df.iterrows():
        gene = row['Gene stable ID']
        term = row['GO term name']
        if term in top10_terms:
            feature_dataset.loc[gene, term] = 1
    return feature_dataset

# Main
ontology_files = [
    r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\GO terms from Biomart\mart_export (9).txt"
]
ontology_file = pd.read_csv(ontology_files[0], sep=',')
aspect_names = ['cellular_component', 'biological_process', 'molecular_function','mitochondrion_localization']
#url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
#graph = obonet.read_obo(url)
#id_to_aspect = {id_: data.get('namespace') for id_, data in graph.nodes(data=True)}
#id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}

for aspect in aspect_names:
    if aspect=='mitochondrion_localization':
        mito_genes_ontology_aspect = ontology_file[ontology_file['GO domain'] == 'cellular_component']
        mito_genes_ontology_aspect = mito_genes_ontology_aspect[mito_genes_ontology_aspect['GO term name'].str.contains('mitoc', case=False, na=False)]
    elif aspect=='cellular_component':
        mito_genes_ontology_aspect = ontology_file[ontology_file['GO domain'] == aspect]
        mito_genes_ontology_aspect = mito_genes_ontology_aspect[
            ~mito_genes_ontology_aspect['GO term name'].str.contains('mitoc|organelle', case=False, na=False)]
    else:
        mito_genes_ontology_aspect = ontology_file[ontology_file['GO domain'] == aspect]
    mito_genes_ontology_aspect = mito_genes_ontology_aspect[['Gene stable ID', 'GO term name']]
    top10_terms = plot_histogram(mito_genes_ontology_aspect, aspect)
    feature_dataset = create_feature_dataset(mito_genes_ontology_aspect, top10_terms)

    #feature_dataset.to_excel(f'{aspect}_top10_binary.xlsx', index=True)