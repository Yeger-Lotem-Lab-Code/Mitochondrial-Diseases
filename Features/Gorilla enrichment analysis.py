import pandas as pd
import numpy as np
from goatools.obo_parser import GODag
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data from TSV files
bio_pro_original = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis\biological_process.xls", sep='\t')
cell_comp_original = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis\cellular_component.xls", sep='\t')
mol_fun_original = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis\molecular_function.xls", sep='\t')

# Create a new column called 'condition' to check the number of genes in each term
bio_pro_original['condition'] = (bio_pro_original['b'] > 46.8) & (bio_pro_original['b'] < 500)
cell_comp_original['condition'] = (cell_comp_original['b'] > 46.8) & (cell_comp_original['b'] < 500)
mol_fun_original['condition'] = (mol_fun_original['b'] > 46.8) & (mol_fun_original['b'] < 500)

# load the manual terms
bio_pro_manual = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis manual choose\gorilla enrichment analysis manual chosen terms.xlsx", sheet_name='bio_pro')
cell_com_manual = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis manual choose\gorilla enrichment analysis manual chosen terms.xlsx", sheet_name='cell_comp')
mol_func_manual = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis manual choose\gorilla enrichment analysis manual chosen terms.xlsx", sheet_name='mol_func')
#filter the original aspects to the manual aspects
bio_pro = bio_pro_original[bio_pro_original['GO Term'].isin(bio_pro_manual['Term'])]
cell_comp = cell_comp_original[cell_comp_original['GO Term'].isin(cell_com_manual['Term'])]
mol_fun = mol_fun_original[mol_fun_original['GO Term'].isin(mol_func_manual['Term'])]

# Load the GO DAG
go_dag = GODag(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis\go-basic.obo")

def get_depth(go_term):
    """Return the depth of a GO term in the hierarchy. If not found, return -1."""
    try:
        term = go_dag[go_term]
        return term.depth
    except KeyError:
        print(f"Warning: GO term {go_term} not found in the GO DAG.")
        return -1
def filter_deepest_terms(df):
    """Filter the dataframe to keep only the deepest terms when there are splits."""
    to_keep = set()
    for term in df.itertuples():
        term_id = term._1
        term_depth = term.depth
        is_deepest = True
        for child in go_dag[term_id].children:
            if child.id in df['GO Term'].values and get_depth(child.id) > term_depth:
                is_deepest = False
                break
        if is_deepest:
            to_keep.add(term_id)
    return df[df['GO Term'].isin(to_keep)]

'''
# Filter, add depth information, and sort by depth and FDR q-value
def process_go_terms(df):
    filtered_terms = df[df['condition']].copy()
    filtered_terms['depth'] = filtered_terms['GO Term'].apply(get_depth)
    filtered_terms = filtered_terms[filtered_terms['depth'] != -1]  # Remove terms with depth -1
    filtered_terms = filtered_terms.sort_values(by=['depth', 'FDR q-value'], ascending=[False, True])
    filtered_terms = filter_deepest_terms(filtered_terms)
    return filtered_terms


# Process each GO aspect
bio_pro_filtered = process_go_terms(bio_pro)
cell_comp_filtered = process_go_terms(cell_comp)
mol_fun_filtered = process_go_terms(mol_fun)

# Get top 10 GO terms for each aspect
top_bio_pro = bio_pro_filtered.head(10)
top_cell_comp = cell_comp_filtered.head(10)
top_mol_fun = mol_fun_filtered.head(10)
'''
# Extract unique genes from the top GO terms
genes = set()
for df in [bio_pro, cell_comp, mol_fun]:
    for gene_list in df['Genes']:
        gene_list = gene_list.replace('[', '').replace(']', '')
        gene_list = gene_list.split(', ')
        for gene in gene_list:
            gene_name = gene.split(' - ', 1)[0].replace(' ', '')
            if gene_name and gene_name[0].isupper() and not gene_name[0].isdigit():
                genes.add(gene_name)

# Create an empty dataframe with genes as rows and top GO terms as columns
# Create an empty dataframe with genes as rows and top GO terms as columns
columns_bio_pro = [f"Biological Process: {row._1} ({row.Description})" for row in bio_pro.itertuples()]
columns_cell_comp = [f"Cellular Component: {row._1} ({row.Description})" for row in cell_comp.itertuples()]
columns_mol_fun = [f"Molecular Function: {row._1} ({row.Description})" for row in mol_fun.itertuples()]
df_genes = pd.DataFrame(index=sorted(genes), columns=columns_bio_pro + columns_cell_comp + columns_mol_fun).fillna(0)

# Populate the dataframe
# Populate the dataframe
def populate_df(df, aspect_prefix):
    for term in df.itertuples():
        go_term = f"{aspect_prefix}: {term._1} ({term.Description})"
        gene_list = term.Genes.replace('[', '').replace(']', '').split(', ')
        for gene in gene_list:
            gene_name = gene.split(' - ', 1)[0].replace(' ', '')
            if gene_name and gene_name[0].isupper() and not gene_name[0].isdigit():
                df_genes.at[gene_name, go_term] = 1

populate_df(bio_pro, 'Biological Process')
populate_df(cell_comp, 'Cellular Component')
populate_df(mol_fun, 'Molecular Function')
df_genes.to_csv(r'C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis manual choose\Gorilla enrichment analysis manual choose.csv')
# Print the dataframe
print(df_genes)
# Plotting number of genes in each term for each aspect
def plot_go_terms(df, aspect, title):
    df['num_genes'] = df['Genes'].apply(lambda genes: len(genes.split(', ')))
    plt.figure(figsize=(12, 8))
    sns.barplot(data=df, x='Description', y='num_genes', palette='viridis')
    plt.xlabel('GO Terms', fontsize=20)
    plt.ylabel('Number of Genes', fontsize=20)
    plt.title(f'Number of Genes in Top 10 {aspect} GO Terms', fontsize=20)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.show()

plot_go_terms(bio_pro, 'Biological Process', 'Top 10 Biological Process GO Terms')
plot_go_terms(cell_comp, 'Cellular Component', 'Top 10 Cellular Component GO Terms')
plot_go_terms(mol_fun, 'Molecular Function', 'Top 10 Molecular Function GO Terms')