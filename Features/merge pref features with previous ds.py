import pandas as pd
import numpy as np

# Load datasets
prev_ds = pd.read_excel(
    r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset.xlsx")
pref_num_elevated_interactors = pd.read_excel(
    r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_elevated_interactors_mito.xlsx")
pref_num_interactors = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_interactors_mito.xlsx")
pref_num_specific_interactions = pd.read_excel(
    r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_specific_interactions_mito.xlsx")
pref_expression = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_expression_mito_genes.xlsx")
pref_development_tissue = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoDevelopment_preferential\development_preferential_values_tissue.csv")
pref_development_time = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoDevelopment_preferential\development_preferential_values_time.csv")

# Reshape development preferential data into wide format with desired column names
pref_development_tissue['New_Column'] = pref_development_tissue.apply(lambda x: f"tissue_{x['Tissue']}_{x['Time_Point']}_pref", axis=1)
pref_development_tissue_wide = pref_development_tissue.pivot(index='Gene', columns='New_Column', values='Preferential_Value_Tissue').reset_index()

pref_development_time['New_Column'] = pref_development_time.apply(lambda x: f"time_{x['Tissue']}_{x['Time_Point']}_pref", axis=1)
pref_development_time_wide = pref_development_time.pivot(index='Gene', columns='New_Column', values='Preferential_Value_Time').reset_index()


# Define a function to merge and drop specific columns
def merge_and_drop(prev_ds, new_df, merge_key='Gene Ensembl ID', substring=None):
    if substring:
        drop_cols = [col for col in prev_ds.columns if substring in col]
    else:
        drop_cols = []

    # Perform the merge
    if substring == '_expression':
        merged_ds = pd.merge(prev_ds, new_df, left_on=merge_key, right_on='ENSG', how='left')
        merged_ds.drop(columns=['ENSG'] + drop_cols, inplace=True)
    else:
        merged_ds = pd.merge(prev_ds, new_df, left_on=merge_key, right_on='Gene', how='left')
        merged_ds.drop(columns=['Gene'] + drop_cols, inplace=True)

    return merged_ds

substrings = ['num_elevated_interactors', 'num_interactors', 'num_specific_interactions', '_expression']

# Merge the datasets and drop previous columns
merged_ds = merge_and_drop(prev_ds, pref_num_elevated_interactors, merge_key='Gene Ensembl ID', substring=substrings[0].replace('_', ''))
merged_ds = merge_and_drop(merged_ds, pref_num_interactors, merge_key='Gene Ensembl ID', substring=substrings[1].replace('_', ''))
merged_ds = merge_and_drop(merged_ds, pref_num_specific_interactions, merge_key='Gene Ensembl ID', substring=substrings[2].replace('_', ''))
merged_ds = merge_and_drop(merged_ds, pref_expression, merge_key='Gene Ensembl ID', substring=substrings[3])

# Merge development preferential values
merged_ds = merge_and_drop(merged_ds, pref_development_tissue_wide, merge_key='Gene Ensembl ID', substring='develop')
merged_ds = merge_and_drop(merged_ds, pref_development_time_wide, merge_key='Gene Ensembl ID', substring=None)

#save the merged dataset
merged_ds.to_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_pref.csv", index=False)


# Print the resulting dataframe
print(merged_ds)