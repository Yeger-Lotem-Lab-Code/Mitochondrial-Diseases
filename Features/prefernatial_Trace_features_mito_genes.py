import pandas as pd
import numpy as np
import statistics
import re

def clean_column_names(df):
    df.columns = [re.sub(r'\W+', '', col) for col in df.columns]
    return df

def filter_genes(data, gene_list):
    # Filter the median expression data based on the provided list of genes
    return data[data['Gene_ID'].isin(gene_list)]

def extract_tissue_features(trace_file, feature_suffix):
    feature_data = trace_file.filter(regex=f'Gene_ID|.*_{feature_suffix}$')
    feature_data.columns = [col.split('_')[0] for col in feature_data.columns]
    return feature_data

def add_whole_tissues_columns(feature, tissue_associations):
    # Clean the column names in median_expression
    feature = clean_column_names(feature)
    previous_columns_names = feature.columns[1:]

    for tissue in tissue_associations['mitoModel'].unique():
        associated_tissues = tissue_associations.loc[tissue_associations['mitoModel'] == tissue, 'expression|num_interactors|num_elevated_interactors|num_specific_interactions'].tolist()
        associated_tissues = [re.sub(r'\W+', '', col).lower() for col in associated_tissues]
        if len(associated_tissues) == 1:
            feature[tissue] = feature[associated_tissues[0]]
        else:
            feature[tissue] = feature[associated_tissues].median(axis=1)
    feature = feature.drop(columns= previous_columns_names)
    return feature

def calculate_median_and_iqr(feature):
    gene_statistics = []

    for _, row in feature.iterrows():
        temp_list = row[1:].astype(float).tolist()
        # Assuming the first column is 'ENSG'
        median_all_tissues = statistics.median(temp_list)
        Q1 = np.quantile(temp_list, 0.25)
        Q3 = np.quantile(temp_list, 0.75)
        IQR = Q3 - Q1
        if IQR == 0:
            IQR = 1
        gene_statistics.append({
            'ENSG': row['Gene'],
            'Median_all_tissues': median_all_tissues,
            'IQR': IQR
        })
    return pd.DataFrame(gene_statistics)

def calculate_preferential_expression(feature, feature_with_whole_columns, gene_statistics, feature_name):
    # Merge the filtered data with gene_statistics on 'ENSG'
    merged_data = pd.merge(feature_with_whole_columns, gene_statistics, left_on='Gene', right_on='ENSG')
    merged_data.drop(columns=['ENSG'], inplace=True)
    # Calculate preferential expression scores
    for tissue in feature_with_whole_columns.columns[1:]:
        if tissue not in merged_data.columns:
            continue
        merged_data[tissue + '_preferential_' + feature_name] = (merged_data[tissue] - merged_data['Median_all_tissues']) / merged_data['IQR']

    # Select relevant columns for output
    output_columns = ['Gene'] + [col for col in merged_data.columns if col.endswith('_preferential_' + feature_name)]
    merged_data[output_columns].to_excel(f'preferential_{feature_name}_mito.xlsx', index=False)


trace_file = pd.read_csv(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\2_TRACE_Dataset_united_development_features.csv")
gene_list = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta\Gene_List_1.xlsx")['Gene Ensembl ID']

# Filter data for the specified genes
filtered_data = filter_genes(trace_file, gene_list)

tissue_associations = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx", 'mitoTrace')
tissue_associations['expression|num_interactors|num_elevated_interactors|num_specific_interactions'] = tissue_associations['expression|num_interactors|num_elevated_interactors|num_specific_interactions'].apply(lambda x: re.sub(r'\W+', '', x))

for feature_suffix in ['num_interactors', 'num_elevated_interactors', 'num_specific_interactions']:
    feature = extract_tissue_features(filtered_data, feature_suffix)
    feature_with_whole_columns = add_whole_tissues_columns(feature, tissue_associations)
    gene_statistics = calculate_median_and_iqr(feature_with_whole_columns)
    calculate_preferential_expression(feature, feature_with_whole_columns, gene_statistics, feature_suffix)