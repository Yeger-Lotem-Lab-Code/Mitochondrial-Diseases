import pandas as pd
import numpy as np
import statistics
import re


def load_data(median_expression_file, tissue_associations_file):
    # Load the median expression data
    median_expression = pd.read_csv(median_expression_file, delimiter=',')
    # Remove trailing period and number from ENSG IDs
    median_expression['ENSG'] = median_expression['ENSG'].str.replace(r'\.\d+$', '', regex=True)

    # Load tissue associations data
    tissue_associations = pd.read_excel(tissue_associations_file, 'mitoTrace')

    return median_expression, tissue_associations

def clean_column_names(df):
    df.columns = [re.sub(r'\W+', '', col) for col in df.columns]
    return df


def filter_genes(data, gene_list):
    # Filter the median expression data based on the provided list of genes
    return data[data['ENSG'].isin(gene_list)]


def add_whole_tissues_columns(median_expression, tissue_associations):
    # Clean the column names in median_expression
    median_expression = clean_column_names(median_expression)
    previous_columns_names = median_expression.columns[1:]

    for tissue in tissue_associations['mitoModel'].unique():
        associated_tissues = tissue_associations.loc[tissue_associations['mitoModel'] == tissue, 'expression|num_interactors|num_elevated_interactors|num_specific_interactions'].tolist()
        associated_tissues = [re.sub(r'\W+', '', col) for col in associated_tissues]
        if len(associated_tissues) == 1:
            median_expression[tissue] = median_expression[associated_tissues[0]]
        else:
            median_expression[tissue] = median_expression[associated_tissues].median(axis=1)
    tissues_in_both = set(median_expression.columns[1:]).intersection(tissue_associations['mitoModel'].unique())
    previous_columns_names_without_tissues_in_both = [col for col in previous_columns_names if col not in tissues_in_both]
    median_expression = median_expression.drop(columns= previous_columns_names_without_tissues_in_both)
    return median_expression


def calculate_median_and_iqr(median_expression):
    gene_statistics = []

    for _, row in median_expression.iterrows():
        temp_list = row[1:].astype(float).tolist()  # Assuming the first two columns are 'ENSG' and 'Name'
        median_all_tissues = statistics.median(temp_list)
        Q1 = np.quantile(temp_list, 0.25)
        Q3 = np.quantile(temp_list, 0.75)
        IQR = Q3 - Q1
        if IQR == 0:
            IQR = 1
        gene_statistics.append({
            'ENSG': row['ENSG'],
            'Median_all_tissues': median_all_tissues,
            'IQR': IQR
        })

    return pd.DataFrame(gene_statistics)


def calculate_preferential_expression(median_expression_with_whole, gene_statistics,output_file):
    # Merge the filtered data with gene_statistics on 'ENSG' and 'Name'
    merged_data = pd.merge(median_expression_with_whole, gene_statistics, on=['ENSG'])

    # Calculate preferential expression scores
    for tissue in median_expression_with_whole.columns[1:]:
        if tissue not in merged_data.columns:
            continue
        merged_data[tissue + '_preferential_expression'] = (merged_data[tissue] - merged_data['Median_all_tissues']) / merged_data['IQR']

    # Select relevant columns for output
    output_columns = ['ENSG'] + [col for col in merged_data.columns if col.endswith('_preferential_expression')]
    merged_data[output_columns].to_excel(output_file, index=False)




if __name__ == "__main__":
    # Define file paths
    median_expression_file = r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\GTEX_median_expression.csv"
    med_exp_file = r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\GTEX_median_expression.csv"
    tissue_associations_file = r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx"
    preferential_scores_file = 'preferential_expression_mito_genes.xlsx'

    # Load data
    median_expression, tissue_associations = load_data(median_expression_file, tissue_associations_file)
    tissue_associations['expression|num_interactors|num_elevated_interactors|num_specific_interactions'] = tissue_associations['expression|num_interactors|num_elevated_interactors|num_specific_interactions'].apply(lambda x: re.sub(r'\W+', '', x))

    # Define the list of genes to process
    gene_list = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta\Gene_List_1.xlsx")['Gene Ensembl ID']  # Replace with your actual list of ENSG IDs

    # Filter data for the specified genes
    filtered_data = filter_genes(median_expression, gene_list)

    # Add "whole" tissues columns
    median_expression_with_whole = add_whole_tissues_columns(filtered_data, tissue_associations)

    # Calculate median expression and IQR for each gene
    gene_statistics = calculate_median_and_iqr(median_expression_with_whole)

    # Calculate preferential expression
    calculate_preferential_expression(median_expression_with_whole, gene_statistics, preferential_scores_file)