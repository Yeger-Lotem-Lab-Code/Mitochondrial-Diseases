import csv
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np


def mininmal_growth(dismito_genes, nondismito_genes, crispr, output_file):
    dismito = list(dismito_genes)
    nondismito = list(nondismito_genes)
    fieldnames = ['Gene', 'Dis', 'Min']

    print(dismito, 'dismito')
    print(nondismito, 'nondismito')
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for rowNum, row in crispr.iterrows():
            write_row = {}
            write_row['Gene'] = row['Gene']

            # Set mito value based on membership in mito and prot_coding lists
            if row['Gene'] in dismito:
                write_row['Dis'] = 'Yes'
            elif row['Gene'] in nondismito:
                write_row['Dis'] = 'No'
            else:
                write_row['Dis'] = ''

            # If 'mito' column has a value, calculate Min
            if write_row['Dis'] != '':
                val_arr = []
                for item in row.index:
                    if item != 'Gene' and pd.notna(row[item]):
                        val_arr.append(float(row[item]))

                if val_arr:
                    write_row['Min'] = min(val_arr)
                else:
                    write_row['Min'] = None

                writer.writerow(write_row)


mitocarta_genes_and_labels = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\Labeling\Labels_by_OMIM_MitoCarta.xlsx")
dismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 1]['Gene']
dismito_genes = list(dismito_genes)
nondismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 0]['Gene']
nondismito_genes = list(nondismito_genes)
crisprPath = r"C:\Users\liorr\Dropbox\Lior\Characterization\CRISPR\CRISPRGeneEffect.csv"
crispr = pd.read_csv(crisprPath, sep = ',')
# Clean up column names and prepare CRISPR data
crispr.columns = [col.split(' ')[0] for col in crispr.columns]
crispr = crispr.T
crispr.insert(loc=0, column='Gene', value=crispr.index)

output_file = 'DepMap_CRISPR_min_score_final_with_Dismito_26-08-24.csv'

mininmal_growth(dismito_genes, nondismito_genes, crispr, output_file)
