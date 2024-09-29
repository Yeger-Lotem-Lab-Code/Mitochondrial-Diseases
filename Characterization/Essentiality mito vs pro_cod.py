import csv
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np

def mininmal_growth(mitoGenes, protein_coding_genes, crispr, output_file):
    mito = list(mitoGenes)
    prot_coding = list(protein_coding_genes)
    fieldnames = ['Gene', 'mito', 'Min']

    print(mito, 'mito')
    print(prot_coding, 'prot')
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for rowNum, row in crispr.iterrows():
            write_row = {}
            write_row['Gene'] = row['Gene']

            # Set mito value based on membership in mito and prot_coding lists
            if row['Gene'] in mito:
                write_row['mito'] = 'Yes'
            elif row['Gene'] in prot_coding:
                write_row['mito'] = 'No'
            else:
                write_row['mito'] = ''

            # If 'mito' column has a value, calculate Min
            if write_row['mito'] != '':
                val_arr = []
                for item in row.index:
                    if item != 'Gene' and pd.notna(row[item]):
                        val_arr.append(float(row[item]))

                if val_arr:
                    write_row['Min'] = min(val_arr)
                else:
                    write_row['Min'] = None

                writer.writerow(write_row)


mitoGenesPath = r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Human.MitoCarta3.0.xls"
mitoGenes = pd.read_excel(mitoGenesPath, 'A Human MitoCarta3.0')['Symbol']
mitoGenes = list(mitoGenes)
protein_coding_genes = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Characterization\PROTEIN CODING GENES.txt", sep=',')['Gene name']
protein_coding_genes = list(protein_coding_genes)

mitoGenes_set = set(mitoGenes)
protein_coding_genes_set = set(protein_coding_genes)
filtered_protein_coding_genes = protein_coding_genes_set - mitoGenes_set
filtered_protein_coding_genes = list(filtered_protein_coding_genes)
filtered_protein_coding_genes = [gene for gene in filtered_protein_coding_genes if gene == gene]

crisprPath = r"C:\Users\liorr\Dropbox\Lior\Characterization\CRISPR\CRISPRGeneEffect.csv"
crispr = pd.read_csv(crisprPath, sep = ',')
# Clean up column names and prepare CRISPR data
crispr.columns = [col.split(' ')[0] for col in crispr.columns]
crispr = crispr.T
crispr.insert(loc=0, column='Gene', value=crispr.index)

output_file = 'DepMap_CRISPR_min_score_final_with_mito_12-2-23.csv'

mininmal_growth(mitoGenes, filtered_protein_coding_genes, crispr, output_file)
