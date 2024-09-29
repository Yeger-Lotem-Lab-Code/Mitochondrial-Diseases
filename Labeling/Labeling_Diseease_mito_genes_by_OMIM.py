import pandas as pd

def process_mito_genes_file(mito_genes_path, omim_file_path, labels_output_path, phenotypes_output_path,diseases_output_path, gene_col, ensg_col):
    # Load mito genes file
    mitoGenes_file = pd.read_excel(mito_genes_path)
    gene_names = mitoGenes_file[[ensg_col, gene_col]]
    gene_names.columns = ['Ensembl ID', 'Gene name']

    # Load morbidMap associations
    headers = ['Phenotype', 'Gene Symbols', 'MIM Number', 'Cyto Location']
    omim_file = pd.read_csv(omim_file_path, sep='\t', comment='#', names=headers)
    omim_file_filtered = omim_file[omim_file['Phenotype'].str.contains(r'\(3\)')]

    phenotypes_results = []
    label_results = []
    disease_gene_associations = {}

    for gene in gene_names['Gene name']:
        # Filter the rows where the gene appears in the Gene Symbols column
        gene_rows = omim_file_filtered[omim_file_filtered['Gene Symbols'].str.contains(gene)]
        phenotypes = gene_rows['Phenotype'].tolist()
        label = 1 if phenotypes else 0
        ensg_id = gene_names.loc[gene_names['Gene name'] == gene, 'Ensembl ID'].values[0]
        phenotypes_results.append({'Gene': gene, 'Phenotypes': phenotypes})
        label_results.append({'Gene': gene, 'ENSG ID': ensg_id, 'Label': label})
        # Collect disease-gene associations
        for phenotype in phenotypes:
            if phenotype not in disease_gene_associations:
                disease_gene_associations[phenotype] = []
            disease_gene_associations[phenotype].append(gene)

    disease_gene_df = pd.DataFrame(
        [{'Disease': disease, 'Mito Genes': ', '.join(genes)} for disease, genes in disease_gene_associations.items()]
    )

    # Convert results to DataFrames
    label_results_df = pd.DataFrame(label_results)
    phenotypes_results_df = pd.DataFrame(phenotypes_results)

    # Export to Excel files
    #label_results_df.to_excel(labels_output_path, index=False)
    #phenotypes_results_df.to_excel(phenotypes_output_path, index=False)
    disease_gene_df.to_excel(diseases_output_path, index=False)
# List of files to process with corresponding columns
mito_genes_files_info = [
    {
        'path': r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta\mitoCarta_only_genes.xlsx",
        'gene_col': 'Gene Symbol',
        'ensg_col': 'First ENSG ID'
    },
    #{
        #'path': r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta plus review\merged_mitochondrial_genes_MitoCarta and review.xlsx",
        #'gene_col': 'Gene Symbol',
        #'ensg_col': 'First ENSG ID'
   # },
    #{
        #'path': r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta plus review after manually curation\all_mito_genes_contigs_patches_exc_with_all_names_final.xlsx",
        #'gene_col': 'Gene name',
        #'ensg_col': 'Ensembl ID'
    #}
]

# Output file paths for each input file
labels_output_paths = [
    'Labels_by_OMIM_MitoCarta.xlsx',
    #'Labels_by_OMIM_MitoCarta + Review.xlsx',
    #'Labels_by_OMIM_MitoCarta + Review after manual curation.xlsx'
]

phenotypes_output_paths = [
    'mitochondrial_genes_and_associated_OMIM_phenotypes_MitoCarta.xlsx',
    #'mitochondrial_genes_and_associated_OMIM_phenotypes_MitoCarta + Review.xlsx',
    #'mitochondrial_genes_and_associated_OMIM_phenotypes_MitoCarta + Review after manual curation.xlsx'
]
diseases_output_paths = [
    'Diseases_with_Associated_Mito_Genes_MitoCarta.xlsx',
    #'Diseases_with_Associated_Mito_Genes_MitoCarta + Review.xlsx',
    #'Diseases_with_Associated_Mito_Genes_MitoCarta + Review after manual curation.xlsx'
]


# OMIM file path
omim_file_path = r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\Labeling\morbidmap_Feb2021.txt"

# Loop through each file and process
for file_info, labels_output, phenotypes_output, diseases_output in zip(mito_genes_files_info, labels_output_paths, phenotypes_output_paths, diseases_output_paths):
    process_mito_genes_file(
        file_info['path'],
        omim_file_path,
        labels_output,
        phenotypes_output,
        diseases_output,
        file_info['gene_col'],
        file_info['ensg_col']
    )

print('Processing completed for all files.')