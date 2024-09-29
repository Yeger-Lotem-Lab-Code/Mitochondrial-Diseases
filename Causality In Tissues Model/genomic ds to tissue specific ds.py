import pandas as pd
import numpy as np


def create_variables(mitoGenes, mitoModel_tissues):
    variables_df = pd.DataFrame(columns=['Gene Name', 'Ensembl_ID', 'Tissue'])
    mitoGenes_series = pd.Series(mitoGenes['Gene_name'].unique())
    mitoEnsemblIDs_series = pd.Series(mitoGenes['Ensembl_ID'].unique())
    genes_names_expanded = mitoGenes_series.repeat(len(mitoModel_tissues))
    genes_ensembl_expanded = mitoEnsemblIDs_series.repeat(len(mitoModel_tissues))
    variables_df['Gene Name'] = genes_names_expanded.reset_index(drop=True)
    variables_df['Ensembl_ID'] = genes_ensembl_expanded.reset_index(drop=True)
    variables_df['Tissue'] = pd.concat([mitoModel_tissues] * len(mitoGenes), ignore_index=True)
    return variables_df


def assign_feature_values(labeled_variables_df, genomic_ds):
    # Iterate over each feature column in the genomic dataset, skipping the first 3 columns
    for column in genomic_ds.columns[3:]:
        # Check if the column contains any non-tissue-specific feature name
        if any(feature in column for feature in non_tissue_features):
            continue  # Skip this column, as it is not tissue-specific

        # Extract tissue and feature names from the column name
        if '_' in column:
            tissue, feature = column.split('_', 1)  # Split only on the first '_'
            tissue = tissue.lower().replace(' ', '').replace('-', '')

            # Check if the tissue from labeled_variables_df appears in the column name
            if feature not in labeled_variables_df.columns:
                labeled_variables_df[feature] = np.nan
            matching_genes = labeled_variables_df['Gene Name'].isin(genomic_ds['Gene_x']) & \
                                 (labeled_variables_df['Tissue'] == tissue)
            labeled_variables_df.loc[matching_genes, feature] = labeled_variables_df.loc[matching_genes, 'Gene Name'].map(
                    genomic_ds.set_index('Gene_x')[column]
                )
    return labeled_variables_df

def assign_development_feature_values(labeled_variables_df, genomic_ds):
    relevantcolumns = ['Gene_x'] + [col for col in genomic_ds.columns if 'prefexp' in col]
    for column in relevantcolumns[1:]:
        type,tissue,time = column.split('_', 2)
        tissue = tissue.lower().replace(' ', '').replace('-', '')
        feature = type+'_'+time
        if feature not in labeled_variables_df.columns:
            labeled_variables_df[feature] = np.nan
        matching_genes = labeled_variables_df['Gene Name'].isin(genomic_ds['Gene_x']) & \
                         (labeled_variables_df['Tissue'] == tissue)
        labeled_variables_df.loc[matching_genes, feature] = labeled_variables_df.loc[matching_genes, 'Gene Name'].map(
            genomic_ds.set_index('Gene_x')[column]
        )
    return labeled_variables_df

def assign_genomic_feature_values(labeled_variables_df, genomic_ds):
    relevantcolumns = [col for col in genomic_ds.columns if any(feat in col for feat in non_tissue_features[1:])]
    matching_genes = labeled_variables_df['Gene Name'].isin(genomic_ds['Gene_x'])
    for column in relevantcolumns[1:]:
        labeled_variables_df.loc[matching_genes, column] = labeled_variables_df.loc[matching_genes, 'Gene Name'].map(
            genomic_ds.set_index('Gene_x')[column]
        )
    return labeled_variables_df





most_confidence_mitoGenes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confident dataset for mitocarta omim positive genes\Most_confidence_genes_manual_labels.xlsx")
mitoModel_tissues = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel_tissues.xlsx")['mitoModel']
mitoModel_tissues = mitoModel_tissues.str.replace(' ', '').str.lower().str.replace('-', '')
most_confidence_variables_df = create_variables(most_confidence_mitoGenes, mitoModel_tissues)
most_confident_labels_unique = most_confidence_mitoGenes.drop_duplicates(subset=['Gene_name', 'Candidate_tissue'])
most_confident_labels_unique['manual label'] = 1
labeled_variables_df = most_confidence_variables_df.merge(most_confident_labels_unique, left_on=['Gene Name', 'Tissue'], right_on=['Gene_name', 'Candidate_tissue'], how='left')
labeled_variables_df['label'] = labeled_variables_df['manual label'].fillna(0).astype(int)
# Add the 'label' column, assigning 1 where there's a match, otherwise 0
labeled_variables_df = labeled_variables_df.drop(['Gene_name', 'Candidate_tissue', 'Ensembl_ID_y', 'manual label'], axis=1)

genomic_ds = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla.csv")
non_tissue_features = ['prefexp','GO:', 'ProteinLength', 'TargetP_Score', 'PGC_Induction_Score',
                       'CoexpressionGnfN50_Score', 'Yeast', 'Rickettsia', 'hg19',
                       'mitoPathway', 'protein enrichment', 'rev_cat','core-variable']
labeled_variables_df = assign_feature_values(labeled_variables_df, genomic_ds)
labeled_variables_df = assign_development_feature_values(labeled_variables_df, genomic_ds)
labeled_variables_df = assign_genomic_feature_values(labeled_variables_df, genomic_ds)
labeled_variables_df.to_csv(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\Most_confidence_genes_manual_labels_with_features.csv", index=False)

