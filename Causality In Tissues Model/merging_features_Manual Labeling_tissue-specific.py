import pandas as pd

# create entire genes-tissues variables
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

def reshape_table(pref_table, gene_col, value_name):
    pref_table_long = pd.melt(pref_table, id_vars=[gene_col], var_name='Tissue_Feature', value_name=value_name)
    pref_table_long[['Tissue', 'Feature']] = pref_table_long['Tissue_Feature'].str.split('_', n=1, expand=True)
    pref_table_long.drop(columns=['Tissue_Feature', 'Feature'], inplace=True)
    pref_table_long['Tissue'] = pref_table_long['Tissue'].str.lower().str.replace(' ', '')
    return pref_table_long

def merge_tissue_specific_features(variables_df, mitoTrace, pref_exp_mito, pref_num_ele_interactors, pref_num_interactors, pref_num_specific_interactions, mito_proact_scores, mitoAge, mitoTSScore, mitoPeakIntensity, mito_development_pref_time, mito_development_pref_tissue):
    variables_df['Tissue'] = variables_df['Tissue'].str.lower().str.replace(' ', '')
    mitoTrace_filtered = mitoTrace.loc[:, ~mitoTrace.columns.str.contains('expression|interactors|interactions|development|mediantipapathways', case=False)]
    final_tissue_specific_df = variables_df.merge(mitoTrace_filtered, left_on=['Ensembl_ID_x', 'Tissue'],
                                                  right_on=['Gene', 'Tissue'], how='left')
    final_tissue_specific_df = final_tissue_specific_df.drop(['Unnamed: 0', 'Gene'], axis=1)
    #change the strings of muscle-skeletal in the tissue column to muscleskeletal in pref_exp_mito, pref_num_ele_interactors, pref_num_interactors, pref_num_specific_interactions

    pref_exp_mito_long = reshape_table(pref_exp_mito, 'ENSG', 'Preferential_Expression')
    pref_exp_mito_long['Tissue'] = pref_exp_mito_long['Tissue'].str.lower().str.replace('-', '')
    pref_num_ele_interactors_long = reshape_table(pref_num_ele_interactors, 'Gene', 'Num_Elevated_Interactors')
    pref_num_ele_interactors_long['Tissue'] = pref_num_ele_interactors_long['Tissue'].str.lower().str.replace('-', '')
    pref_num_interactors_long = reshape_table(pref_num_interactors, 'Gene', 'Num_Interactors')
    pref_num_interactors_long['Tissue'] = pref_num_interactors_long['Tissue'].str.lower().str.replace('-', '')
    pref_num_specific_interactions_long = reshape_table(pref_num_specific_interactions, 'Gene', 'Num_Specific_Interactions')
    pref_num_specific_interactions_long['Tissue'] = pref_num_specific_interactions_long['Tissue'].str.lower().str.replace('-', '')

    final_tissue_specific_df = final_tissue_specific_df.merge(pref_exp_mito_long, left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['ENSG', 'Tissue'], how='left')
    final_tissue_specific_df.drop(columns=['ENSG'], inplace=True)

    final_tissue_specific_df = final_tissue_specific_df.merge(pref_num_ele_interactors_long,
                                                              left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['Gene', 'Tissue'], how='left')
    final_tissue_specific_df.drop(columns=['Gene'], inplace=True)

    final_tissue_specific_df = final_tissue_specific_df.merge(pref_num_interactors_long,
                                                              left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['Gene', 'Tissue'], how='left')
    final_tissue_specific_df.drop(columns=['Gene'], inplace=True)

    final_tissue_specific_df = final_tissue_specific_df.merge(pref_num_specific_interactions_long,
                                                              left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['Gene', 'Tissue'], how='left')
    final_tissue_specific_df.drop(columns=['Gene'], inplace=True)

    proact_scores_long = reshape_table(mito_proact_scores, 'Unnamed: 0', 'mean proact score')
    proact_scores_long['Tissue'] = proact_scores_long['Tissue'].str.lower().str.replace('-', '')
    final_tissue_specific_df = final_tissue_specific_df.merge(proact_scores_long, left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['Unnamed: 0', 'Tissue'], how='left')
    final_tissue_specific_df.drop(columns=['Unnamed: 0'], inplace=True)

    mitoAge['Tissue'] = mitoAge['Tissue'].str.lower().str.replace(' ', '').str.replace('-', '')
    final_tissue_specific_df = final_tissue_specific_df.merge(mitoAge, left_on=['Gene Name', 'Tissue'],
                                                              right_on=['index', 'Tissue'], how='left')
    final_tissue_specific_df = final_tissue_specific_df.drop(['Unnamed: 0', 'index'], axis=1)
    mitoTSScore['Tissue'] = mitoTSScore['Tissue'].str.lower().str.replace(' ', '').str.replace('-', '')
    final_tissue_specific_df = final_tissue_specific_df.merge(mitoTSScore, left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['index', 'Tissue'], how='left')
    final_tissue_specific_df = final_tissue_specific_df.drop(['index'], axis=1)
    mitoPeakIntensity['Tissue'] = mitoPeakIntensity['Tissue'].str.lower().str.replace(' ', '').str.replace('-', '')
    final_tissue_specific_df = final_tissue_specific_df.merge(mitoPeakIntensity, left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['Gene', 'Tissue'], how='left')
    final_tissue_specific_df = final_tissue_specific_df.drop(['Gene'], axis=1)

    # Reshape the development preferential tables
    mito_development_pref_time_long = pd.melt(mito_development_pref_time, id_vars=['Gene'], var_name='Tissue_Time',
                                              value_name='Preferential_Value_Time')
    mito_development_pref_time_long[['Type', 'Tissue', 'Time_Point', 'Suffix']] = mito_development_pref_time_long[
        'Tissue_Time'].str.split('_', expand=True)
    mito_development_pref_time_long = mito_development_pref_time_long.drop(columns=['Type', 'Suffix'])
    mito_development_pref_time_long['Tissue'] = mito_development_pref_time_long['Tissue'].str.lower().str.replace(' ',
                                                                                                                  '')

    mito_development_pref_tissue_long = pd.melt(mito_development_pref_tissue, id_vars=['Gene'], var_name='Tissue_Time',
                                                value_name='Preferential_Value_Tissue')
    mito_development_pref_tissue_long[['Type', 'Tissue', 'Time_Point', 'Suffix']] = mito_development_pref_tissue_long[
        'Tissue_Time'].str.split('_', expand=True)
    mito_development_pref_tissue_long = mito_development_pref_tissue_long.drop(columns=['Type', 'Suffix'])
    mito_development_pref_tissue_long['Tissue'] = mito_development_pref_tissue_long['Tissue'].str.lower().str.replace(
        ' ', '')

    # Pivot the development preferential values to wide format
    mito_development_pref_time_wide = mito_development_pref_time_long.pivot_table(index=['Gene', 'Tissue'],
                                                                                  columns='Time_Point',
                                                                                  values='Preferential_Value_Time',
                                                                                  aggfunc='first').reset_index()
    mito_development_pref_time_wide.columns = [f'Time_{col}' if col not in ['Gene', 'Tissue'] else col for col in
                                               mito_development_pref_time_wide.columns]

    mito_development_pref_tissue_wide = mito_development_pref_tissue_long.pivot_table(index=['Gene', 'Tissue'],
                                                                                      columns='Time_Point',
                                                                                      values='Preferential_Value_Tissue',
                                                                                      aggfunc='first').reset_index()
    mito_development_pref_tissue_wide.columns = [f'Tissue_{col}' if col not in ['Gene', 'Tissue'] else col for col in
                                                 mito_development_pref_tissue_wide.columns]

    final_tissue_specific_df = final_tissue_specific_df.merge(mito_development_pref_time_wide,
                                                              left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['Gene', 'Tissue'], how='left')
    final_tissue_specific_df.drop(columns=['Gene'], inplace=True)

    final_tissue_specific_df = final_tissue_specific_df.merge(mito_development_pref_tissue_wide,
                                                              left_on=['Ensembl_ID_x', 'Tissue'],
                                                              right_on=['Gene', 'Tissue'], how='left')
    final_tissue_specific_df.drop(columns=['Gene'], inplace=True)

    return final_tissue_specific_df

def merge_genomic_features(final_tissue_specific_df, mitoCarta_dataset, mitoCoreVar, mitoTSCategories,
                           mito_disCategories,  GO_terms_enriched, mitoPathways):
    mitoCarta_dataset['Tissue'] = mitoCarta_dataset['Tissue'].str.lower().str.replace(' ', '')
    final_whole_dataset = final_tissue_specific_df.merge(mitoCarta_dataset, left_on=['Ensembl_ID_x', 'Tissue'],
                                                         right_on=['Gene', 'Tissue'], how='left')
    final_whole_dataset = final_whole_dataset.drop(['Unnamed: 0', 'Gene'], axis=1)

    mitoCoreVar['Tissue'] = mitoCoreVar['Tissue'].str.lower().str.replace(' ', '')
    final_whole_dataset = final_whole_dataset.merge(mitoCoreVar, left_on=['Ensembl_ID_x', 'Tissue'],
                                                    right_on=['Gene', 'Tissue'], how='left')
    final_whole_dataset = final_whole_dataset.drop(['Gene'], axis=1)

    mitoTSCategories['Tissue'] = mitoTSCategories['Tissue'].str.lower().str.replace(' ', '')
    final_whole_dataset = final_whole_dataset.merge(mitoTSCategories, left_on=['Ensembl_ID_x', 'Tissue'],
                                                    right_on=['Gene', 'Tissue'], how='left')
    final_whole_dataset = final_whole_dataset.drop(['Gene'], axis=1)

    final_whole_dataset = final_whole_dataset.merge(mito_disCategories, left_on=['Ensembl_ID_x'],
                                                    right_on=['Unnamed: 0'], how='left')
    final_whole_dataset = final_whole_dataset.drop(['Unnamed: 0'], axis=1)

    #final_whole_dataset = final_whole_dataset.merge(mitoLocaliztion, left_on=['Ensembl_ID_x'],
     #                                               right_on=['Unnamed: 0'], how='left')
    #final_whole_dataset = final_whole_dataset.drop(['Unnamed: 0'], axis=1)

    final_whole_dataset = final_whole_dataset.merge(GO_terms_enriched, left_on=['Ensembl_ID_x'],
                                                    right_on=['ENSG'], how='left')
    final_whole_dataset = final_whole_dataset.drop(['ENSG', 'Gene'], axis=1)

    final_whole_dataset = final_whole_dataset.merge(mitoPathways, left_on=['Gene Name'],
                                                    right_on=['Gene.1'], how='left')
    final_whole_dataset = final_whole_dataset.drop(['Gene', 'Gene.1'], axis=1)
    return final_whole_dataset

# Load the datasets
most_confidence_mitoGenes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confident dataset for mitocarta omim positive genes\Most_confidence_genes_manual_labels.xlsx")
mitoModel_tissues = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel_tissues.xlsx")['mitoModel']
mitoModel_tissues = mitoModel_tissues.str.replace(' ', '').str.lower().str.replace('-', '')
most_confidence_variables_df = create_variables(most_confidence_mitoGenes, mitoModel_tissues)

# Add the manual labels
most_confident_labels = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confident dataset for mitocarta omim positive genes\Most_confidence_genes_manual_labels.xlsx")
most_confident_labels['manual label'] = 1
most_confident_labels_unique = most_confident_labels.drop_duplicates(subset=['Gene_name', 'Candidate_tissue'])
labeled_variables_df = most_confidence_variables_df.merge(most_confident_labels_unique, left_on=['Gene Name', 'Tissue'], right_on=['Gene_name', 'Candidate_tissue'], how='left')
labeled_variables_df['label'] = labeled_variables_df['manual label'].fillna(0).astype(int)

# Add the 'label' column, assigning 1 where there's a match, otherwise 0
labeled_variables_df = labeled_variables_df.drop(['Gene_name', 'Candidate_tissue', 'Ensembl_ID_y', 'manual label'], axis=1)

# Import tissue-specific features
mitoTrace = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_dataset.xlsx")
pref_exp_mito = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_expression_mito_genes.xlsx")
pref_num_ele_interactors = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_elevated_interactors_mito.xlsx")
pref_num_interactors = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_interactors_mito.xlsx")
pref_num_specific_interactions = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_specific_interactions_mito.xlsx")
mito_proact_scores = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\ProAct_mean_scores_causality in tissues.csv")

mitoAge = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoAgeCo.xlsx")
mitoTSScore = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTS_scores_dataset.xlsx")
mitoPeakIntensity = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoPeak_Intensity_dataset.xlsx")

# Import development preferential values
mito_development_pref_time = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoDevelopment_preferential\final_time_preferential_values.csv")
mito_development_pref_tissue = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoDevelopment_preferential\final_tissue_preferential_values.csv")

# Import genomic features
#mitoCarta_dataset = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\mitocarta nominal\mitoCarta_final_dataset.xlsx")
#mitoCoreVar = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mito_core_variable_final_dataset.xlsx")
#mitoTSCategories = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\TS category\mitoTS_categories_dataset_binary.xlsx")
#mito_disCategories = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\category\category_binary.xlsx")
#mitoLocaliztion = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\mitochondrial localization\mitochondrial_localizations_top10_binary.xlsx")
#mitoGoSlim_cellular_component = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\GO slim cellular_component\cellular_component_top10_binary.xlsx")
#mitoGoSlim_biological_process = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\GO slim biological_process\biological_process_top10_binary.xlsx")
#mitoGoSlim_molecular_function = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\GO slim molecular_function\molecular_function_top10_binary.xlsx")

#GO_terms_enriched = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Panther enrichment analysis - mito vs prot cod\GO_terms_enriched_mitocarta_genes vs prot cod.csv")
#mitoPathways = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\mitoPathway\mitoPathways_dataset_binary.xlsx")

# Perform merges for most confident dataset
final_tissue_specific_df = merge_tissue_specific_features(labeled_variables_df, mitoTrace, pref_exp_mito, pref_num_ele_interactors, pref_num_interactors, pref_num_specific_interactions,mito_proact_scores, mitoAge, mitoTSScore, mitoPeakIntensity, mito_development_pref_time, mito_development_pref_tissue)
genomic_dataset = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla.csv")
#change gene ensg id for gene that has gene name SCO2
genomic_dataset.loc[genomic_dataset['Gene_x'] == 'SCO2', 'ENSG ID'] = 'ENSG00000130489'
#filter genomic dataset to columns that do not include 'preferential_expression', 'num_elevated_interactors', 'num_interactors', 'num_specific_interactions', 'prefexp', 'diffnetmed', 'egene', 'paralogsratiohighestidentity', 'mediantipapathways','developmentcv', 'mediancrisprscore', 'medianrnaiscore','TS_score', 'age coefficient', 'peak_intensity'
genomic_dataset_filtered = genomic_dataset.loc[:, ~genomic_dataset.columns.str.contains('preferential_expression|num_elevated_interactors|num_interactors|num_specific_interactions|prefexp|diffnetmed|egene|paralogsratiohighestidentity|proact|developmentcv|mediancrisprscore|medianrnaiscore|TS score|age coefficient|peak_intensity', case=False)]


final_tissue_specific_df = final_tissue_specific_df.merge(genomic_dataset_filtered, left_on=['Ensembl_ID_x'], right_on=['ENSG ID'], how='left')
final_tissue_specific_df = final_tissue_specific_df.drop(['ENSG ID', 'Gene_x', 'Label'], axis=1)

#final_whole_dataset = merge_genomic_features(final_tissue_specific_df, mitoCarta_dataset, mitoCoreVar, mitoTSCategories, mito_disCategories, GO_terms_enriched, mitoPathways)

# Perform the same steps for least confident dataset
# least_confidence_mitoGenes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Least confident dataset for mitocarta omim positive genes\Least_confidence_genes_manual_labels.xlsx")
# least_confidence_variables_df = create_variables(least_confidence_mitoGenes, mitoModel_tissues)
# least_confident_labels = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Least confident dataset for mitocarta omim positive genes\Least_confidence_genes_manual_labels.xlsx")
# least_confident_labels['manual label'] = 1
# least_confident_labels_unique = least_confident_labels.drop_duplicates(subset=['Gene_name', 'Candidate_tissue'])
# labeled_least_variables_df = least_confidence_variables_df.merge(least_confident_labels_unique, left_on=['Gene Name', 'Tissue'], right_on=['Gene_name', 'Candidate_tissue'], how='left')
# labeled_least_variables_df['label'] = labeled_least_variables_df['manual label'].fillna(0).astype(int)
# labeled_least_variables_df = labeled_least_variables_df.drop(['Gene_name', 'Candidate_tissue', 'Ensembl_ID_y', 'manual label'], axis=1)

# final_tissue_specific_least_df = merge_tissue_specific_features(labeled_least_variables_df, mitoTrace, pref_exp_mito, pref_num_ele_interactors, pref_num_interactors, pref_num_specific_interactions, mitoAge, mitoTSScore, mitoPeakIntensity, mito_development_pref_time, mito_development_pref_tissue)
# final_whole_least_dataset = merge_genomic_features(final_tissue_specific_least_df, mitoCarta_dataset, mitoCoreVar, mitoTSCategories, mito_disCategories, mitoLocaliztion, mitoGoSlim_cellular_component, mitoGoSlim_biological_process, mitoGoSlim_molecular_function, mitoPathways)

# Save the final datasets to Excel
final_tissue_specific_df.to_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confidence ds for mitocarta omim positive with pref\Most_confidence_manual_labels_with_features_dataset_new_all omim positive mitocarta_withprefGorilla.csv", index=False)
# final_whole_least_dataset.to_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Least confident dataset for mitocarta omim positive genes\Least_confidence_manual_labels_with_features_dataset_new_all omim positive mitocarta.xlsx", index=False)

print('meow')