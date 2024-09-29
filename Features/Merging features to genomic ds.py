import pandas as pd
import numpy as np


OMIM_Labels = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\Labeling\Labels_by_OMIM_MitoCarta.xlsx")
#tissues associations file
tissue_associations = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx", sheet_name=0)
#Add labels to the final dataset
'''
def merge_features(feature, final_ds, key_columns):
    feature = feature[['Gene'] + list(set(feature.columns) - set(['Gene']))]
    final_ds = pd.merge(final_ds, feature, left_on='ENSG', right_on='Gene', how='left')
    final_ds.drop(columns=['Gene'], inplace=True)
    return final_ds
'''
# Load datasets
# Preferential features
pref_exp = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_expression_mito_genes.xlsx")
pref_num_elevated_interactors = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_elevated_interactors_mito.xlsx")
pref_num_interactors = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_interactors_mito.xlsx")
pref_num_specific_interactions = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoTrace_preferential\preferential_num_specific_interactions_mito.xlsx")
pref_development_tissue = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoDevelopment_preferential\final_tissue_preferential_values.csv")
pref_development_time = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\mitoDevelopment_preferential\final_time_preferential_values.csv")
#merge all pref features with OMIM labels
OMIM_Labels = pd.merge(OMIM_Labels, pref_exp, left_on='ENSG ID', right_on='ENSG', how='left')
OMIM_Labels.drop(columns=['ENSG'], inplace=True)
OMIM_Labels = pd.merge(OMIM_Labels, pref_num_elevated_interactors, left_on='ENSG ID', right_on='Gene', how='left')
OMIM_Labels.drop(columns=['Gene_y'], inplace=True)

OMIM_Labels = pd.merge(OMIM_Labels, pref_num_interactors, left_on='ENSG ID', right_on='Gene', how='left')
OMIM_Labels.drop(columns=['Gene'], inplace=True)
OMIM_Labels = pd.merge(OMIM_Labels, pref_num_specific_interactions, left_on='ENSG ID', right_on='Gene', how='left')
OMIM_Labels.drop(columns=['Gene'], inplace=True)
OMIM_Labels = pd.merge(OMIM_Labels, pref_development_tissue, left_on='ENSG ID', right_on='Gene', how='left')
OMIM_Labels.drop(columns=['Gene'], inplace=True)
OMIM_Labels = pd.merge(OMIM_Labels, pref_development_time, left_on='ENSG ID', right_on='Gene', how='left')
OMIM_Labels.drop(columns=['Gene'], inplace=True)

#mitoTrace
mitotrace_file = pd.read_excel(r"C:\Users\liorr\Documents\pythonProject2\mitoTrace_genomic.xlsx")
#filter mitotrace file to only columns that the relevant features strings appear in
OMIM_Labels = pd.merge(OMIM_Labels, mitotrace_file, left_on='ENSG ID', right_on='Gene_ID', how='left')
OMIM_Labels.drop(columns=['Gene_ID'], inplace=True)
#ProAct scores
ProAct_scores = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\ProAct_mean_scores_causality in tissues.csv")
OMIM_Labels = pd.merge(OMIM_Labels, ProAct_scores, left_on='ENSG ID', right_on='Unnamed: 0', how='left')
OMIM_Labels.drop(columns=['Unnamed: 0'], inplace=True)

#GO terms enriched
GO_terms_enriched = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\Gorilla enrichment analysis manual choose\Gorilla enrichment analysis manual choose.csv")
OMIM_Labels = pd.merge(OMIM_Labels, GO_terms_enriched, left_on='Gene_x', right_on='Unnamed: 0', how='left')
OMIM_Labels.drop(columns=['Unnamed: 0'], inplace=True)
#mitoCarta features
mitocarta = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Human.MitoCarta3.0.xls", sheet_name='A Human MitoCarta3.0')
def check_in_ensembl_id(row, mitocarta_df):
    matching_row = mitocarta_df[mitocarta_df['EnsemblGeneID_mapping_version_20200130'].apply(lambda x: row['ENSG ID'] in str(x).split('|'))]
    if not matching_row.empty:
        return matching_row.iloc[0]
    return pd.Series([np.nan] * len(mitocarta.columns), index=mitocarta.columns)

# Apply the check_in_ensembl_id function
mitocarta_matches = OMIM_Labels.apply(lambda row: check_in_ensembl_id(row, mitocarta), axis=1)
mitocarta_matches_df = pd.DataFrame(mitocarta_matches)
mitocarta_matches_df.reset_index(drop=True, inplace=True)
OMIM_Labels.reset_index(drop=True, inplace=True)
specific_columns = ['ProteinLength', 'TargetP_Score', 'PGC_Induction_Score', 'CoexpressionGnfN50_Score']
mitocarta_specific = mitocarta_matches_df[specific_columns]
# Merge the specific columns with OMIM_Labels by gene
OMIM_Labels = pd.merge(OMIM_Labels, mitocarta_specific, left_index=True, right_index=True)



peak_intensity_columns = ['cerebrum_total_peak_intensity_log10', 'cerebellum_total_peak_intensity_log10',
                          'brainstem_total_peak_intensity_log10', 'spinalcord_total_peak_intensity_log10',
                          'kidney_total_peak_intensity_log10', 'liver_total_peak_intensity_log10',
                          'heart_total_peak_intensity_log10', 'skeletalmuscle_total_peak_intensity_log10',
                          'adipose_total_peak_intensity_log10', 'smallintestine_total_peak_intensity_log10',
                          'largeintestine_total_peak_intensity_log10', 'stomach_total_peak_intensity_log10',
                          'placenta_total_peak_intensity_log10', 'testis_total_peak_intensity_log10']

mitocarta_peak_intensity = mitocarta_matches_df[['EnsemblGeneID_mapping_version_20200130'] + peak_intensity_columns]

# DataFrame to store the selected and processed peak intensity columns
processed_peak_intensity = pd.DataFrame()

# Loop over the tissue associations and retain the relevant peak intensity columns
for i, row in tissue_associations.iterrows():
    peak_intensity_tissues = row['Peak Intensity']
    mito_model_tissue = row['mitoModel']

    if pd.isna(peak_intensity_tissues):
        continue

    associated_tissues = peak_intensity_tissues.split('|')

    # Retain only the columns matching the mitoModel tissues
    existing_tissues = [tissue for tissue in associated_tissues if tissue in mitocarta_peak_intensity.columns]

    if existing_tissues:
        if len(existing_tissues) == 1:
            # Rename the column to match mitoModel tissue
            mitocarta_peak_intensity.rename(columns={existing_tissues[0]: f"{mito_model_tissue}_peak_intensity"},
                                            inplace=True)
            processed_peak_intensity[f"{mito_model_tissue}_peak_intensity"] = mitocarta_peak_intensity[
                f"{mito_model_tissue}_peak_intensity"]
        elif len(existing_tissues) > 1:
            # Calculate the mean across the matching tissues
            processed_peak_intensity[f"{mito_model_tissue}_peak_intensity"] = mitocarta_peak_intensity[
                existing_tissues].mean(axis=1)

# Merge processed peak intensity columns with OMIM_Labels
OMIM_Labels = pd.merge(OMIM_Labels, processed_peak_intensity, left_index=True, right_index=True)


columns_to_span = ['YeastMitoHomolog_Score', 'RickettsiaHomolog_Score', 'hg19_Chromosome']

for column in columns_to_span:
    unique_values = mitocarta_matches[column].dropna().unique()
    for value in unique_values:
        new_column_name = f"{column}_{value}"
        OMIM_Labels[new_column_name] = np.where(mitocarta_matches[column] == value, 1, 0)

mitopathway = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\mitoPathway\mitoPathways_dataset_binary.xlsx")
OMIM_Labels = pd.merge(OMIM_Labels, mitopathway, left_on='Gene_x', right_on='Gene.1', how='left')
OMIM_Labels.drop(columns=['Gene','Gene.1'], inplace=True)


#TS features
TS_scores = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Tissue-Specifity score\NIHMS1624446-supplement-2.xlsx", sheet_name='G protein TS score', skiprows=2)
not_relevant_columns_TS = TS_scores.columns[1:4]
TS_scores.drop(columns=not_relevant_columns_TS, inplace=True)
for column in TS_scores.columns[1:]:
    if TS_scores[column].dtype == object:  # Only apply to object (string) columns
        TS_scores[column] = TS_scores[column].apply(lambda x: x.split(';')[0] if isinstance(x, str) else x)
processed_TS_scores = pd.DataFrame()
processed_TS_scores['ensembl_id'] = TS_scores['ensembl_id']  # Add the original ensembl_id column
for i, row in tissue_associations.iterrows():
    mitoModel_tissue = row['mitoModel']
    ts_tissues = row['TS Score']
    if ts_tissues == ts_tissues:
        associated_tissues = ts_tissues.split('|')
        ts_score_columns = TS_scores[associated_tissues].apply(pd.to_numeric, errors='coerce')
        if len(associated_tissues) == 1:
            new_column_name = f"{mitoModel_tissue}_TS score"
            processed_TS_scores[new_column_name] = ts_score_columns[associated_tissues[0]]
        elif len(associated_tissues) > 1:
            new_column_name = f"{mitoModel_tissue}_TS score"
            processed_TS_scores[new_column_name] = ts_score_columns.mean(skipna=True,axis=1)

OMIM_Labels = pd.merge(OMIM_Labels, processed_TS_scores, left_on='ENSG ID', right_on=processed_TS_scores['ensembl_id'], how='left')
OMIM_Labels.drop(columns=['ensembl_id'], inplace=True)
OMIM_Labels.drop_duplicates(inplace=True)
TS_categories = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Tissue-Specifity score\NIHMS1624446-supplement-2.xlsx", sheet_name='H comparison to other studies', skiprows=2)
unique_values_prot_enrichment = TS_categories['Our_enrichment_category'].dropna().unique()
for value in unique_values_prot_enrichment:
    new_column_name = f"protein enrichment_{value}"
    OMIM_Labels[new_column_name] = 0
for gene in OMIM_Labels['ENSG ID']:
    if gene in TS_categories['ensembl_id'].values:
        cat_value = TS_categories[TS_categories['ensembl_id'] == gene]['Our_enrichment_category'].values[0]
        if cat_value == cat_value:
            OMIM_Labels.loc[OMIM_Labels['ENSG ID'] == gene, f"protein enrichment_{cat_value}"] = 1
#Age Coefficients
age_coefficient_mouse = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Age Coefficient\41586_2020_2496_MOESM8_ESM (1).xlsx",sheet_name='DGE_result.tissue.age_coef')
age_coefficient_mouse_human_orthologes = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Age Coefficient\human_mouse_orthologs.xls")
#filter orthologes to only mitocarta genes
age_coefficient_mouse_human_orthologes_mito = age_coefficient_mouse_human_orthologes[age_coefficient_mouse_human_orthologes['Gene stable ID'].isin(OMIM_Labels['ENSG ID'])]
#Filter age coefficient mouse to only genes from age_coefficient_mouse_human_orthologes_mito
age_coefficient_mouse = age_coefficient_mouse[age_coefficient_mouse['Unnamed: 0'].isin(age_coefficient_mouse_human_orthologes_mito['Mouse gene name'])]
#add for each gene the Gene stable ID from age_coefficient_mouse_human_orthologes_mito in the age_coefficient_mouse dataframe
age_coefficient_human = pd.merge(age_coefficient_mouse, age_coefficient_mouse_human_orthologes_mito, left_on='Unnamed: 0', right_on='Mouse gene name', how='left')
age_coefficient_human.drop(columns=['Mouse gene name', 'Unnamed: 0', 'Mouse gene stable ID' ], inplace=True)
processed_age_coefficient = pd.DataFrame()
processed_age_coefficient['Gene stable ID'] = age_coefficient_human['Gene stable ID']
processed_age_coefficient['Gene name'] = age_coefficient_human['Gene name']

# Map tissues from the age coefficient file to mitoModel tissues using the tissue_associations file
for i, row in tissue_associations.iterrows():
    age_coef_tissues = row['age coefficient']
    mito_model_tissue = row['mitoModel']
    if age_coef_tissues==age_coef_tissues:
        associated_tissues = age_coef_tissues.split('|')
        ageco_columns = age_coefficient_human[associated_tissues].apply(pd.to_numeric, errors='coerce')

        if len(associated_tissues) == 1:
            new_column_name = f"{mito_model_tissue}_age coefficient"
            processed_age_coefficient[new_column_name] = ageco_columns[associated_tissues[0]]
        elif len(associated_tissues) > 1:
            new_column_name = f"{mito_model_tissue}_age coefficient"
            processed_age_coefficient[new_column_name] = ageco_columns.mean(skipna=True, axis=1)

numeric_columns = processed_age_coefficient.select_dtypes(include=[np.number]).columns.tolist() + ['Gene stable ID']
age_coefficient_human_grouped = processed_age_coefficient[numeric_columns].groupby('Gene stable ID').median().reset_index()

OMIM_Labels = pd.merge(OMIM_Labels, age_coefficient_human_grouped, left_on='ENSG ID', right_on='Gene stable ID', how='left')
OMIM_Labels.drop(columns=['Gene stable ID'], inplace=True)

#core-variable
core_var = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\core-variable\mito_genes_core-variable.xlsx")
#turn variable in the Status column of core var to 1 and core to 0
core_var['Status'] = np.where(core_var['Status'] == 'variable', 1, 0)
#rename Status column to core-variable
core_var.rename(columns={'Status': 'core-variable'}, inplace=True)
OMIM_Labels = pd.merge(OMIM_Labels, core_var, left_on='ENSG ID', right_on='Unnamed: 0', how='left')
OMIM_Labels.drop(columns=['Unnamed: 0', 'Gene_Name'], inplace=True)
#review category
category = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\creating dataset\top10 categorial features\category\category_binary.xlsx")
category.rename(columns=lambda x: f"{x}_rev_cat" if x != 'Unnamed: 0' else x, inplace=True)
OMIM_Labels = pd.merge(OMIM_Labels, category, left_on='ENSG ID', right_on='Unnamed: 0', how='left')
OMIM_Labels.drop(columns=['Unnamed: 0'], inplace=True)
#check for unique rows
OMIM_Labels.drop_duplicates(inplace=True)
OMIM_Labels.to_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla.csv", index=False)
print('meow')