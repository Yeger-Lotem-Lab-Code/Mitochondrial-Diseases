import pandas as pd

df = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\Most confidence tissue-specific dataset with only positive tissues.csv")

# Function to impute missing values and track genes that needed imputation
def impute_missing_values(df):
    one_columns = [col for col in df.columns if 'egene' in col]
    zero_columns = [
        col for col in df.columns if 'preferential_expression' in col or 'num_elevated_interactors' in col or
                                     'num_interactors' in col or 'num_specific_interactions' in col or 'diff_net_med' in col or
                                     'paralogs_ratio_highest_identity' in col or 'development_cv' in col or
                                     'prefexp' in col or 'TS score' in col or 'GO:' in col
    ]

    mode_columns = ['ProteinLength',
                    'TargetP_Score', 'PGC_Induction_Score', 'CoexpressionGnfN50_Score', 'core-variable'
                    ]
    # Append to the mode_columns columns where those strings appear in-'neutralmediancrisprscore', 'neutralmedianrnaiscore'
    mode_columns.extend([col for col in df.columns if
                         'crispr' in col or 'rnai' in col or 'proact' in col or 'peak_intensity' in col or 'age coefficient' in col or 'rev_cat' in col or 'protein enrichment' in col])
    imputation_counts = pd.DataFrame(columns=['Feature', 'Imputations', 'Gene-Tissue_Pairs_Needed_Imputation', 'Imputation_Type'])
# Impute one values
    for column in one_columns:
        if column in df.columns:
            value = 1
            missing_pairs = df.loc[df[column].isna(), ['Gene Name', 'Tissue']].apply(tuple, axis=1).tolist()
            num_imputations = len(missing_pairs)
            df[column].fillna(value, inplace=True)
            imputation_counts = pd.concat([imputation_counts, pd.DataFrame({'Feature': [column], 'Imputations': [num_imputations], 'Gene-Tissue_Pairs_Needed_Imputation': [missing_pairs], 'Imputation_Type': ['One']})])

    # Impute zero values
    for column in zero_columns:
        if column in df.columns:
            missing_pairs = df.loc[df[column].isna(), ['Gene Name', 'Tissue']].apply(tuple, axis=1).tolist()
            num_imputations = len(missing_pairs)
            df[column].fillna(0, inplace=True)
            imputation_counts = pd.concat([imputation_counts, pd.DataFrame({'Feature': [column], 'Imputations': [num_imputations], 'Gene-Tissue_Pairs_Needed_Imputation': [missing_pairs], 'Imputation_Type': ['Zero']})])

    # Impute mode values
    for column in mode_columns:
        if column in df.columns:
            mode_value = df[column].mode()[0]
            missing_pairs = df.loc[df[column].isna(), ['Gene Name', 'Tissue']].apply(tuple, axis=1).tolist()
            num_imputations = len(missing_pairs)
            df[column].fillna(mode_value, inplace=True)
            imputation_counts = pd.concat([imputation_counts, pd.DataFrame({'Feature': [column], 'Imputations': [num_imputations], 'Gene-Tissue_Pairs_Needed_Imputation': [missing_pairs], 'Imputation_Type': ['Mode']})])

    return df, imputation_counts

df_imputed, imputation_counts = impute_missing_values(df)
df_imputed.to_csv(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\Most confidence tissue-specific dataset with only positive tissues_imputed_dataset.csv", index=False)
imputation_counts.to_csv(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\Most confidence tissue-specific dataset with only positive tissues_imputation_counts.csv", index=False)
print('meow')
