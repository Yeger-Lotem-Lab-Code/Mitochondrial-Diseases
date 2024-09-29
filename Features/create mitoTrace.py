import pandas as pd
import re


def clean_column_names(df):
    df.columns = [re.sub(r'\W+', '', col) for col in df.columns]
    return df


def escape_special_chars(text):
    # Escape special characters in feature names for regex
    return re.escape(text)


def extract_tissue_features(trace_file, feature_suffix):
    # Escape special characters in the feature_suffix
    escaped_suffix = escape_special_chars(feature_suffix)

    # Filter columns based on the escaped feature suffix
    feature_data = trace_file.filter(regex=f'Gene_ID|.*_{escaped_suffix}$')

    # Remove suffix from the column names
    feature_data.columns = [col.split('_')[0] if '_' in col else col for col in feature_data.columns]

    return feature_data


def add_whole_tissues_columns(feature, tissue_associations, feature_suffix):
    # Clean the column names in feature
    feature = clean_column_names(feature)

    previous_columns_names = feature.columns[1:]  # Exclude 'Gene_ID'

    for tissue in tissue_associations['mitoModel'].unique():
        associated_tissues = tissue_associations.loc[
            tissue_associations['mitoModel'] == tissue, feature_suffix].dropna().tolist()

        if associated_tissues:
            # Clean the associated tissues column names
            associated_tissues = [re.sub(r'\W+', '', col).lower() for col in associated_tissues]

            # Ensure that associated_tissues columns exist in the feature DataFrame
            existing_tissues = [col for col in associated_tissues if col in feature.columns]

            if len(existing_tissues) == 1:
                # If only one matching tissue, copy the column
                feature[tissue] = feature[existing_tissues[0]]
            elif len(existing_tissues) > 1:
                # If multiple matching tissues, compute the median
                feature[tissue] = feature[existing_tissues].median(axis=1)

    # Drop the original tissue columns
    feature = feature.drop(columns=previous_columns_names)
    feature.columns = [col if col == 'Gene_ID' else f"{col}_{feature_suffix}" for col in feature.columns]
    return feature


tissue_associations = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx", 'mitoTrace')

mitochondrial_genes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta plus review\merged_mitochondrial_genes_MitoCarta and review.xlsx")['First ENSG ID']
Trace = pd.read_csv(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Aneuploidy\5_TRACE_Dataset_with_RNAi.csv")

#filter trace to mitochondrial genes only
mitoTrace = Trace[Trace['Gene_ID'].isin(mitochondrial_genes)]
relevant_features = ['Gene_ID','diff_net_med', 'egene', 'paralogs_ratio_highest_identity', 'development_cv', '(neutral) median crispr score', '(neutral) median rnai score']
matching_columns = [col for col in mitoTrace.columns if any(feature in col for feature in relevant_features)]
mitoTrace = mitoTrace[matching_columns]

merged_features = mitoTrace[['Gene_ID']].copy()

for feature_suffix in relevant_features[1:]:
    feature = extract_tissue_features(mitoTrace, feature_suffix)
    feature_with_whole_columns = add_whole_tissues_columns(feature, tissue_associations, feature_suffix)
    merged_features = pd.merge(merged_features, feature_with_whole_columns, left_on='Gene_ID', right_on='Gene_'+f"{feature_suffix}", how='left')
    merged_features.drop(columns=['Gene_'+f"{feature_suffix}"], inplace=True)
merged_features.to_excel('mitoTrace_genomic.xlsx', index=False)


print('meow')

