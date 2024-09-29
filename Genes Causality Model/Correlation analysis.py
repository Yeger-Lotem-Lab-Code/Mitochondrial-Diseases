import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

def infer_feature_type(name):
    # Convert feature name to lowercase for consistent comparisons
    name = name.lower()

    # Initialize a variable to store the inferred feature type
    new_name = ''

    # Define conditions for matching feature types
    if 'preferential_expression' in name:
        new_name = 'Tissue Preferential Expression'
    elif 'preferential_num_interactors' in name or 'preferential_num_elevated_interactors' in name or 'preferential_num_specific_interactions' in name or 'diffnetmed' in name:
        new_name = 'Tissue PPIs (Preferential)'
    elif 'proact' in name:
        new_name = 'Tissue ProACT Scores'
    elif '(neutral) median crispr score' in name or '(neutral) median rnai score' in name:
        new_name = 'Tissue Essentiality Scores'
    elif 'peak_intensity' in name:
        new_name = 'Tissue Peak Intensity'
    elif 'egene' in name:
        new_name = 'Tissue eQTL'
    elif 'paralogs_ratio_highest_identity' in name:
        new_name = 'Tissue Paralogs Relationships'
    elif 'ts score' in name:
        new_name = 'Tissue TS Score'
    elif 'age coefficient' in name:
        new_name = 'Tissue Age Coefficient'
    elif 'time' in name:
        new_name = 'Tissue Development Time'
    elif 'tissue' in name:
        new_name = 'Tissue Development Tissue'
    elif 'development_cv' in name:
        new_name = 'Tissue Development CV'
    # Return the inferred feature type
    return new_name

def extract_tissue_name(name):
    return name.split('_')[0]

def calculate_spearman_corr(df, feature_columns):
    # Calculate the Spearman correlation matrix
    corr_matrix = df[feature_columns].corr(method='spearman')
    return corr_matrix

def plot_correlation_heatmap(corr_matrix, title):
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=False, cmap='coolwarm', center=0, vmin=-1, vmax=1)
    plt.title(title, fontsize=16)
    plt.show()

df = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv")
feature_types = {col: infer_feature_type(col) for col in df.columns}
dftrans = df.T
dftrans['feature_type'] = dftrans.index.map(feature_types)

tissue_specific_feature_types = ['Tissue Preferential Expression', 'Tissue PPIs (Preferential)',
                                 'Tissue ProACT Scores', 'Tissue Essentiality Scores',
                                 'Tissue Peak Intensity', 'Tissue eQTL',
                                 'Tissue Paralogs Relationships', 'Tissue TS Score']

for feature_type in tissue_specific_feature_types:
    feature_rows = dftrans[dftrans['feature_type'] == feature_type].index.tolist()
    feature_rows_cleaned = [extract_tissue_name(row) for row in feature_rows]

    if len(feature_rows) > 1:
        # Calculate Spearman correlation
        corr_matrix = calculate_spearman_corr(df, feature_rows)
        corr_matrix.index = feature_rows_cleaned
        corr_matrix.columns = feature_rows_cleaned
        # Plot heatmap
        plot_correlation_heatmap(corr_matrix, f'Spearman Correlation for {feature_type} Across Tissues')
# Features to consider for the non-tissue-specific correlation analysis
non_tissue_specific_features = [
    'ProteinLength',
    'mitoPathway_Metabolism',
    'Cellular Component: GO:0043233 (organelle lumen)',
    'mitoPathway_OXPHOS'
]
corr_matrix_non_tissue = calculate_spearman_corr(df, non_tissue_specific_features)
plot_correlation_heatmap(corr_matrix_non_tissue, 'Spearman Correlation for Non-Tissue-Specific Features')


