import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ks_2samp


# Load your dataset
dataset = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv")

# Define the top 10 influential features identified from SHAP analysis
top_10_features = [
     'ProteinLength',
    'Colon1_mean proact',
    'mitoPathway_Metabolism',
    'Cellular Component: GO:0043233 (organelle lumen)',
    'Stomach_TS score',
    'Whole Blood_mean proact',
    'Whole Esophagus_mean proact',
    'Artery2_mean proact',
    'mitoPathway_OXPHOS',
    'Pituitary_mean proact'
]

# Separate the dataset into causal and non-causal groups
causal_genes = dataset[dataset['Label'] == 1]
non_causal_genes = dataset[dataset['Label'] == 0]

# Plot the distribution of each top 10 feature for both groups
plt.figure(figsize=(20, 30))
significant_features = []

for i, feature in enumerate(top_10_features, 1):
    plt.subplot(5, 2, i)  # 5 rows, 2 columns
    median_causal = causal_genes[feature].median()
    median_non_causal = non_causal_genes[feature].median()

    sns.kdeplot(causal_genes[feature], label=f'Causal (median={median_causal:.2f})', shade=True, color='blue')
    sns.kdeplot(non_causal_genes[feature], label=f'Non-Causal (median={median_non_causal:.2f})', shade=True, color='red')

    # Add median lines
    plt.axvline(median_causal, color='blue', linestyle='--')
    plt.axvline(median_non_causal, color='red', linestyle='--')

    # KS2 test
    ks_stat, p_value = ks_2samp(causal_genes[feature], non_causal_genes[feature])
    if p_value < 0.05:
        significant_features.append(feature)
    plt.title(f'Distribution of {feature}')
    plt.xlabel(feature)
    plt.ylabel('Density')
    plt.legend()

plt.tight_layout()
plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\5 fold cross validation - ds with pref and gorilla\top_10_features_distribution.png")
plt.show()

# Print significant features
if significant_features:
    print("Significant features (p < 0.05):")
    for feature in significant_features:
        print(feature)
else:
    print("No significant features found.")