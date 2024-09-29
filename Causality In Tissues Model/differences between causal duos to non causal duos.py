import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp

dataset_path = r"C:\Users\liorr\Dropbox\Lior\MITOCH~1\MANUAL~1\MOSTCO~4\MOSTCO~2.CSV"
df = pd.read_csv(dataset_path).iloc[:, 1:]# Columns: 'gene', 'tissue', 'label', 'Preferential_Expression', 'TS_score', 'Time_Childhood', ...
for column in df.columns:
    if column not in ['Gene Name', 'Ensembl_ID_x', 'Tissue', 'label']:  # Skip the label and ID columns
        max_value = df[column].replace([np.inf, -np.inf], np.nan).max()
        median_value = df[column].median()
        df[column].replace([np.inf, -np.inf], max_value, inplace=True)
        df[column].fillna(median_value, inplace=True)
top_10_features = [
    'Num_Elevated_Interactors',
    'TS score',
    'Preferential_Expression',
    'Age Coefficient',
    'egene',
    'Time_Childhood',
    'diffnetmed',
    'Time_Adolescent',
    'Time_Elderly',
    'Peak_Intensity'
]

# Separate the dataset into causal and non-causal groups
causal_duos = df[df['label'] == 1]
non_causal_duos = df[df['label'] == 0]

# Plot the distribution of each top 10 feature for both groups
plt.figure(figsize=(20, 30))

significant_features = []

for i, feature in enumerate(top_10_features, 1):
    plt.subplot(5, 2, i)  # 5 rows, 2 columns

    median_causal = causal_duos[feature].median()
    median_non_causal = non_causal_duos[feature].median()

    sns.kdeplot(causal_duos[feature], label=f'Causal (median={median_causal:.2f})', shade=True, color='blue')
    sns.kdeplot(non_causal_duos[feature], label=f'Non-Causal (median={median_non_causal:.2f})', shade=True, color='red')

    # Add median lines
    plt.axvline(median_causal, color='blue', linestyle='--')
    plt.axvline(median_non_causal, color='red', linestyle='--')

    # KS2 test
    ks_stat, p_value = ks_2samp(causal_duos[feature], non_causal_duos[feature])
    if p_value < 0.05:
        significant_features.append(feature)

    plt.title(f'Distribution of {feature}')
    plt.xlabel(feature)
    plt.ylabel('Density')
    plt.legend()

plt.tight_layout()
plt.savefig(
    r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confidence ds for mitocarta omim positive with pref\Comparison_Plots_prefGorilla\top_10_features_distribution.png")
plt.show()

# Print significant features
if significant_features:
    print("Significant features (p < 0.05):")
    for feature in significant_features:
        print(feature)
else:
    print("No significant features found.")
