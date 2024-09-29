import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Import data
genomic_ds = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv")
# Create df of mito genes and labels from genomic_ds
mitocarta_genes_and_labels = genomic_ds[['ENSG ID', 'Label']]
# Create df of mito genes and labels
dismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 1]['ENSG ID']
non_dismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 0]['ENSG ID']

TS_scores = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Tissue-Specifity score\NIHMS1624446-supplement-2.xlsx", sheet_name='G protein TS score', skiprows=2)
TS_scores = TS_scores.rename(columns={'ensembl_id': 'ENSG ID'})
not_relevant_columns_TS = TS_scores.columns[1:4]
TS_scores.drop(columns=not_relevant_columns_TS, inplace=True)
for column in TS_scores.columns[1:]:
    if TS_scores[column].dtype == object:  # Only apply to object (string) columns
        TS_scores[column] = TS_scores[column].apply(lambda x: x.split(';')[0] if isinstance(x, str) else x)

# Filter ts scores file to only mito genes
TS_scores = TS_scores[TS_scores['ENSG ID'].isin(genomic_ds['ENSG ID'])]
TS_features = TS_scores.columns
p_values = []
significance_markers = []
dismito_medians = []
nondismito_medians = []
tissues = []

for feature in TS_features[1:]:
    # Extract the tissue name from the column name (e.g., 'Brain_num_interactors' -> 'Brain')
    tissues.append(feature)

    # Calculate median expression for each group (dismito and nondismito) for the current PPI feature
    dismito_values = TS_scores[TS_scores['ENSG ID'].isin(dismito_genes)][feature].dropna()
    nondismito_values = TS_scores[TS_scores['ENSG ID'].isin(non_dismito_genes)][feature].dropna()

    # Convert values to numeric and drop non-numeric entries
    dismito_values = pd.to_numeric(dismito_values, errors='coerce').dropna()
    nondismito_values = pd.to_numeric(nondismito_values, errors='coerce').dropna()

    # Store median values for plotting
    dismito_medians.append(dismito_values.median())
    nondismito_medians.append(nondismito_values.median())

    # Perform the Mann-Whitney U test
    stat, p_value = mannwhitneyu(dismito_values, nondismito_values, alternative='two-sided')
    p_values.append(p_value)

    # Print the results of the test
    print(f"Tissue: {feature}, Test: Mann-Whitney U, Group 1 (Dismito): n={len(dismito_values)}, "
          f"Group 2 (Non-Dismito): n={len(nondismito_values)}, p-value: {p_value}")

# Adjust p-values for multiple testing correction using FDR
_, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

# Determine significance markers after correction
for p_value in corrected_p_values:
    if p_value < 0.05:
        significance_markers.append('*')
    else:
        significance_markers.append('')

# Prepare data for plotting median expression of PPI features across tissues
data = pd.DataFrame({
    'Tissue': tissues,
    'Dismito': dismito_medians,
    'Non-Dismito': nondismito_medians,
})
data_melted = data.melt(id_vars='Tissue', var_name='Gene Type', value_name='Median Value')

# Plot median expression for current PPI type across tissues
sns.set_palette("colorblind")
plt.figure(figsize=(12, 18))  # Increased height for horizontal plot
ax = sns.barplot(y='Tissue', x='Median Value', hue='Gene Type', data=data_melted, orient='h')

# Set legend and title
plt.legend(loc='upper right', bbox_to_anchor=(1, 1), ncol=1)
plt.title('Tissue Specificity Score Across Tissues for Dismito vs Non-Dismito Genes', pad=20)

# Add grid lines
plt.grid(axis='x', linestyle='--', alpha=0.7)

# Add asterisks above the bars
for i, (tissue, marker) in enumerate(zip(tissues, significance_markers)):
    if marker:
        max_value = data_melted[data_melted['Tissue'] == tissue]['Median Value'].max()
        ax.text(max_value + 0.1, i, marker, ha='center', va='bottom', color='black', fontsize=12)

# Rotate x-axis labels and show plot
plt.xlabel('Median TS Score')
plt.ylabel('Tissue')
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig(
    fr"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Tissue specificity scores across tissues dis mito vs non dis mito.png")

plt.show()
plt.close()

# Create a DataFrame to store tissues and their corrected p-values
significant_tissues = pd.DataFrame({
    'Tissue': tissues,
    'Corrected p-value': corrected_p_values
})

# Filter tissues with significant corrected p-values (p < 0.05)
significant_tissues = significant_tissues[significant_tissues['Corrected p-value'] < 0.05]

# Sort by corrected p-value in descending order (from max to min)
significant_tissues_sorted = significant_tissues.sort_values(by='Corrected p-value', ascending=False)

# Print the tissues and their corrected p-values
print("Tissues with significant corrected p-values (FDR < 0.05) sorted from max to min:")
print(significant_tissues_sorted.to_string(index=False))