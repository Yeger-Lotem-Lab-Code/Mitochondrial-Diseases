import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, ks_2samp
from statsmodels.stats.multitest import multipletests

# Take MitoCarta genes and protein-coding genes and calculate the median expression for each tissue
mitocarta_genes_and_labels = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\Labeling\Labels_by_OMIM_MitoCarta.xlsx")
dismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 1]['ENSG ID']
nondismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 0]['ENSG ID']

# Load GTEx dataset
gtex_file = r"C:\Users\liorr\Dropbox\Lior\Characterization\GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
df_gtex = pd.read_csv(gtex_file, sep='\t', skiprows=2)
df_gtex_filtered = df_gtex[df_gtex.iloc[:, 2:].ge(1).any(axis=1)]
df_gtex_filtered['Name'] = df_gtex_filtered['Name'].apply(lambda x: x.split('.')[0])

# Extract relevant data from GTEx for these gene sets
dismito_expression = df_gtex_filtered[df_gtex_filtered['Name'].isin(dismito_genes)]
nondismito_expression = df_gtex_filtered[df_gtex_filtered['Name'].isin(nondismito_genes)]

tissues = df_gtex.columns[2:]

# List to store the p-values and significance markers
p_values = []
significance_markers = []

# Calculate median expression and perform Mann-Whitney U tests
for tissue in tissues:
    dismito_exp_values = dismito_expression[tissue].dropna()
    nondismito_exp_values = nondismito_expression[tissue].dropna()

    # Perform the Mann-Whitney U test
    stat, p_value = mannwhitneyu(dismito_exp_values, nondismito_exp_values, alternative='two-sided')
    p_values.append(p_value)

    # Print the results of the test
    print(f"Tissue: {tissue}, Test: Mann-Whitney U, Group 1 (Dismito): n={len(dismito_exp_values)}, "
          f"Group 2 (Non-Dismito): n={len(nondismito_exp_values)}, p-value: {p_value}")

# Adjust p-values using Bonferroni correction
_, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

# Determine significance markers after correction
for p_value in corrected_p_values:
    if p_value < 0.05:
        significance_markers.append('*')
    else:
        significance_markers.append('')
# Calculate median expression per tissue for each group
dismito_median_expression = dismito_expression.iloc[:, 2:].median(axis=0)
nondismito_median_expression = nondismito_expression.iloc[:, 2:].median(axis=0)
# Prepare data for plotting median expression
data = pd.DataFrame({
    'Tissue': tissues,
    'Dismito': dismito_median_expression,
    'Non-Dismito': nondismito_median_expression,
})
data_melted = data.melt(id_vars='Tissue', var_name='Gene Type', value_name='Expression')
sns.set_palette("colorblind")
plt.figure(figsize=(15, 6))

ax = sns.barplot(x='Tissue', y='Expression', hue='Gene Type', data=data_melted)

# Set legend to the upper left
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)

# Set the title with some padding and move the asterisks above the bars
plt.title('Disease Mitochondrial Genes vs Non Disease Mitochondrial Genes Expression Across Tissues (≥1 TPM)', pad=20)

# Add asterisks above the bars
for i, marker in enumerate(significance_markers):
    if marker:
        ax.text(i, max(data_melted['Expression']) + 1, marker, ha='center', va='bottom', color='black', fontsize=12)

# Rotate x-axis labels
plt.xticks(rotation=90)
plt.xlabel('Tissue')
plt.ylabel('Median Expression Level (TPM)')

plt.tight_layout(rect=[0, 0, 0.85, 1])
#plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Expression across tissues dis mito vs non dis mito.png")
# Adjust layout to prevent overlap with the legend
plt.show()
plt.close()

# --------- PART 2: Distribution of Genes by Number of Expressing Tissues with KS Test --------- #
# Count how many tissues express each gene ≥1 TPM
df_gtex_filtered['Expressing Tissues Count'] = df_gtex_filtered.iloc[:, 2:].ge(1).sum(axis=1)

bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
labels = ['1-5', '6-10', '11-15', '16-20', '21-25', '26-30', '31-35', '36-40', '41-45', '46-50', '51-56']

df_gtex_filtered['Binned'] = pd.cut(df_gtex_filtered['Expressing Tissues Count'], bins=bins, labels=labels, right=True)
# Separate into different categories: DisMito and non-DisMito
dismito_distribution = df_gtex_filtered[df_gtex_filtered['Name'].isin(dismito_genes)]['Binned'].value_counts(sort=False)
nondismito_distribution = df_gtex_filtered[df_gtex_filtered['Name'].isin(nondismito_genes)]['Binned'].value_counts(sort=False)

# Calculate percentage of genes in each bin
dismito_distribution_percent = (dismito_distribution / dismito_distribution.sum()) * 100
nondismito_distribution_percent = (nondismito_distribution / nondismito_distribution.sum()) * 100


# Perform Kolmogorov-Smirnov Test on the distributions
ks_stat, ks_p_value = ks_2samp(dismito_distribution_percent, nondismito_distribution_percent)

print(f"Kolmogorov-Smirnov test between Dis-Mitochondrial and non Dis-Mitochondrial Genes: KS Stat={ks_stat}, p-value={ks_p_value}")

# Get the same color palette used in the first plot
palette = sns.color_palette("colorblind")

# Plotting the distribution as percentage
plt.figure(figsize=(10, 6))
bar_width = 0.35
r1 = range(len(labels))

# Bar plots for mitochondrial and protein-coding genes with same colors as first plot
plt.bar(r1, dismito_distribution_percent.values, color=palette[0], width=bar_width, label='Dismito')
plt.bar([x + bar_width for x in r1], nondismito_distribution_percent.values, color=palette[1], width=bar_width, label='Non-Dismito')

# Set x-ticks and labels
plt.xticks([r + bar_width / 2 for r in r1], labels, rotation=45)
plt.xlabel('Expressing tissues (#)')
plt.ylabel('Expressed genes (% of total)')
plt.title('Distribution of Genes by Number of Expressing Tissues (≥1 TPM)')

# Add legend and show plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
plt.tight_layout()
#plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Distribution of Disease and Non Disease Mitochondrial Genes by Number of Expressing Tissues.png")
plt.show()
plt.close()


