import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, ks_2samp
from statsmodels.stats.multitest import multipletests

# Take MitoCarta genes and protein-coding genes and calculate the median expression for each tissue
mitocarta_genes = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Human.MitoCarta3.0.xls", sheet_name='A Human MitoCarta3.0')['EnsemblGeneID_mapping_version_20200130']
mitocarta_genes = mitocarta_genes.apply(lambda x: x.split('|')[0])

protein_coding_genes = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Characterization\PROTEIN CODING GENES.txt", sep=',')
protein_coding_genes = protein_coding_genes[~protein_coding_genes['Gene stable ID'].isin(mitocarta_genes)]

# Load GTEx dataset
gtex_file = r"C:\Users\liorr\Dropbox\Lior\Characterization\GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
df_gtex = pd.read_csv(gtex_file, sep='\t', skiprows=2)
df_gtex_filtered = df_gtex[df_gtex.iloc[:, 2:].ge(1).any(axis=1)]
df_gtex_filtered['Name'] = df_gtex_filtered['Name'].apply(lambda x: x.split('.')[0])

# Extract relevant data from GTEx for these gene sets
mito_expression = df_gtex_filtered[df_gtex_filtered['Name'].isin(mitocarta_genes)]
protein_expression = df_gtex_filtered[df_gtex_filtered['Name'].isin(protein_coding_genes['Gene stable ID'])]

# Prepare data for statistical tests and plotting for median expression
tissues = df_gtex.columns[2:]  # Assuming the first columns are gene information

# List to store the p-values and significance markers
p_values = []
significance_markers = []

# Calculate median expression and perform Mann-Whitney U tests
for tissue in tissues:
    mito_exp_values = mito_expression[tissue].dropna()
    protein_exp_values = protein_expression[tissue].dropna()

    # Perform the Mann-Whitney U test
    stat, p_value = mannwhitneyu(mito_exp_values, protein_exp_values, alternative='two-sided')
    p_values.append(p_value)

    # Print the results of the test
    print(f"Tissue: {tissue}, Test: Mann-Whitney U, Group 1 (Mito): n={len(mito_exp_values)}, "
          f"Group 2 (Protein Coding): n={len(protein_exp_values)}, p-value: {p_value}")

# Adjust p-values using Bonferroni correction
_, corrected_p_values, _, _ = multipletests(p_values, method='bonferroni')

# Determine significance markers after correction
for p_value in corrected_p_values:
    if p_value < 0.05:
        significance_markers.append('*')
    else:
        significance_markers.append('')

# Calculate median expression per tissue for each group
mito_median_expression = mito_expression.iloc[:, 2:].median(axis=0)
protein_median_expression = protein_expression.iloc[:, 2:].median(axis=0)

# Prepare data for plotting median expression
data = pd.DataFrame({
    'Tissue': tissues,
    'Mitochondrial Genes': mito_median_expression.values,
    'Protein Coding Genes': protein_median_expression.values
})
data_melted = data.melt(id_vars='Tissue', var_name='Gene Type', value_name='Expression')

# Set seaborn color palette for the first plot
sns.set_palette("colorblind")

# Plotting median expression
plt.figure(figsize=(15, 6))
ax = sns.barplot(x='Tissue', y='Expression', hue='Gene Type', data=data_melted)

# Set legend to the upper left
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)

# Set the title with some padding and move the asterisks above the bars
plt.title('Mitochondrial Genes vs Protein Coding Genes Expression Across Tissues (≥1 TPM)', pad=20)

# Add asterisks above the bars
for i, marker in enumerate(significance_markers):
    if marker:
        ax.text(i, max(data_melted['Expression']) + 1, marker, ha='center', va='bottom', color='black', fontsize=12)

# Rotate x-axis labels
plt.xticks(rotation=90)
plt.xlabel('Tissue')
plt.ylabel('Median Expression Level (TPM)')

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Expression across tissues mito vs pro_cod.png")
# Adjust layout to prevent overlap with the legend
plt.show()
plt.close()
# --------- PART 2: Distribution of Genes by Number of Expressing Tissues with KS Test --------- #

# Count how many tissues express each gene ≥1 TPM
df_gtex_filtered['Expressing Tissues Count'] = df_gtex_filtered.iloc[:, 2:].ge(1).sum(axis=1)

bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
labels = ['1-5', '6-10', '11-15', '16-20', '21-25', '26-30', '31-35', '36-40', '41-45', '46-50', '51-56']

df_gtex_filtered['Binned'] = pd.cut(df_gtex_filtered['Expressing Tissues Count'], bins=bins, labels=labels, right=True)

# Separate into different categories: mitochondrial and protein-coding
mito_distribution = df_gtex_filtered[df_gtex_filtered['Name'].isin(mitocarta_genes)]['Binned'].value_counts(sort=False)
protein_distribution = df_gtex_filtered[df_gtex_filtered['Name'].isin(protein_coding_genes['Gene stable ID'])]['Binned'].value_counts(sort=False)

# Calculate percentage of genes in each bin
mito_distribution_percent = (mito_distribution / mito_distribution.sum()) * 100
protein_distribution_percent = (protein_distribution / protein_distribution.sum()) * 100

# Perform Kolmogorov-Smirnov Test on the distributions
ks_stat, ks_p_value = ks_2samp(mito_distribution_percent, protein_distribution_percent)
print(f"Kolmogorov-Smirnov test between Mitochondrial and Protein-Coding Genes: KS Stat={ks_stat}, p-value={ks_p_value}")

# Get the same color palette used in the first plot
palette = sns.color_palette("colorblind")

# Plotting the distribution as percentage
plt.figure(figsize=(10, 6))
bar_width = 0.35
r1 = range(len(labels))

# Bar plots for mitochondrial and protein-coding genes with same colors as first plot
plt.bar(r1, mito_distribution_percent.values, color=palette[0], width=bar_width, label='Mitochondrial Genes')
plt.bar([x + bar_width for x in r1], protein_distribution_percent.values, color=palette[1], width=bar_width, label='Protein Coding Genes')

# Set x-ticks and labels
plt.xticks([r + bar_width / 2 for r in r1], labels, rotation=45)
plt.xlabel('Expressing tissues (#)')
plt.ylabel('Expressed genes (% of total)')
plt.title('Distribution of Genes by Number of Expressing Tissues (≥1 TPM)')

# Add legend and show plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
plt.tight_layout()
plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Distribution of Genes by Number of Expressing Tissues.png")
plt.show()
plt.close()
