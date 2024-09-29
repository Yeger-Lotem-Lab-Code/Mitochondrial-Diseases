import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as stats


#import manual curated list of mito genes
LiorsGenes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Lior's genes.xlsx", sheet_name='Sheet1')
AvivsGenes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Aviv's genes.xlsx", sheet_name='Sheet1')
OrisGenes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Ori's genes.xlsx", sheet_name='Sheet1')
# Load the dictionaries for tissue name mappings
LiorsDict = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Lior's genes.xlsx", sheet_name='Sheet2')
AvivsDict = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Aviv's genes.xlsx", sheet_name='Sheet2')
OrisDict = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Ori's genes.xlsx", sheet_name='Sheet2')
# Create dictionaries from the 'lior name' and 'model name' columns
LiorMapping = dict(zip(LiorsDict['lior name'].str.lower(), LiorsDict['model name']))
AvivMapping = dict(zip(AvivsDict['aviv name'].str.lower(), AvivsDict['model name']))
OriMapping = dict(zip(OrisDict['ori name'].str.lower(), OrisDict['model name']))

# Concatenate all lists
allGenes = pd.concat([LiorsGenes, AvivsGenes, OrisGenes], ignore_index=True)
# Merge the dictionaries into one mapping (keys are in lowercase)
tissue_mapping = {**LiorMapping, **AvivMapping, **OriMapping}
# Manually add mappings that were not covered
tissue_mapping.update({
    'skeletal - muscle': 'Muscle - Skeletal',
    'whole colon': 'Whole Colon'
})

# Convert the Candidate_tissue column to lowercase to ensure case-insensitive matching
allGenes['Candidate_tissue'] = allGenes['Candidate_tissue'].str.lower().str.strip()
# Replace the Candidate_tissue values with the corresponding model names
allGenes['Candidate_tissue'] = allGenes['Candidate_tissue'].replace(tissue_mapping)
mitocartaOMIMpositivegenesnotinreview = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\MitoCarta OMIM positive genes not in review.xlsx", sheet_name='Sheet1')
mitocartaOMIMpositivegenesnotinreview_tissues_assoc = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\MitoCarta OMIM positive genes not in review.xlsx", sheet_name='Sheet2')
mitocartaOMIMpositivegenesnotinreview_dict_tissues = {}
for idx, row in mitocartaOMIMpositivegenesnotinreview_tissues_assoc.iterrows():
    mitocartaOMIMpositivegenesnotinreview_name = row['lior name']
    model_name = row['model name']
    mitocartaOMIMpositivegenesnotinreview_dict_tissues[mitocartaOMIMpositivegenesnotinreview_name] = model_name

mitocartaOMIMpositivegenesnotinreview.replace({"Candidate_tissue": mitocartaOMIMpositivegenesnotinreview_dict_tissues}, inplace=True)
allGenes = pd.concat([allGenes, mitocartaOMIMpositivegenesnotinreview], ignore_index=True)
#filter to only genes that appear in mitocarta
mitocarta_genes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta\Gene_List_1.xlsx")['Gene Name']
allGenes = allGenes[allGenes['Gene_name'].isin(mitocarta_genes)]
allGenes.drop_duplicates(inplace=True)
#remove rows with confidence level lower than 3
allGenes = allGenes[allGenes['Confidence_level_inflicted'] >= 3]
#remove rows where tissues are not in the model
model_tissues = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx", sheet_name='Sheet1')['mitoModel'].str.strip()
allGenes = allGenes[allGenes['Candidate_tissue'].isin(model_tissues)]

# ------------------- FIRST PLOT: Number of Diseases Associated with Tissues -------------------
# Group by tissue and count the unique diseases, filtering only for model tissues
tissue_disease_count = allGenes.groupby('Candidate_tissue')['Disease_ID'].nunique().reset_index()
# Rename columns for clarity
tissue_disease_count.columns = ['Tissue', 'Number of Diseases']
# Sort by Number of Diseases in descending order
tissue_disease_count = tissue_disease_count.sort_values(by='Number of Diseases', ascending=False)

plt.figure(figsize=(16, 8))
plt.bar(tissue_disease_count['Tissue'], tissue_disease_count['Number of Diseases'], color='red')
# Formatting the plot to match the example given
plt.title('Number of Diseases Associated with Tissues', fontsize=16)
plt.ylabel('Number of Diseases', fontsize=14)
plt.xticks(rotation=90, ha='center', fontsize=10)
plt.tight_layout()
#plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Number of Diseases Associated with Tissues.png")
plt.show()

# ------------------- SECOND PLOT: Number of Diseases Grouped by Number of Inflicted Tissues -------------------
# Group by disease and count the number of unique tissues it affects
disease_tissue_count = allGenes.groupby('Disease_ID')['Candidate_tissue'].nunique().reset_index()
disease_tissue_count.columns = ['Disease_ID', 'Number_of_Tissues']

# Count the number of diseases that affect a certain number of tissues
tissue_count_distribution = disease_tissue_count['Number_of_Tissues'].value_counts().reset_index()
tissue_count_distribution.columns = ['Number_of_Tissues', 'Number_of_Diseases']
tissue_count_distribution = tissue_count_distribution.sort_values(by='Number_of_Tissues')

# Plot the number of diseases grouped by the number of tissues they inflict
plt.figure(figsize=(8, 6))
plt.bar(tissue_count_distribution['Number_of_Tissues'], tissue_count_distribution['Number_of_Diseases'], color='gray')

# Formatting the plot to match the example given
plt.title('Number of Diseases Grouped by Number of Inflicted Tissues', fontsize=16)
plt.xlabel('Number of Inflicted Tissues', fontsize=14)
plt.ylabel('Number of Diseases', fontsize=14)
plt.xticks(tissue_count_distribution['Number_of_Tissues'])
plt.tight_layout()

# Save and display the plot
#plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Number of Diseases Grouped by Inflicted Tissues.png")
plt.show()

# ------------------- Third PLOT: Gene Expression (TPM) in Inflicted vs Non-inflicted Tissues -------------------
gtex_file = r"C:\Users\liorr\Dropbox\Lior\Characterization\GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
gtex_data = pd.read_csv(gtex_file, sep='\t', skiprows=2)
expression_data = gtex_data.set_index('Name')
model_tissues = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx", sheet_name='mitoTrace')

# Split the expression data indices with '.' and take the first part
expression_data.index = expression_data.index.str.split('.').str[0]

# Filter the expression data to only include genes that are in the allGenes dataframe
expression_data = expression_data[expression_data.index.isin(allGenes['Ensembl_ID'])]

# Create a dictionary with mitomodel tissues as keys and the corresponding columns in the expression data as values
expression_df = pd.DataFrame()
for tissue in model_tissues['mitoModel']:
    corr_tissues = model_tissues[model_tissues['mitoModel'] == tissue]['expression|num_interactors|num_elevated_interactors|num_specific_interactions'].values
    expression_df[tissue] = expression_data[expression_data.columns.intersection(corr_tissues)].mean(axis=1)

# Identify inflicted tissues
inflicted_tissues = list(tissue_disease_count['Tissue'].unique())

# Calculate inflicted and non-inflicted median expression values for each gene
for rowNum, row in expression_df.iterrows():
    val_inflicted = []
    val_non_inflicted = []
    for item in row.index:
        if item in inflicted_tissues:
            val_inflicted.append(row[item])
        else:
            val_non_inflicted.append(row[item])
    # Add the values to the dataframe
    expression_df.loc[rowNum, 'inflicted'] = np.mean(val_inflicted)
    expression_df.loc[rowNum, 'non_inflicted'] = np.mean(val_non_inflicted)

# Create a combined DataFrame for plotting
plot_data = expression_df[['inflicted', 'non_inflicted']].reset_index()
plot_data = plot_data.melt(id_vars='Name', var_name='Group', value_name='Gene Expression (TPM)')

upper_threshold = plot_data['Gene Expression (TPM)'].quantile(0.9)
filtered_data = plot_data[plot_data['Gene Expression (TPM)'] <= upper_threshold]
inflicted_expression = filtered_data[filtered_data['Group'] == 'inflicted']['Gene Expression (TPM)']
non_inflicted_expression = filtered_data[filtered_data['Group'] == 'non_inflicted']['Gene Expression (TPM)']

# Perform the Mann-Whitney U Test
stat, p_value = stats.mannwhitneyu(inflicted_expression, non_inflicted_expression, alternative='two-sided')

# Print the test result
print(f"Mann-Whitney U Test Statistic: {stat}")
print(f"P-value: {p_value}")

# Plot using seaborn boxplot
plt.figure(figsize=(8, 6))
sns.boxplot(x='Group', y='Gene Expression (TPM)', data=filtered_data, palette=['#D62728', '#7F7F7F'])

# Add plot title and labels
plt.title('Gene Expression (TPM) in Inflicted vs Non-inflicted Tissues', fontsize=16)
plt.ylabel('Gene expression (TPM)', fontsize=14)

# Display the number of tissues in each group using len() instead of .shape
inflicted_count = len(inflicted_tissues)
non_inflicted_count = expression_df.shape[1] - inflicted_count
plt.xticks([0, 1], [f'Inflicted\n(n = {inflicted_count:,})', f'Non-inflicted\n(n = {non_inflicted_count:,})'])

# Save and show the plot
plt.tight_layout()
#plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Inflicted_vs_Non_Inflicted_Gene_Expression.png")
plt.show()

# ------------------- Fourth PLOT: Number of genes associated to each tissue -------------------
# Count the number of genes associated with each tissue
gene_tissues_count = allGenes.groupby('Candidate_tissue')['Gene_name'].nunique().reset_index()
# Rename columns for clarity
gene_tissues_count.columns = ['Tissue', 'Number of Genes']
# Sort by Number of Diseases in descending order
gene_tissues_count = gene_tissues_count.sort_values(by='Number of Genes', ascending=False)

plt.figure(figsize=(16, 8))
plt.bar(gene_tissues_count['Tissue'], gene_tissues_count['Number of Genes'], color='blue')
# Formatting the plot to match the example given
plt.title('Number of Genes Associated with Tissues', fontsize=16)
plt.ylabel('Number of Genes', fontsize=14)
plt.xticks(rotation=90, ha='center', fontsize=10)
plt.tight_layout()
plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Number of Genes Associated with Tissues.png")
plt.show()