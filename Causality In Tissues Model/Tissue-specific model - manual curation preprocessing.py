import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#preprocess aviv's dataset
aviv_dataset = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Aviv's genes.xlsx", sheet_name='Sheet1')
aviv_tissues_assoc = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Aviv's genes.xlsx", sheet_name='Sheet2')
aviv_dict_tissues = {}
for idx, row in aviv_tissues_assoc.iterrows():
    aviv_name = row['aviv name']
    model_name = row['model name']
    aviv_dict_tissues[aviv_name] = model_name


aviv_dataset.replace({"Candidate_tissue": aviv_dict_tissues}, inplace=True)

#preprocess Ori's dataset
ori_dataset = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Ori's genes.xlsx", sheet_name='Sheet1')
ori_tissues_assoc = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Ori's genes.xlsx", sheet_name='Sheet2')
ori_dict_tissues = {}
for idx, row in ori_tissues_assoc.iterrows():
    ori_name = row['ori name']
    model_name = row['model name']
    ori_dict_tissues[ori_name] = model_name

ori_dataset.replace({"Candidate_tissue": ori_dict_tissues}, inplace=True)

#preprocess Lior's dataset
lior_dataset = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Lior's genes.xlsx", sheet_name='Sheet1')
lior_tissues_assoc = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Lior's genes.xlsx", sheet_name='Sheet2')
lior_dict_tissues = {}
for idx, row in lior_tissues_assoc.iterrows():
    lior_name = row['lior name']
    model_name = row['model name']
    lior_dict_tissues[lior_name] = model_name

lior_dataset.replace({"Candidate_tissue": lior_dict_tissues}, inplace=True)

#preprocess mitocarta omim positive genes not in review
mitocartaOMIMpositivegenesnotinreview = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\MitoCarta OMIM positive genes not in review.xlsx", sheet_name='Sheet1')
mitocartaOMIMpositivegenesnotinreview_tissues_assoc = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\MitoCarta OMIM positive genes not in review.xlsx", sheet_name='Sheet2')
mitocartaOMIMpositivegenesnotinreview_dict_tissues = {}
for idx, row in mitocartaOMIMpositivegenesnotinreview_tissues_assoc.iterrows():
    mitocartaOMIMpositivegenesnotinreview_name = row['lior name']
    model_name = row['model name']
    mitocartaOMIMpositivegenesnotinreview_dict_tissues[mitocartaOMIMpositivegenesnotinreview_name] = model_name

mitocartaOMIMpositivegenesnotinreview.replace({"Candidate_tissue": mitocartaOMIMpositivegenesnotinreview_dict_tissues}, inplace=True)




merged_dataset = pd.concat([aviv_dataset, ori_dataset, lior_dataset,mitocartaOMIMpositivegenesnotinreview], ignore_index=True)
merged_dataset['Candidate_tissue'] = merged_dataset['Candidate_tissue'].str.replace(' ', '').str.lower().str.replace('-', '')
mitoModel_tissues = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitoModel_tissues.xlsx")['mitoModel']
mitoModel_tissues = mitoModel_tissues.str.replace(' ', '').str.lower().str.replace('-', '')
merged_dataset = merged_dataset[merged_dataset['Candidate_tissue'].isin(mitoModel_tissues)]
#filter to only genes that appear in mitocarta
mitocarta_genes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta\Gene_List_1.xlsx")['Gene Name']
merged_dataset = merged_dataset[merged_dataset['Gene_name'].isin(mitocarta_genes)]


most_confidence_dataset = merged_dataset[merged_dataset['Confidence_level_inflicted'] == 3]
least_confidence_dataset = merged_dataset[merged_dataset['Confidence_level_inflicted'].isin([3, 2])]
most_confidence_positive_genes = list(most_confidence_dataset['Gene_name'].unique())
least_confidence_positive_genes = list(least_confidence_dataset['Gene_name'].unique())

most_confidence_dataset = most_confidence_dataset[['Ensembl_ID','Gene_name', 'Candidate_tissue']]
least_confidence_dataset = least_confidence_dataset[['Ensembl_ID', 'Gene_name', 'Candidate_tissue']]


most_confidence_dataset.to_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confident dataset for mitocarta omim positive genes\Most_confidence_genes_manual_labels.xlsx", index=False)
least_confidence_dataset.to_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Least confident dataset for mitocarta omim positive genes\Least_confidence_genes_manual_labels.xlsx", index=False)

plt.figure(figsize=(12, 6))
sns.countplot(y='Candidate_tissue', data=most_confidence_dataset, order=most_confidence_dataset['Candidate_tissue'].value_counts().index)
plt.title('Number of Genes Associated with Each Tissue (Most Confidence)')
plt.xlabel('Number of Genes')
plt.ylabel('Tissue')
plt.tight_layout()
plt.show()

# Plot the number of genes associated with each tissue in the least_confidence_dataset
plt.figure(figsize=(12, 6))
sns.countplot(y='Candidate_tissue', data=least_confidence_dataset, order=least_confidence_dataset['Candidate_tissue'].value_counts().index)
plt.title('Number of Genes Associated with Each Tissue (Least Confidence)')
plt.xlabel('Number of Genes')
plt.ylabel('Tissue')
plt.tight_layout()
plt.show()

#Plot also for tissues in the OMIM positive genes, HPO labeled dataset
OMIMpositivegenesHPOlabeleddataset = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\OMIM Positive HPO labeled\MitoCarta positive genes by OMIM, HPO labeled dataset.xlsx")
filtered_df = OMIMpositivegenesHPOlabeleddataset[OMIMpositivegenesHPOlabeleddataset['Label'] == 1]

tissue_label_counts = OMIMpositivegenesHPOlabeleddataset[OMIMpositivegenesHPOlabeleddataset['Label'] == 1].groupby('Tissue').size().reset_index(name='Count')
plt.figure(figsize=(12, 6))
sns.countplot(y='Tissue', data=filtered_df, order=filtered_df['Tissue'].value_counts().index)
plt.title('Number of Genes with Label 1 in Each Tissue')
plt.xlabel('Number of Genes')
plt.ylabel('Tissue')
plt.tight_layout()
plt.show()
print('meow')

