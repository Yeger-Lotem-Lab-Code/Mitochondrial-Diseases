import pandas as pd
#Least confident dataset
'''
dataset_Least_Conf = pd.read_excel(r"C:\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Least confident dataset for mitocarta omim positive genes\Least_confidence_manual_labels_with_features_dataset_new_all omim positive mitocarta.xlsx")
tissues_with_positive_labels_Least = dataset_Least_Conf[dataset_Least_Conf['label'] == 1]['Tissue'].unique()
filtered_subset_df_Least = dataset_Least_Conf[dataset_Least_Conf['Tissue'].isin(tissues_with_positive_labels_Least)]
num_positive_labels_Least = filtered_subset_df_Least[filtered_subset_df_Least['label'] == 1].shape[0]
print(f'Number of positive labels in Least confident dataset: {num_positive_labels_Least}')
filtered_subset_df_Least.to_excel(r"C:\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Least confident dataset for mitocarta omim positive genes\Least confident tissue-specific dataset with only positive tissues_new_all omim positive mitocarta.xlsx")
'''

#Most confident dataset
dataset_Most_Conf = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\Most_confidence_genes_manual_labels_with_features.csv")
tissues_with_positive_labels_Most = dataset_Most_Conf[dataset_Most_Conf['label'] == 1]['Tissue'].unique()
filtered_subset_df_Most = dataset_Most_Conf[dataset_Most_Conf['Tissue'].isin(tissues_with_positive_labels_Most)]
num_positive_labels_Most = filtered_subset_df_Most[filtered_subset_df_Most['label'] == 1].shape[0]
print(f'Number of positive labels in Most confident dataset: {num_positive_labels_Most}')
filtered_subset_df_Most.to_csv(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\Most confidence tissue-specific dataset with only positive tissues.csv")

print('meow')
