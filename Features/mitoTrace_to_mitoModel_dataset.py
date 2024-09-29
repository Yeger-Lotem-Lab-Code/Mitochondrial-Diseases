import pandas as pd
import numpy as np


def create_tissues_file_dict(features_list, mito_tissues,tissues_file):
    tissues_file_dict = {}
    for feature in features_list:
        tissues_file_dict[feature] = {}
        for mito_tissue in mito_tissues:
            tissues_file_dict[feature][mito_tissue] = tissues_file.loc[tissues_file['mitomodel']==mito_tissue, feature]
    return tissues_file_dict


features_file = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\features.xlsx")
dataset_dict = {'Gene': [], 'Tissue':[]}
features_list = []
for f in features_file['feature type'][17:31]:
    dataset_dict[f.lower().replace(' ', '').replace('-', '').replace('_', '').replace('(', '').replace(')', '')] = []
    features_list.append(f.lower().replace(' ', '').replace('-', '').replace('_', '').replace('(', '').replace(')', ''))
tissues_file = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\creating dataset\trace\tissues_file_features_split_tissues_preprocessed.xlsx")
mito_tissues = tissues_file['mitomodel'].unique()
mito_tissues = [x for x in mito_tissues if str(x) != 'nan']
tissues_file_dict = create_tissues_file_dict(features_list, mito_tissues,tissues_file)
mito_Trace = pd.read_excel(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\creating dataset\trace\mitoTrace.xlsx")
for rowNum, row in mito_Trace.iterrows():
    for tissue in mito_tissues:
        dataset_dict['Gene'].append(row['geneid'])
        dataset_dict['Tissue'].append(tissue)
        for feature in features_list:
            feature_tissues = tissues_file_dict[feature][tissue]
            values_list = []
            for f_t in feature_tissues:
                if f_t==f_t:
                    value = row[f_t+feature]
                    if value==value:
                        values_list.append(value)
            if 3 < len(values_list):
                dataset_dict[feature].append(np.nanmedian(values_list))
            elif 0 < len(values_list) <= 3:
                dataset_dict[feature].append(np.nanmean(values_list))
            else:
                dataset_dict[feature].append("")

mitoTrace_dataset = pd.DataFrame.from_dict(dataset_dict)
mitoTrace_dataset.to_excel('mitoTrace_dataset.xlsx')
print('meow')