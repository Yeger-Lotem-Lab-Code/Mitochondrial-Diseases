import pandas as pd
import numpy as np

# Load the data
file_path = r"C:\Users\liorr\Dropbox\Lior\development feature\Human_rpkm.txt"
gene_list_path = r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta\Gene_List_1.xlsx"
tissue_associations_path = r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx"

# Read the first line to get the column names
with open(file_path, 'r') as file:
    columns_line = file.readline().strip()
    columns = columns_line.split('" "')
    columns[0] = columns[0].replace('Names "', '')
    columns[-1] = columns[-1].replace('"', '')

# Read the rest of the data into a DataFrame
data = pd.read_csv(file_path, sep='\t', skiprows=1, header=None)

# Split the first column into ENSG and expression values
data_split = data[0].str.split(' ', n=1, expand=True)
data_split.columns = ['ENSG', 'Expression_Values']

# Split the expression values into separate columns
expression_values = data_split['Expression_Values'].str.split(' ', expand=True)
expression_values.columns = columns[1:]

# Combine the ENSG column with the expression values
data_cleaned = pd.concat([data_split['ENSG'], expression_values], axis=1)

# Melt the data to long format
data_long = pd.melt(data_cleaned, id_vars=['ENSG'], var_name='Tissue_Time', value_name='Expression_Value')

# Convert Expression_Value to numeric
data_long['Expression_Value'] = pd.to_numeric(data_long['Expression_Value'], errors='coerce')

# Split 'Tissue_Time' into 'Tissue' and 'Time_Point'
data_long[['Tissue', 'Time_Point', 'Replicate']] = data_long['Tissue_Time'].str.extract(r'([A-Za-z]+)\.(\d+\w+|\w+)\.(\d+)')

# Drop the unnecessary 'Replicate' column
data_long.drop(columns=['Tissue_Time', 'Replicate'], inplace=True)

def categorize_time_point(tp):
    if 'wpc' in tp:
        return 'Antenatal'
    elif tp == 'newborn':
        return 'Neonatal'
    elif tp in ['infant', 'toddler']:
        return 'Infancy'
    elif tp == 'school':
        return 'Childhood'
    elif tp == 'teenager':
        return 'Adolescent'
    elif tp in ['youngTeenager', 'oldTeenager']:
        return 'Adolescent'
    elif tp in ['youngAdult', 'youngMidAge', 'adult']:
        return 'Adult'
    elif tp in ['olderMidAge', 'Senior', 'senior']:
        return 'Elderly'
    return None  # This should ideally not happen

data_long['Broad_Time_Point'] = data_long['Time_Point'].apply(categorize_time_point)

# Ensure there are no missing values in Broad_Time_Point
missing_categories = data_long[data_long['Broad_Time_Point'].isnull()]
if not missing_categories.empty:
    print("Warning: Some time points could not be categorized")
    print(missing_categories)

# Aggregate data by taking the median expression value for each gene in each tissue and broad time point
data_aggregated = data_long.groupby(['ENSG', 'Tissue', 'Broad_Time_Point']).agg({'Expression_Value': 'median'}).reset_index()

# Filter the data for the list of mitochondrial genes
gene_list = pd.read_excel(gene_list_path)['Gene Ensembl ID']
data_filtered = data_aggregated[data_aggregated['ENSG'].isin(gene_list)]

# Calculate preferential values
def calculate_preferential_values(data, gene_list):
    'filter gene_list to the genes from gene list that appear in the data'
    gene_list = [gene for gene in gene_list if gene in data['ENSG'].unique()]

    preferential_values_tissue = []
    preferential_values_time = []

    # Calculate preferential values across tissues at a given time point
    for tp in data['Broad_Time_Point'].unique():
        subset_tp = data[data['Broad_Time_Point'] == tp]
        for t in subset_tp['Tissue'].unique():
            subset_tp_t = subset_tp[subset_tp['Tissue'] == t]

            for gene in gene_list:
                val_gt_list = subset_tp_t[subset_tp_t['ENSG'] == gene]['Expression_Value'].values
                if len(val_gt_list) == 0:
                    continue
                val_gt = val_gt_list[0]
                val_gt_all = subset_tp[subset_tp['ENSG'] == gene]['Expression_Value']

                median_all = np.median(val_gt_all)
                iqr_all = np.percentile(val_gt_all, 75) - np.percentile(val_gt_all, 25)

                pref_val = (val_gt - median_all) / iqr_all if iqr_all != 0 else 0
                preferential_values_tissue.append((gene, t, tp, pref_val))

    # Calculate preferential values across time points for a given tissue
    for t in data['Tissue'].unique():
        subset_t = data[data['Tissue'] == t]
        for tp in subset_t['Broad_Time_Point'].unique():
            subset_t_tp = subset_t[subset_t['Broad_Time_Point'] == tp]

            for gene in gene_list:
                val_gtp_list = subset_t_tp[subset_t_tp['ENSG'] == gene]['Expression_Value'].values
                if len(val_gtp_list) == 0:
                    continue
                val_gtp = val_gtp_list[0]
                val_gtp_all = subset_t[subset_t['ENSG'] == gene]['Expression_Value']

                median_all = np.median(val_gtp_all)
                iqr_all = np.percentile(val_gtp_all, 75) - np.percentile(val_gtp_all, 25)

                pref_val = (val_gtp - median_all) / iqr_all if iqr_all != 0 else 0
                preferential_values_time.append((gene, t, tp, pref_val))

    return pd.DataFrame(preferential_values_tissue,\
                        columns=['Gene', 'Tissue', 'Time_Point', 'Preferential_Value_Tissue']), \
           pd.DataFrame(preferential_values_time, columns=['Gene', 'Tissue', 'Time_Point', 'Preferential_Value_Time'])

# Calculate preferential values
preferential_values_tissue, preferential_values_time = calculate_preferential_values(data_filtered, gene_list)
#preferential_values_tissue.to_csv('development_preferential_values_tissue.csv', index=False)
#preferential_values_time.to_csv('development_preferential_values_time.csv', index=False)

# Load tissue associations
tissue_associations = pd.read_excel(tissue_associations_path, sheet_name='development')
tissue_associations['development'] = tissue_associations['development'].str.lower()
tissue_associations['mitoModel'] = tissue_associations['mitoModel'].str.lower()
preferential_values_tissue['Tissue'] = preferential_values_tissue['Tissue'].str.lower()
preferential_values_time['Tissue'] = preferential_values_time['Tissue'].str.lower()

# Merge preferential values with tissue associations
merged_tissue = preferential_values_tissue.merge(tissue_associations, left_on='Tissue', right_on='development')
merged_time = preferential_values_time.merge(tissue_associations, left_on='Tissue', right_on='development')

# Pivot the data
pivot_tissue = merged_tissue.pivot_table(index='Gene', columns=['mitoModel', 'Time_Point'], values='Preferential_Value_Tissue')
pivot_tissue.columns = [f"Tissue_{col[0]}_{col[1]}_prefexp" for col in pivot_tissue.columns]
pivot_tissue.reset_index(inplace=True)

pivot_time = merged_time.pivot_table(index='Gene', columns=['mitoModel', 'Time_Point'], values='Preferential_Value_Time')
pivot_time.columns = [f"Time_{col[0]}_{col[1]}_prefexp" for col in pivot_time.columns]
pivot_time.reset_index(inplace=True)

# Save the final DataFrame to CSV files
pivot_tissue.to_csv('final_tissue_preferential_values.csv', index=False)
pivot_time.to_csv('final_time_preferential_values.csv', index=False)