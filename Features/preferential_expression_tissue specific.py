import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr

# Function to calculate the model's prediction values
def calculate_predictions(algorithm, X, y, genes_list, num_folds=5):
    skf = StratifiedGroupKFold(n_splits=num_folds, shuffle=True, random_state=42)
    predictions = np.zeros(len(y))

    for train_idx, val_idx in skf.split(X, y, groups=genes_list):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        # Fit the algorithm
        algorithm.fit(X_train, y_train)
        predictions[val_idx] = algorithm.predict_proba(X_val)[:, 1]

    return predictions

def na_to_median(col):
    '''This function takes a column and fills N/A and infinity values in the column with the median value of the column'''
    col = col.replace([np.inf, -np.inf], np.nan).fillna(col.median())
    return col

# Load the dataset
dataset_path = r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confidence ds for mitocarta omim positive with pref\Most_confidence_manual_labels_with_features_dataset_new_all omim positive mitocarta_withpref.csv"
dataset = pd.read_csv(dataset_path).iloc[:, 1:]

# Extract features and target
features = dataset.iloc[:, 3:]  # Assuming features start from the 4th column
features = features.apply(na_to_median, axis=0)
target = dataset['label']  # Replace 'label' with your actual target column name
genes_list = dataset['Ensembl_ID_x']

# Define the RandomForestClassifier
model = RandomForestClassifier()

# Calculate the model's prediction values
predictions = calculate_predictions(model, features, target, genes_list)

# Add the predictions and preferential expression values to the dataset
dataset['predictions'] = predictions
dataset['preferential_expression'] = features['Preferential_Expression']

# Divide the gene-tissue pairs based on their labels
# Split the dataset by label
group_0 = dataset[dataset['label'] == 0]
group_1 = dataset[dataset['label'] == 1]

# Plot the preferential expression distribution
plt.figure(figsize=(10, 6))
sns.kdeplot(group_0['preferential_expression'], label='Label 0', shade=True)
sns.kdeplot(group_1['preferential_expression'], label='Label 1', shade=True)
plt.title('Preferential Expression Distribution by Label')
plt.xlabel('Preferential Expression')
plt.ylabel('Density')
plt.legend()
plt.xlim(-1000, 1000)
plt.ylim(0, 0.009)
plt.show()

# Conduct a Mann-Whitney U test
u_stat, p_value = mannwhitneyu(group_0['preferential_expression'], group_1['preferential_expression'])
print(f'Mann-Whitney U test results: U-statistic = {u_stat}, p-value = {p_value}')

# Calculate the correlation and create a scatter plot with a trend line
model_predictions = dataset['predictions']
preferential_expression = dataset['preferential_expression']

correlation, p_value = spearmanr(model_predictions, preferential_expression)
print(f'Spearman correlation: {correlation}, p-value: {p_value}')
dataset['log_preferential_expression'] = np.log1p(dataset['preferential_expression'])
plt.figure(figsize=(12, 8))
sns.regplot(data=dataset, x='log_preferential_expression', y='predictions', scatter_kws={'s': 10, 'alpha': 0.5})
plt.title('Correlation between Model Predictions and Preferential Expression', fontsize=16)
plt.xlabel('Log of Preferential Expression', fontsize=14)
plt.ylabel('Model Prediction', fontsize=14)
plt.grid(True)
plt.show()