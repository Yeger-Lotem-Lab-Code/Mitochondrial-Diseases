import pandas as pd
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset
df = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_pref.csv")

# Replace inf values with NaN and fill NaN with the median of each column
df.replace([np.inf, -np.inf], np.nan, inplace=True)
X = df.iloc[:, 3:].apply(lambda x: x.fillna(x.median()))

# Target variable
y = df['Label_OMIM']

# Initialize the XGBoost classifier
xgb = XGBClassifier()

# Perform cross-validation to get AUC scores
auc_scores = cross_val_score(xgb, X, y, scoring='roc_auc', cv=5)
pr_scores = cross_val_score(xgb, X, y, scoring='average_precision', cv=5)

# Print average AUC and PR scores
print(f"Average AUC: {auc_scores.mean():.2f}")
print(f"Average PR: {pr_scores.mean():.2f}")

# Fit the model
xgb.fit(X, y)

# Get feature importances
feature_importances = xgb.feature_importances_

# Create a DataFrame for features and their importances
features = pd.DataFrame({
    'Feature': X.columns,
    'Importance': feature_importances
})

# Sort features by importance
features = features.sort_values(by='Importance', ascending=False)

# Plot feature importances
'''
plt.figure(figsize=(10, 30))  # Adjust the size as needed
plt.barh(features['Feature'], features['Importance'])
plt.xlabel('Importance', fontsize=14)
plt.ylabel('Feature', fontsize=14)
plt.title('Feature Importances from XGBoost', fontsize=14)
plt.gca().invert_yaxis()
plt.xticks(fontsize=12)
plt.yticks(fontsize=6)
plt.show()
'''
# Define the number of top features to select
N = list(range(100, 400, 10))  # A range from 130 to 200 with a step of 10

auc_results = []
pr_results = []

for n in N:
    # Select the top n features
    top_features = features.head(n)['Feature']

    # Create a new dataset with only the top n features
    X_top = X[top_features]

    # Perform cross-validation to get AUC and PR scores for the top n features
    auc_scores_top = cross_val_score(xgb, X_top, y, scoring='roc_auc', cv=5)
    pr_scores_top = cross_val_score(xgb, X_top, y, scoring='average_precision', cv=5)

    # Store the results
    auc_results.append(auc_scores_top.mean())
    pr_results.append(pr_scores_top.mean())

    # Print average AUC and PR scores for the top n features
    print(f"Average AUC with top {n} features: {auc_scores_top.mean():.2f}")
    print(f"Average PR with top {n} features: {pr_scores_top.mean():.2f}")

# Plotting the results
plt.figure(figsize=(14, 7))

plt.subplot(1, 2, 1)
plt.plot(N, auc_results, marker='o')
plt.xlabel('Number of Features', fontsize=14)
plt.ylabel('AUC', fontsize=14)
plt.title('AUC vs Number of Features', fontsize=14)

plt.subplot(1, 2, 2)
plt.plot(N, pr_results, marker='o')
plt.xlabel('Number of Features', fontsize=14)
plt.ylabel('PR', fontsize=14)
plt.title('PR vs Number of Features', fontsize=14)

plt.tight_layout()
plt.show()