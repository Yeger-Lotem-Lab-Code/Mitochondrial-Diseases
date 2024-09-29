import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold, cross_val_score
import matplotlib.pyplot as plt
import numpy as np

# Load your dataset
# Make sure to replace 'your_dataset.csv' with the path to your dataset
df = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Causality in tissues\Manual Curation\Most confidence ds for mitocarta omim positive with pref\Most confident tissue-specific dataset with only positive tissues_new_all omim positive mitocarta with pref.csv").iloc[:,1:]
df.replace([np.inf, -np.inf], np.nan, inplace=True)
# Split dataset into features and target variable
X = df.iloc[:,4:]
X = X.apply(lambda x: x.fillna(x.median()))
y = df['label']
groups = df['Ensembl_ID_x']
group_kfold = GroupKFold(n_splits=5)

# Create a Random Forest classifier
rf = RandomForestClassifier(random_state=42)

auc_scores = cross_val_score(rf, X, y, groups=groups, cv=group_kfold, scoring='roc_auc')
pr_scores = cross_val_score(rf, X, y, groups=groups, cv=group_kfold, scoring='average_precision')

print(f"Average AUC: {auc_scores.mean():.2f}")
print(f"Average PR: {pr_scores.mean():.2f}")

rf.fit(X, y)
# Get feature importances
feature_importances = rf.feature_importances_

# Create a DataFrame for the feature importances
features = pd.DataFrame({
    'Feature': X.columns,
    'Importance': feature_importances
})
# Sort features by importance
features = features.sort_values(by='Importance', ascending=False)

# Plot feature importances
plt.figure(figsize=(30, 30))
plt.barh(features['Feature'], features['Importance'])
plt.xlabel('Importance')
plt.ylabel('Feature')
plt.title('Feature Importances from Random Forest')
plt.gca().invert_yaxis()
plt.show()

N = [105, 106,107,108,109,110] # Number of top features to select
for n in N:
    top_features = features.head(n)['Feature']

# Create a new dataset with only the top N features
    X_top = X[top_features]

    auc_scores_top = cross_val_score(rf, X_top, y, groups=groups, cv=group_kfold, scoring='roc_auc')
    pr_scores_top = cross_val_score(rf, X_top, y, groups=groups, cv=group_kfold, scoring='average_precision')

    print(f"Average AUC with top {n} features: {auc_scores_top.mean():.2f}")
    print(f"Average PR with top {n} features: {pr_scores_top.mean():.2f}")