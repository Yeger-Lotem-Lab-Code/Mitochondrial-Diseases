import pandas as pd
import numpy as np
from sklearn.model_selection import cross_val_score, GroupKFold
from sklearn.metrics import make_scorer, precision_recall_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import shap
import matplotlib.pyplot as plt
from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.inspection import permutation_importance, PartialDependenceDisplay


# Load your dataset
df = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\MITOCH~1\MANUAL~1\MOSTCO~4\MOSTCO~2.CSV").iloc[:, 1:]  # Uncomment and modify this line to load your dataset

for col in df.columns[4:]:
    df[col].fillna(df[col].median(), inplace=True)
for col in df.columns[4:]:
    if np.isinf(df[col]).any():
        max_non_inf = df.loc[~np.isinf(df[col]), col].max()
        df[col].replace(np.inf, max_non_inf, inplace=True)

# Assuming 'label' is your target column and 'Ensembl_ID_x' is the gene identifier
X = df.drop(['label', 'Ensembl_ID_x', 'Gene Name', 'Tissue'], axis=1)
y = df['label']
genes_list = df['Ensembl_ID_x']

# Scale the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Calculate VIF for each feature and remove those with high VIF
vif_data = pd.DataFrame()
vif_data["feature"] = X.columns
vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(len(X.columns))]
high_vif_features = vif_data[vif_data["VIF"] > 10]["feature"]
X_vif = X.drop(high_vif_features, axis=1)
X_scaled_vif = scaler.fit_transform(X_vif)

# Define the RandomForest model with class weights
best_model = RandomForestClassifier(random_state=42, class_weight='balanced')

# Define evaluation metrics
scoring = ['accuracy', 'precision', 'recall', 'f1']

# Define cross-validation strategy with stratification by genes_list
cv = GroupKFold(n_splits=5)

# Function to calculate precision-recall AUC
def pr_auc_score(y_true, y_pred):
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    return auc(recall, precision)

pr_auc = make_scorer(pr_auc_score, needs_proba=True)

# Evaluate the model
results = {}
scores = cross_val_score(best_model, X_scaled_vif, y, cv=cv, scoring='accuracy', groups=genes_list)
pr_auc_scores = cross_val_score(best_model, X_scaled_vif, y, cv=cv, scoring=pr_auc, groups=genes_list)
mean_pr_auc = np.mean(pr_auc_scores)
results['Random Forest'] = {
    'accuracy': np.mean(scores),
    'pr_auc': mean_pr_auc
}

# Display results
for model_name, metrics in results.items():
    print(f"{model_name}:")
    print(f"  Accuracy: {metrics['accuracy']:.4f}")
    print(f"  PR AUC: {metrics['pr_auc']:.4f}")

# Fit the model on the entire dataset
best_model.fit(X_scaled_vif, y)

# Calculate SHAP values
explainer = shap.Explainer(best_model, X_scaled_vif)
shap_values = explainer(X_scaled_vif)

# Plot SHAP summary
shap.summary_plot(shap_values, X_scaled_vif, feature_names=X_vif.columns)
plt.title('SHAP Summary Plot for Random Forest (Scaled and Reduced Features)')
plt.show()

# Plot Permutation Importances
result = permutation_importance(best_model, X_scaled_vif, y, n_repeats=10, random_state=42, n_jobs=2)
sorted_idx = result.importances_mean.argsort()
plt.boxplot(result.importances[sorted_idx].T, vert=False, labels=X_vif.columns[sorted_idx])
plt.title("Permutation Importances (train set)")
plt.tight_layout()
plt.show()

# Plot Partial Dependence Plots
features = list(range(X_vif.shape[1]))
PartialDependenceDisplay.from_estimator(best_model, X_scaled_vif, features, grid_resolution=50)
plt.show()