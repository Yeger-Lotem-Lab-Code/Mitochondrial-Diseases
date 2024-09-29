import pandas as pd
import shap
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from xgboost import XGBClassifier

# Step 1: Read the dataset
data = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv")

# Separate features and labels
X = data.iloc[:, 3:]  # Features start from the fourth column
y = data.iloc[:, 2]   # Labels are in the third column

# Step 2: Train the model and calculate SHAP values
model = XGBClassifier(eval_metric='logloss',
                      colsample_bytree=np.float64(0.898027900926573),
                      learning_rate=np.float64(0.02738831867639583),
                      max_depth=6, n_estimators=175,
                      subsample=np.float64(0.9315111447751979))
model.fit(X, y)

explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)
shap_values_df = pd.DataFrame(shap_values, columns=X.columns)

# Step 3: Normalize SHAP values for the top most contributing features
top_features = shap_values_df.abs().mean().nlargest(10).index
xgb_shap = shap_values_df[top_features].values

# Normalizing SHAP values: subtract min, divide by max, subtract median
xgb_shap_normed = xgb_shap.copy()
for i in range(xgb_shap_normed.shape[1]):
    xgb_shap_normed[:, i] -= xgb_shap_normed[:, i].min()
    xgb_shap_normed[:, i] /= xgb_shap_normed[:, i].max()
    xgb_shap_normed[:, i] -= np.median(xgb_shap_normed[:, i])

# Step 4: Compute distance matrix using squared Euclidean distance for hierarchical clustering
D = pdist(xgb_shap_normed, metric='sqeuclidean')
linkage_matrix = linkage(D, method='complete')

# Step 5: Perform clustering and plot dendrogram with heatmap
fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(18, 10), gridspec_kw={'width_ratios': [1, 4]})

# Plot dendrogram
dendro = dendrogram(linkage_matrix, no_labels=True, orientation='left', ax=ax_dendro, color_threshold=2)
ax_dendro.set_title('Dendrogram of Hierarchical Clustering', fontsize=14)
ax_dendro.set_xlabel('Distance', fontsize=12)
ax_dendro.axhline(y=2, color='r', linestyle='--')
# Get the order of the genes from the dendrogram
ordered_indices = dendro['leaves']

# Plot heatmap
im = ax_heatmap.imshow(xgb_shap_normed[ordered_indices, :], aspect='auto', cmap='coolwarm', vmin=-1, vmax=1)
ax_heatmap.set_xticks(np.arange(len(top_features)))
ax_heatmap.set_xticklabels(top_features, rotation=45, ha='right', fontsize=12)
ax_heatmap.set_title('Heatmap of Normalized SHAP Values', fontsize=14)
ax_heatmap.spines['right'].set_visible(False)
ax_heatmap.spines['left'].set_visible(False)
ax_heatmap.spines['top'].set_visible(False)
ax_heatmap.spines['bottom'].set_visible(False)

# Add color bar for better interpretation
cbar = plt.colorbar(im, ax=ax_heatmap, label='Normalized SHAP Values')
cbar.ax.tick_params(labelsize=10)

plt.tight_layout()
plt.show()

# Step 6: Find clusters and extract genes for a specific range
# Divide the data into clusters using the threshold (e.g., 3 clusters)
clusters = fcluster(linkage_matrix, t=5, criterion='maxclust')
data['Cluster'] = clusters

for cluster_num in range(1, 6):
    cluster_genes = data[data['Cluster'] == cluster_num]['Gene_x']  # Adjust column name if different
    print(f"Cluster {cluster_num} genes:")
    for gene in cluster_genes:
        print(gene)
    print('number of genes in cluster:', len(cluster_genes))
    print("\n" + "-"*50 + "\n")
