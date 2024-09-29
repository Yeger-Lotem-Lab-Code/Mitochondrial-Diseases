import pandas as pd
from xgboost import XGBClassifier
import shap
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.preprocessing import MinMaxScaler


def infer_feature_type(name):
    # Convert feature name to lowercase for consistent comparisons
    name = name.lower()

    # Initialize a variable to store the inferred feature type
    new_name = 'gene specific'

    # Define conditions for matching feature types
    if 'preferential_expression' in name:
        new_name = 'Tissue Preferential Expression'
    elif 'preferential_num_interactors' in name or 'preferential_num_elevated_interactors' in name or 'preferential_num_specific_interactions' in name or 'diff_net_med' in name:
        new_name = 'Tissue PPIs (Preferential)'
    elif 'proact' in name:
        new_name = 'Tissue ProACT Scores'
    elif '(neutral) median crispr score' in name or '(neutral) median rnai score' in name:
        new_name = 'Tissue Essentiality Scores'
    elif 'peak_intensity' in name:
        new_name = 'Tissue Peak Intensity'
    elif 'egene' in name:
        new_name = 'Tissue eQTL'
    elif 'paralogs_ratio_highest_identity' in name:
        new_name = 'Tissue Paralogs Relationships'
    elif 'ts score' in name:
        new_name = 'Tissue TS Score'
    elif 'age coefficient' in name:
        new_name = 'Tissue Age Coefficient'
    elif 'time' in name:
        new_name = 'Tissue Development Time'
    elif 'tissue' in name:
        new_name = 'Tissue Development Tissue'
    elif 'development_cv' in name:
        new_name = 'Tissue Development CV'
    # Return the inferred feature type
    return new_name

def infer_origin_tissue(name):

    name = name.lower()

    new_name = 'gene specific'

    if 'heart1' in name:
        new_name = 'Heart1'
    elif 'heart2' in name:
        new_name = 'Heart2'
    elif 'whole heart' in name:
        new_name = 'Whole Heart'
    elif 'brain1' in name:
        new_name = 'Brain1'
    elif 'brain0' in name:
        new_name = 'Brain0'
    elif 'brain2' in name:
        new_name = 'Brain2'
    elif 'brain3' in name:
        new_name = 'Brain3'
    elif 'brain4' in name:
        new_name = 'Brain4'
    elif 'brain5' in name:
        new_name = 'Brain5'
    elif 'whole brain' in name:
        new_name = 'Whole Brain'
    elif 'kidney1' in name:
        new_name = 'Kidney1'
    elif 'kidney2' in name:
        new_name = 'Kidney2'
    elif 'whole kidney' in name:
        new_name = 'Whole Kidney'
    elif 'esophagus1' in name:
        new_name = 'Esophagus1'
    elif 'esophagus2' in name:
        new_name = 'Esophagus2'
    elif 'esophagus3' in name:
        new_name = 'Esophagus3'
    elif 'whole esophagus' in name:
        new_name = 'Whole Esophagus'
    elif 'spleen' in name:
        new_name = 'Spleen'
    elif 'pituitary' in name:
        new_name = 'Pituitary'
    elif 'adrenal' in name:
        new_name = 'Adrenal Gland'
    elif 'skin1' in name:
        new_name = 'Skin1'
    elif 'skin2' in name:
        new_name = 'Skin2'
    elif 'whole skin' in name:
        new_name = 'Whole Skin'
    elif 'uterus' in name:
        new_name = 'Uterus'
    elif 'artery1' in name:
        new_name = 'Artery1'
    elif 'artery2' in name:
        new_name = 'Artery2'
    elif 'artery3' in name:
        new_name = 'Artery3'
    elif 'whole artery' in name:
        new_name = 'Whole Artery'
    elif 'muscle' in name:
        new_name = 'Muscle'
    elif 'stomach' in name:
        new_name = 'Stomach'
    elif 'small intestine' in name:
        new_name = 'Small Intestine'
    elif 'pancreas' in name:
        new_name = 'Pancreas'
    elif 'thyroid' in name:
        new_name = 'Thyroid'
    elif 'vagina' in name:
        new_name = 'Vagina'
    elif 'ovary' in name:
        new_name = 'Ovary'
    elif 'adipose1' in name:
        new_name = 'Adipose1'
    elif 'adipose2' in name:
        new_name = 'Adipose2'
    elif 'whole adipose' in name:
        new_name = 'Whole Adipose'
    elif 'liver' in name:
        new_name = 'Liver'
    elif 'cervix1' in name:
        new_name = 'Cervix1'
    elif 'cervix2' in name:
        new_name = 'Cervix2'
    elif 'whole cervix' in name:
        new_name = 'Whole Cervix'
    elif 'nerve' in name:
        new_name = 'Nerve'
    elif 'testis' in name:
        new_name = 'Testis'
    elif 'blood' in name:
        new_name = 'Blood'
    elif 'salivary' in name:
        new_name = 'Minor Salivary Gland'
    elif 'colon1' in name:
        new_name = 'Colon1'
    elif 'colon2' in name:
        new_name = 'Colon2'
    elif 'whole colon' in name:
        new_name = 'Whole Colon'
    elif 'lung' in name:
        new_name = 'Lung'
    elif 'fallopian' in name:
        new_name = 'Fallopian Tube'
    elif 'fibroblasts' in name:
        new_name = 'Fibroblasts'
    elif 'breast' in name:
        new_name = 'Breast'
    elif 'bladder' in name:
        new_name = 'Bladder'
    elif 'prostate' in name:
        new_name = 'Prostate'

    return new_name

dataset_path = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv"
Dataset = pd.read_csv(dataset_path)
X = Dataset.iloc[:, 3:]
y = Dataset.iloc[:,2]

model = XGBClassifier(use_label_encoder=False, eval_metric='logloss',colsample_bytree = np.float64(0.898027900926573), learning_rate= np.float64(0.02738831867639583), max_depth= 6, n_estimators= 175, subsample= np.float64(0.9315111447751979))
model.fit(X, y)
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)
shap_df = pd.DataFrame(shap_values, columns=X.columns)
# Create mappings for feature types and origin tissues based on column names
feature_types = {col: infer_feature_type(col) for col in shap_df.columns}
origin_tissues = {col: infer_origin_tissue(col) for col in shap_df.columns}
shap_sum = np.abs(shap_values).mean(axis=0)

d = {'column_name': X.columns.tolist(), 'shap_importance': shap_sum.tolist()}
importance_df = pd.DataFrame(d)
importance_df = importance_df.sort_values('shap_importance', ascending=False)
for row_name in importance_df.index.values:
    feature_name = importance_df.loc[row_name, 'column_name']
    feature_type = infer_feature_type(feature_name)
    origin_tissue = infer_origin_tissue(feature_name)
    importance_df.loc[row_name, 'feature_type'] = feature_type
    importance_df.loc[row_name, 'origin_tissue'] = origin_tissue

shap_sum = importance_df.loc[:, 'shap_importance'].sum(axis=0)

feature_type_importance = {}
origin_tissue_importance = {}

for feature in importance_df['feature_type'].unique():
    feature_sum = importance_df.loc[importance_df['feature_type'] == feature, 'shap_importance'].sum()
    feature_type_importance[feature] = feature_sum / shap_sum

for tissue in importance_df['origin_tissue'].unique():
    tissue_sum = importance_df.loc[importance_df['origin_tissue'] == tissue, 'shap_importance'].sum()
    origin_tissue_importance[tissue] = tissue_sum / shap_sum

# Convert the dictionaries to DataFrames for easier plotting
feature_type_importance_df = pd.DataFrame(list(feature_type_importance.items()), columns=['Feature Type', 'Normalized SHAP Importance'])
feature_type_importance_df = feature_type_importance_df.iloc[1:,:]
origin_tissue_importance_df = pd.DataFrame(list(origin_tissue_importance.items()), columns=['Origin Tissue', 'Normalized SHAP Importance'])
origin_tissue_importance_df = origin_tissue_importance_df.iloc[1:,:]

# Plot the normalized SHAP importance by feature type
plt.figure(figsize=(12, 8))
sns.barplot(x='Normalized SHAP Importance', y='Feature Type', data=feature_type_importance_df.sort_values(by='Normalized SHAP Importance', ascending=False))
plt.title('Normalized SHAP Importance by Feature Type')
plt.xlabel('Normalized SHAP Importance')
plt.ylabel('Feature Type')
plt.tight_layout()
plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\Aggregation Analysis\Aggregated_SHAP_Feature_Type.png")
plt.show()

# Plot the normalized SHAP importance by origin tissue
plt.figure(figsize=(12, 8))
sns.barplot(x='Normalized SHAP Importance', y='Origin Tissue', data=origin_tissue_importance_df.sort_values(by='Normalized SHAP Importance', ascending=False))
plt.title('Normalized SHAP Importance by Origin Tissue')
plt.xlabel('Normalized SHAP Importance')
plt.ylabel('Origin Tissue')
plt.tight_layout()
plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\Aggregation Analysis\Aggregated_SHAP_Tissue.png")
plt.show()



