import os
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
import shap
import numpy as np
import matplotlib.pyplot as plt

# Function to extract top misclassified genes and perform SHAP analysis
def find_top_misclassified_pairs(algorithm, genes_list, X, y, dataset, num_folds=5):
    skf = StratifiedGroupKFold(n_splits=num_folds, shuffle=True, random_state=42)
    misclassified_pairs_all_folds = []

    # Iterate over each fold
    for fold_num, (train_idx, val_idx) in enumerate(skf.split(X, y, groups=genes_list)):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        # Fit the model
        algorithm.fit(X_train, y_train)

        # Predict probabilities on the validation set
        y_pred_proba = algorithm.predict_proba(X_val)[:, 1]

        # Find misclassified genes
        val_dataset = dataset.iloc[val_idx].copy()
        val_dataset['Predicted_Probability'] = y_pred_proba
        val_dataset['True_Label'] = y_val.values
        misclassified = val_dataset[
            (val_dataset['True_Label'] == 0) & (val_dataset['Predicted_Probability'] > 0.5)].copy()
        misclassified['Fold'] = fold_num

        # Append misclassified genes from each fold
        misclassified_pairs_all_folds.append(misclassified)

    # Concatenate all misclassified genes from all folds
    if misclassified_pairs_all_folds:
        misclassified_pairs_all_folds = pd.concat(misclassified_pairs_all_folds)
        # Sort by predicted probability and select the top 3 misclassified genes
        top_misclassified_pairs = misclassified_pairs_all_folds.sort_values(by='Predicted_Probability',
                                                                            ascending=False).head(3)
    else:
        top_misclassified_pairs = pd.DataFrame()

    return misclassified_pairs_all_folds, top_misclassified_pairs



def generate_shap_plots_per_pair(algorithm, X, y, top_misclassified_pairs, results_dir):
    # Loop through the top 3 misclassified genes
    for idx in top_misclassified_pairs.index:
        gene_name = top_misclassified_pairs.loc[idx, 'Gene Name']
        tissue = top_misclassified_pairs.loc[idx, 'Tissue']
        # Create train and test sets: test set is the specific gene, train set is all others
        X_test = X.loc[[idx]]
        X_train = X.drop(idx)
        y_train = y.drop(idx)

        # Fit model on the training set excluding the specific gene
        model = algorithm.fit(X_train, y_train)

        # Initialize SHAP Explainer for the fitted model
        explainer = shap.TreeExplainer(model)
        test_shap_values = explainer.shap_values(X_test)
        expected_value = explainer.expected_value

        # SHAP Summary Plot for the specific gene
        plt.figure(figsize=(10, 8))
        shap.summary_plot(test_shap_values[:,:,1], X_test, max_display=20, show=False)
        plt.savefig(f"{results_dir}/shap_summary_plot_gene_{gene_name},{tissue}.png", dpi=150, bbox_inches='tight', pad_inches=0)
        plt.close()

        # SHAP Decision Plot for the specific gene
        plt.figure(figsize=(10, 8))
        shap.decision_plot(expected_value[1], test_shap_values[:,:,1], X_test, show=False)
        plt.savefig(f"{results_dir}/shap_decision_plot_gene_{gene_name},{tissue}.png", dpi=150, bbox_inches='tight',
                    pad_inches=0)
        plt.close()




def main():

    # Dataset loading
    dataset_path = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\MOSTCO~3_imputed_dataset.csv"
    dataset = pd.read_csv(dataset_path).iloc[:, 1:]
    genes_list = dataset['Ensembl_ID_x']  # Adjust this to the column containing gene identifiers
    features = dataset.iloc[:, 4:]  # Adjust based on your actual data structure
    target = dataset['label']  # Replace 'label' with your actual target column name
    algorithm = RandomForestClassifier(
        n_estimators=50,
        min_samples_split=2,
        min_samples_leaf=2,
        max_features='sqrt',
        max_depth=10,
        bootstrap=False
    )
    # Extract the top misclassified genes and perform SHAP analysis consistently
    misclassified_pairs_all_folds, top_misclassified_pairs = find_top_misclassified_pairs(algorithm, genes_list, features, target, dataset, num_folds=5)
    if not top_misclassified_pairs.empty:
        # Generate SHAP plots for the top 3 misclassified genes
        generate_shap_plots_per_pair(algorithm, features, target, top_misclassified_pairs,
                                     results_dir=r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\Proof of Concept")

        print("All Misclassified pairs Across Folds:")
        print(misclassified_pairs_all_folds[['Gene Name', 'Tissue', 'Predicted_Probability', 'Fold']])

        print("\nTop 3 Misclassified Genes:")
        print(top_misclassified_pairs[['Gene Name', 'Tissue', 'Predicted_Probability']])
    else:
        print("No misclassified genes found.")

if __name__ == "__main__":
    main()