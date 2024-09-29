import os
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier


def find_top_misclassified_genes(algorithm, X, y, dataset, num_folds=5):
    skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=42)
    misclassified_genes_all_folds = []

    # Iterate over each fold
    for fold_num, (train_idx, val_idx) in enumerate(skf.split(X, y)):
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
        misclassified_genes_all_folds.append(misclassified)

    # Concatenate all misclassified genes from all folds
    if misclassified_genes_all_folds:
        misclassified_genes_all_folds = pd.concat(misclassified_genes_all_folds)
        # Sort by predicted probability and select the top 3 misclassified genes
        top_misclassified_genes = misclassified_genes_all_folds.sort_values(by='Predicted_Probability',
                                                                            ascending=False).head(3)
    else:
        top_misclassified_genes = pd.DataFrame()

    return misclassified_genes_all_folds, top_misclassified_genes


def generate_shap_plots_per_gene(algorithm, X, y, top_misclassified_genes, results_dir):
    # Loop through the top 3 misclassified genes
    for idx in top_misclassified_genes.index:
        gene_name = top_misclassified_genes.loc[idx, 'Gene_x']

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
        shap.summary_plot(test_shap_values, X_test, max_display=20, show=False)
        plt.savefig(f"{results_dir}/shap_summary_plot_gene_{gene_name}.png", dpi=150, bbox_inches='tight', pad_inches=0)
        plt.close()

        # SHAP Decision Plot for the specific gene
        plt.figure(figsize=(10, 8))
        shap.decision_plot(expected_value, test_shap_values, X_test, show=False)
        plt.savefig(f"{results_dir}/shap_decision_plot_gene_{gene_name}.png", dpi=150, bbox_inches='tight',
                    pad_inches=0)
        plt.close()


# Main Execution
def main():
    # Algorithm configuration
    algorithm = XGBClassifier(colsample_bytree=0.898, learning_rate=0.027, max_depth=6, n_estimators=175,
                              subsample=0.932, random_state=42)

    # Dataset loading
    dataset_path = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv"
    dataset = pd.read_csv(dataset_path)

    features = dataset.iloc[:, 3:]  # Adjust the range if necessary
    target = dataset.iloc[:, 2]  # Assuming the third column is the target

    # Identify top misclassified genes
    misclassified_genes_all_folds, top_misclassified_genes = find_top_misclassified_genes(algorithm, features, target,
                                                                                          dataset, num_folds=5)

    if not top_misclassified_genes.empty:
        # Generate SHAP plots for the top 3 misclassified genes
        generate_shap_plots_per_gene(algorithm, features, target, top_misclassified_genes,
                                     results_dir=r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\Proof of Concept")

        print("All Misclassified Genes Across Folds:")
        print(misclassified_genes_all_folds[['Gene_x', 'Predicted_Probability', 'Fold']])

        print("\nTop 3 Misclassified Genes:")
        print(top_misclassified_genes[['Gene_x', 'Predicted_Probability']])
    else:
        print("No misclassified genes found.")


if __name__ == "__main__":
    main()