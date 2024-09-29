import os
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier
from sklearn.metrics import auc, precision_recall_curve, roc_curve, average_precision_score
import matplotlib.pyplot as plt
import numpy as np

def calculate_roc_pr_auc(algorithm, X, y, num_folds=5):
    skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=42)
    roc_curve_data = []
    pr_curve_data = []
    top_misclassified_genes = []

    for train_idx, val_idx in skf.split(X, y):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        # Fit the algorithm
        algorithm.fit(X_train, y_train)
        y_pred_proba = algorithm.predict_proba(X_val)[:, 1]

        # Calculate ROC curve and AUC
        fpr, tpr, _ = roc_curve(y_val, y_pred_proba)
        roc_auc = auc(fpr, tpr)
        roc_curve_data.append((fpr, tpr, roc_auc))

        # Calculate precision-recall curve and AUC
        precision, recall, _ = precision_recall_curve(y_val, y_pred_proba)
        pr_auc = average_precision_score(y_val, y_pred_proba)
        pr_curve_data.append((precision, recall, pr_auc))

        # Store misclassified genes
        val_dataset = dataset.iloc[val_idx].copy()
        val_dataset['Predicted_Probability'] = y_pred_proba
        val_dataset['True_Label'] = y_val

        # Filter where true label is 0 but predicted probability is high (>0.5)
        misclassified = val_dataset[(val_dataset['True_Label'] == 0) & (val_dataset['Predicted_Probability'] > 0.5)]
        top_misclassified_genes.append(misclassified)

    # Concatenate all folds' misclassified data and sort by predicted probability
    top_misclassified_genes = pd.concat(top_misclassified_genes)
    top_misclassified_genes = top_misclassified_genes.sort_values(by='Predicted_Probability', ascending=False).head(10)

    return roc_curve_data, pr_curve_data, top_misclassified_genes


# Algorithm configuration
algorithms = {
    "XGBoost_tuned": XGBClassifier(
                                   colsample_bytree=np.float64(0.898027900926573),
                                   learning_rate=np.float64(0.02738831867639583),
                                   max_depth=6, n_estimators=175, subsample=np.float64(0.9315111447751979))
}

# Dataset loading
dataset_path = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv"
dataset = pd.read_csv(dataset_path)
features = dataset.iloc[:, 3:]  # Adjust the range if necessary
target = dataset.iloc[:, 2]  # Assuming the third column is the target
dataset_results = {algorithm_name: {"roc_curve": [], "pr_curve": [], "top_misclassified_genes": []} for algorithm_name in algorithms}

# Model execution and result collection
for algorithm_name, algorithm in algorithms.items():

    # Calculate ROC, PR AUC, and misclassified genes
    roc_curve_data, pr_curve_data, top_misclassified_genes = calculate_roc_pr_auc(algorithm, features, target, num_folds=5)

    # Store results
    dataset_results[algorithm_name]["roc_curve"] = roc_curve_data
    dataset_results[algorithm_name]["pr_curve"] = pr_curve_data
    dataset_results[algorithm_name]["top_misclassified_genes"] = top_misclassified_genes

    # Plot ROC and PR curves for all folds
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f'Cross-Validation Results for {algorithm_name}', fontsize=16)

    # Plot ROC curve
    for i, (fpr, tpr, roc_auc) in enumerate(roc_curve_data):
        ax1.plot(fpr, tpr, label=f'Fold {i+1} (AUC = {roc_auc:.2f})')
    ax1.plot([0, 1], [0, 1], linestyle='--', color='black')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_title('ROC Curve')
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.legend(loc='lower right')
    ax1.grid(True)

    # Plot PR curve
    for i, (precision, recall, pr_auc) in enumerate(pr_curve_data):
        ax2.plot(recall, precision, label=f'Fold {i + 1} (PR AUC = {pr_auc:.2f})')

    # Calculate expected precision-recall for random
    total_positive_instances = sum(target)
    total_instances = len(target)
    expected_prc = total_positive_instances / total_instances

    # Add dashed line for expected precision-recall
    ax2.axhline(y=expected_prc, color='black', linestyle='--', label=f'Expected PRC = {expected_prc:.2f}')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_title('Precision-Recall Curve')
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.legend(loc='lower right')
    ax2.grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save the figure
    plt.show()
    plt.close()

    # Display the top 10 misclassified genes
    print(
        f"Top 10 misclassified genes for {algorithm_name} (labeled as negative but predicted as high probability of being causal):")
    print(top_misclassified_genes[['Gene_x', 'Predicted_Probability']])

    # Plot ROC curves for each algorithm for this dataset
plt.figure(figsize=(10, 8))
for algorithm_name in algorithms:
    mean_fpr = np.linspace(0, 1, 100)
    mean_tpr = np.mean([np.interp(mean_fpr, fpr, tpr) for fpr, tpr, _ in dataset_results[algorithm_name]["roc_curve"]],
                       axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, label=f'{algorithm_name} (AUC = {mean_auc:.2f})')

plt.plot([0, 1], [0, 1], linestyle='--', color='black')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()
plt.close()

# Create box plots for AUC scores for this dataset
auc_scores = [dataset_results[algorithm_name]["roc_curve"] for algorithm_name in algorithms]
auc_scores = [[roc_auc for _, _, roc_auc in roc_curve] for roc_curve in auc_scores]

plt.figure(figsize=(10, 6))
algorithm_names = list(algorithms.keys())
plt.boxplot(auc_scores, labels=algorithm_names)
plt.title('AUC Scores Comparison')
plt.ylim(0, 1)
plt.xlabel('Algorithms')
plt.ylabel('AUC Score')
plt.grid(True)
plt.show()
plt.close()

# Show results
for algorithm_name in algorithms:
    print(f"Results for {algorithm_name}:")
    print("AUC scores:", [roc_auc for _, _, roc_auc in dataset_results[algorithm_name]["roc_curve"]])
    print("Precision-Recall scores:", [pr_auc for _, _, pr_auc in dataset_results[algorithm_name]["pr_curve"]])
    print(f"Top 10 misclassified genes (labeled as negative but predicted as causal) for {algorithm_name}:")
    print(dataset_results[algorithm_name]["top_misclassified_genes"][['Gene_x', 'Predicted_Probability']])