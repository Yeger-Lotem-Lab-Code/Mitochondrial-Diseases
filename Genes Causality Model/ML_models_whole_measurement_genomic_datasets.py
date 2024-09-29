import os
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve, roc_curve, average_precision_score
import matplotlib.pyplot as plt
import numpy as np

# Define a function to calculate ROC and PR curves for an algorithm
def calculate_roc_pr_auc(algorithm, X, y, num_folds=5):
    skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=42)
    roc_curve_data = []
    pr_curve_data = []

    for train_idx, val_idx in skf.split(X, y):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

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

    return roc_curve_data, pr_curve_data

# Define algorithms
algorithms = {
    "Logistic Regression": LogisticRegression(max_iter=2000),
    "XGBoost": XGBClassifier(),
    "Gradient Boosting": GradientBoostingClassifier(),
    "Random Forest": RandomForestClassifier(),
    "XGBoost_tuned": XGBClassifier(use_label_encoder=False, eval_metric='logloss',colsample_bytree = np.float64(0.898027900926573), learning_rate= np.float64(0.02738831867639583), max_depth= 6, n_estimators= 175, subsample= np.float64(0.9315111447751979))
}

# Define the folder structure
base_folder = r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality"
label_folders = ['OMIM Labeled']
mid_labels_folder = ['_OMIM']
gene_lists_names = ['MitoCarta']
folder_to_save = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\5 fold cross validation - ds with pref and gorilla"

# Iterate through the directories and read datasets
results = {}
'''
for label_folder, mid_label in zip(label_folders, mid_labels_folder):
    for gene_list_name in gene_lists_names:
        folder_path = os.path.join(base_folder, label_folder, gene_list_name + mid_label)
        for filename in os.listdir(folder_path):
            if filename.endswith(".csv") and "imputed" in filename:
'''
dataset_path = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv"
dataset = pd.read_csv(dataset_path)
#filter dataset without Stomach_TS score column
#dataset = dataset.drop(columns=['Stomach_TS score'])
features = dataset.iloc[:, 3:]  # Adjust the range if necessary
target = dataset.iloc[:, 2]  # Assuming the third column is the target
dataset_results = {algorithm_name: {"roc_curve": [], "pr_curve": []} for algorithm_name in algorithms}

                # Iterate over algorithms
for algorithm_name, algorithm in algorithms.items():

    roc_curve_data, pr_curve_data = calculate_roc_pr_auc(algorithm, features, target, num_folds=5)

                    # Store curve data for each algorithm
    dataset_results[algorithm_name]["roc_curve"] = roc_curve_data
    dataset_results[algorithm_name]["pr_curve"] = pr_curve_data

                    # Plot ROC and PR curves for all folds on the same graph
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f'Cross-Validation Results for {algorithm_name}', fontsize=16)

                    # Plot ROC curve for all folds
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

    # Plot PR curve for all folds
    for i, (precision, recall, pr_auc) in enumerate(pr_curve_data):
        ax2.plot(recall, precision, label=f'Fold {i + 1} (PR AUC = {pr_auc:.2f})')

    # Calculate expected precision-recall for random
    total_positive_instances = sum(target)
    total_instances = len(target)
    expected_prc = total_positive_instances / total_instances

    # Add the dashed line for expected precision-recall on the Precision-Recall curve
    ax2.axhline(y=expected_prc, color='black', linestyle='--',
                label=f'Expected PRC = {expected_prc:.2f}')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_title('Precision-Recall Curve')
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.legend(loc='lower right')
    ax2.grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save the figure
    plt.savefig(folder_to_save + f'/{algorithm_name}_ROC_PR.png')
    plt.show()
    plt.close()

    # Calculate the expected value of PRC for this dataset
    total_positive_instances = sum(target)
    total_instances = len(target)
    expected_prc = total_positive_instances / total_instances

    print(f"Expected value of PRC: {expected_prc:.2f}")

# Store results for this dataset
    results = dataset_results

# Prepare to save additional comparison plots
    comparison_folder = os.path.join(base_folder, 'Comparison_Plots')
    os.makedirs(comparison_folder, exist_ok=True)

# Plot Precision-Recall Scores Comparison for this dataset
pr_scores = [dataset_results[algorithm_name]["pr_curve"] for algorithm_name in algorithms]
pr_scores = [[pr_auc for _, _, pr_auc in pr_curve] for pr_curve in pr_scores]
algorithm_names = list(algorithms.keys())

plt.figure(figsize=(10, 6))
plt.boxplot(pr_scores, labels=algorithm_names)
plt.title('Precision-Recall Scores Comparison')
plt.ylim(0, 1)
plt.xlabel('Algorithms')
plt.ylabel('Precision-Recall Score')
plt.axhline(y=expected_prc, color='r', linestyle='--', label=f'Expected PRC = {expected_prc:.2f}')
plt.legend()
plt.grid(True)
boxplot_path = os.path.join(comparison_folder, '_Precision_Recall_Comparison.png')
plt.show()
plt.savefig(folder_to_save + '/Precision_Recall_Comparison.png')
plt.close()

                # Plot ROC curves for each algorithm for this dataset
plt.figure(figsize=(10, 8))
for algorithm_name in algorithm_names:
    mean_fpr = np.linspace(0, 1, 100)
    mean_tpr = np.mean([np.interp(mean_fpr, fpr, tpr) for fpr, tpr, _ in dataset_results[algorithm_name]["roc_curve"]], axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, label=f'{algorithm_name} (AUC = {mean_auc:.2f})')

plt.plot([0, 1], [0, 1], linestyle='--', color='black')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(f'ROC Curve')
plt.legend(loc='lower right')
plt.grid(True)
roc_comparison_path = os.path.join(comparison_folder, '_ROC_Curves_Comparison.png')
plt.savefig(folder_to_save + '/_ROC_Curves_Comparison.png')
plt.show()
plt.close()

                # Create box plots for AUC scores for this dataset
auc_scores = [dataset_results[algorithm_name]["roc_curve"] for algorithm_name in algorithms]
auc_scores = [[roc_auc for _, _, roc_auc in roc_curve] for roc_curve in auc_scores]

plt.figure(figsize=(10, 6))
plt.boxplot(auc_scores, labels=algorithm_names)
plt.title('AUC Scores Comparison')
plt.ylim(0, 1)
plt.xlabel('Algorithms')
plt.ylabel('AUC Score')
plt.grid(True)
auc_comparison_path = os.path.join(comparison_folder, '_AUC_Scores_Comparison.png')
plt.savefig(folder_to_save + '/_AUC_Scores_Comparison.png')
plt.show()
plt.close()

                # Show results
for algorithm_name in algorithms:
    print(f"Results for {algorithm_name}:")
    print("AUC scores:", [roc_auc for _, _, roc_auc in dataset_results[algorithm_name]["roc_curve"]])
    print("Precision-Recall scores:", [pr_auc for _, _, pr_auc in dataset_results[algorithm_name]["pr_curve"]])