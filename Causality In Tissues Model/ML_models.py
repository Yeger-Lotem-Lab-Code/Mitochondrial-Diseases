import os
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve, roc_curve, average_precision_score
import matplotlib.pyplot as plt
import numpy as np

# Define a function to calculate ROC and PR curves for an algorithm
def calculate_roc_pr_auc(algorithm, X, y, genes_list, num_folds):
    skf = StratifiedGroupKFold(n_splits=num_folds, shuffle=True, random_state=42)
    roc_curve_data = []
    pr_curve_data = []

    for train_idx, val_idx in skf.split(X, y, groups=genes_list):
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



# Define algorithms and specify which can handle NaNs directly
algorithms = {
    "Logistic Regression": {"model": LogisticRegression()},
    "XGBoost": {"model": XGBClassifier()},
    "Gradient Boosting": {"model": GradientBoostingClassifier()},
    "Random Forest": {"model": RandomForestClassifier()},
     "Random Forest tuned": {"model": RandomForestClassifier( n_estimators= 50, min_samples_split= 2, min_samples_leaf= 2, max_features= 'sqrt', max_depth= 40, bootstrap= False)}
}

# Read the dataset
dataset_path = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\MOSTCO~3_imputed_dataset.csv"
dataset = pd.read_csv(dataset_path).iloc[:, 1:]
genes_list = dataset['Ensembl_ID_x']
# Extract features and target
#mitocartafeatures = ['Peak_Intensity', 'ProteinLength', 'PGC_Induction_Score', 'CoexpressionGnfN50_Score', 'RickettsiaHomolog_Score_Homolog','RickettsiaHomolog_Score_NoHomolog', 'RickettsiaHomolog_Score_Ortholog', 'TargetP_Score',
#                     'YeastMitoHomolog_Score_HomologMitoHighConf', 'YeastMitoHomolog_Score_HomologMitoLowConf', 'YeastMitoHomolog_Score_NoMitoHomolog', 'YeastMitoHomolog_Score_OrthologMitoHighConf', 'YeastMitoHomolog_Score_OrthologMitoLowConf',
#                     'mitoPathway_OXPHOS', 'mitoPathway_Mitochondrial central dogma', 'mitoPathway_Metabolism', 'mitoPathway_Protein import, sorting and homeostasis','mitoPathway_Mitochondrial dynamics and surveillance', 'mitoPathway_Small molecule transport','mitoPathway_Signaling' ]  # Assuming features start from the 5th column
#trace_features = ['diffnetmed', 'egene', 'paralogsratiohighestidentity',
#       'mediantipapathways', 'neutralmediancrisprscore',
#       'neutralmedianrnaiscore', 'Preferential_Expression',
#       'Num_Elevated_Interactors', 'Num_Interactors',
#       'Num_Specific_Interactions']
#Go_features = dataset.columns[80:119]
#development_features = dataset.columns[17:31]
#other_features =dataset.columns[67:71]
#not_important_features = ['hg19_Chromosome_chr15', 'hg19_Chromosome_chr13',
#                              'mitochondrial small ribosomal subunit',
#                             'hg19_Chromosome_chr12',
#                              'hg19_Chromosome_chr4',
#                              'hg19_Chromosome_chr20',
#                              'Metabolism of substrates / Metabolism of cofactors',
#                              'hg19_Chromosome_chr14',
#                              'hg19_Chromosome_chr22',
#                              'hg19_Chromosome_chr7',
#                              'mitochondrial ribosome',
#                              'mitochondrial large ribosomal subunit',
#                              'hg19_Chromosome_chr18',
#                              'Metabolism of toxic compounds',
#                              'hg19_Chromosome_chr21',
#                              'OXPHOS subunits, assembly factors, and electron carriers',
#                              'hg19_Chromosome_chrM',
#                              'Metabolism of cofactors'
#                              ]
#trace_develop_mitocarta_features = list(development_features) + trace_features +mitocartafeatures

#features = dataset.loc[:, 4:]  # Assuming features start from the 5th column
features = dataset.iloc[:, 4:]
#features = features.drop(columns=not_important_features)
target = dataset['label']  # Replace 'score' with your actual target column name

# Store results
results = {algorithm_name: {"roc_curve": [], "pr_curve": []} for algorithm_name in algorithms}

# Define the folder structure for saving plots
base_folder = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues"
comparison_folder = os.path.join(base_folder, '5 fold cross-validation comparison')
os.makedirs(comparison_folder, exist_ok=True)

# Iterate over algorithms
for algorithm_name, algorithm_info in algorithms.items():
    model = algorithm_info["model"]

    # Call the function to calculate ROC and PR curves for each algorithm
    roc_curve_data, pr_curve_data = calculate_roc_pr_auc(model, features, target, genes_list, num_folds=5)

    # Store curve data for each algorithm
    results[algorithm_name]["roc_curve"] = roc_curve_data
    results[algorithm_name]["pr_curve"] = pr_curve_data

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
    ax2.axhline(y=expected_prc, color='black', linestyle='--', label=f'Expected PRC = {expected_prc:.2f}')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_title('Precision-Recall Curve')
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.legend(loc='upper right')
    ax2.grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # Save the figure
    save_path = os.path.join(comparison_folder, f'{algorithm_name}_ROC_PR.png')
    plt.savefig(save_path)
    plt.show()
    plt.close()
    y_pred_proba_full = model.predict_proba(features)[:, 1]

    # Create a DataFrame with actual labels and predicted probabilities
    results_df = pd.DataFrame({'Actual': target, 'Predicted_Probability': y_pred_proba_full})

    # Filter cases where the actual label is 0 and sort by predicted probability in descending order
    high_pred_not_labeled = results_df[results_df['Actual'] == 0].sort_values(by='Predicted_Probability',
                                                                              ascending=False)

    # Print the top 10 variables (features) with highest prediction values but not labeled as 1
    top_10_high_pred_not_labeled = high_pred_not_labeled.head(10)
    print(f"Top 10 high prediction values by {algorithm_name} but not labeled as 1:")
    print(top_10_high_pred_not_labeled)

# Plot Precision-Recall Scores Comparison
pr_scores = [results[algorithm_name]["pr_curve"] for algorithm_name in algorithms]
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
boxplot_path = os.path.join(comparison_folder, 'Precision_Recall_Comparison.png')
plt.savefig(boxplot_path)
plt.show()
plt.close()

# Plot ROC curves for each algorithm
plt.figure(figsize=(10, 8))
for algorithm_name in algorithm_names:
    mean_fpr = np.linspace(0, 1, 100)
    mean_tpr = np.mean([np.interp(mean_fpr, fpr, tpr) for fpr, tpr, _ in results[algorithm_name]["roc_curve"]], axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, label=f'{algorithm_name} (AUC = {mean_auc:.2f})')

plt.plot([0, 1], [0, 1], linestyle='--', color='black')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve Comparison')
plt.legend(loc='lower right')
plt.grid(True)
roc_comparison_path = os.path.join(comparison_folder, 'ROC_Curves_Comparison.png')
plt.savefig(roc_comparison_path)
plt.show()
plt.close()

# Create box plots for AUC scores
auc_scores = [results[algorithm_name]["roc_curve"] for algorithm_name in algorithms]
auc_scores = [[roc_auc for _, _, roc_auc in roc_curve] for roc_curve in auc_scores]

plt.figure(figsize=(10, 6))
plt.boxplot(auc_scores, labels=algorithm_names)
plt.title('AUC Scores Comparison')
plt.ylim(0, 1)
plt.xlabel('Algorithms')
plt.ylabel('AUC Score')
plt.grid(True)
auc_comparison_path = os.path.join(comparison_folder, 'AUC_Scores_Comparison.png')
plt.savefig(auc_comparison_path)
plt.show()
plt.close()

# Show results
for algorithm_name in algorithms:
    print(f"Results for {algorithm_name}:")
    print("AUC scores:", [roc_auc for _, _, roc_auc in results[algorithm_name]["roc_curve"]])
    print("Precision-Recall scores:", [pr_auc for _, _, pr_auc in results[algorithm_name]["pr_curve"]])

# Calculate and print the expected value of PRC
total_positive_instances = sum(target)
total_instances = len(target)
expected_prc = total_positive_instances / total_instances
print(f"Expected value of PRC: {expected_prc:.2f}")