import pandas as pd
import numpy as np
from sklearn.base import clone
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV
from sklearn.metrics import precision_recall_curve, roc_auc_score, make_scorer, auc
from scipy.stats import uniform, randint
import xgboost as xgb

# Load your dataset
dataset = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv")

X = dataset.drop(columns=['Label', 'ENSG ID', 'Gene_x'])
y = dataset['Label']

# Split data into training and validation sets
X_train_full, X_valid, y_train_full, y_valid = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

# Reset index
X_train_full.reset_index(drop=True, inplace=True)
y_train_full.reset_index(drop=True, inplace=True)

# Define the XGBoost model
model = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', early_stopping_rounds=10)

# Define the parameter grid for hyperparameter tuning
param_dist = {
    'n_estimators': randint(50, 200),
    'learning_rate': uniform(0.01, 0.2),
    'max_depth': randint(3, 7),
    'subsample': uniform(0.6, 0.4),
    'colsample_bytree': uniform(0.6, 0.4)
}

# Define the scoring metrics
def pr_auc_score(y_true, y_pred_proba):
    precision, recall, _ = precision_recall_curve(y_true, y_pred_proba)
    return auc(recall, precision)

scoring = {
    'AUC': 'roc_auc',
    'PR_AUC': make_scorer(pr_auc_score, needs_proba=True)
}

# Create the StratifiedKFold object
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)

# Perform the randomized search
random_search = RandomizedSearchCV(estimator=model, param_distributions=param_dist, scoring=scoring, refit='AUC', cv=cv, n_iter=20, verbose=1, n_jobs=-1)
random_search.fit(X_train_full, y_train_full, eval_set=[(X_valid, y_valid)], verbose=False)

# Get the best parameters and scores
best_params = random_search.best_params_
best_auc_score = random_search.best_score_
best_pr_auc_score = random_search.cv_results_['mean_test_PR_AUC'][random_search.best_index_]

print("Best parameters found: ", best_params)
print("Best AUC score: ", best_auc_score)
print("Best PR AUC score: ", best_pr_auc_score)

# Get the best model and update it with early stopping parameters
best_model = xgb.XGBClassifier(
    **best_params,
    use_label_encoder=False,
    eval_metric='logloss',
    early_stopping_rounds=10
)

# Define a function to fit the model and evaluate
def fit_and_score(estimator, X_train, X_test, y_train, y_test):
    """Fit the estimator on the train set and score it on both sets"""
    estimator.fit(X_train, y_train, eval_set=[(X_test, y_test)], verbose=False)
    train_score = roc_auc_score(y_train, estimator.predict_proba(X_train)[:, 1])
    test_score = roc_auc_score(y_test, estimator.predict_proba(X_test)[:, 1])
    return estimator, train_score, test_score

# Evaluate the best model on the entire dataset using custom cross-validation with early stopping
auc_scores = []
pr_auc_scores = []

for train_idx, test_idx in cv.split(X_train_full, y_train_full):
    X_train, X_test = X_train_full.iloc[train_idx], X_train_full.iloc[test_idx]
    y_train, y_test = y_train_full.iloc[train_idx], y_train_full.iloc[test_idx]

    est, train_score, test_score = fit_and_score(clone(best_model), X_train, X_test, y_train, y_test)

    auc_scores.append(test_score)
    pr_auc_scores.append(pr_auc_score(y_test, est.predict_proba(X_test)[:, 1]))

print(f"Mean AUC score from cross-validation: {np.mean(auc_scores)}")
print(f"Mean PR AUC score from cross-validation: {np.mean(pr_auc_scores)}")