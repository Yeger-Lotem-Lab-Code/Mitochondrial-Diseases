import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold, RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, make_scorer

# Load dataset
dataset = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\MOSTCO~3_imputed_dataset.csv").iloc[:,1:]

# Process dataset
def process_dataset(data):
    # Separate features, labels, and groups
    '''
    not_important_features = ['hg19_Chromosome_chr15', 'hg19_Chromosome_chr13',
                              'mitochondrial small ribosomal subunit',
                              'hg19_Chromosome_chr12',
                              'hg19_Chromosome_chr4',
                              'hg19_Chromosome_chr20',
                              'Metabolism of substrates / Metabolism of cofactors',
                              'hg19_Chromosome_chr14',
                              'hg19_Chromosome_chr22',
                              'hg19_Chromosome_chr7',
                              'mitochondrial ribosome',
                              'mitochondrial large ribosomal subunit',
                              'hg19_Chromosome_chr18',
                              'Metabolism of toxic compounds',
                              'hg19_Chromosome_chr21',
                              'OXPHOS subunits, assembly factors, and electron carriers',
                              'hg19_Chromosome_chrM',
                              'Metabolism of cofactors']
     '''
    X = data.iloc[:, 4:]
    #X = X.drop(columns=not_important_features)
    #X.replace([np.inf, -np.inf], np.nan, inplace=True)
    #X = X.apply(lambda x: x.fillna(x.median()))
    y = data['label']
    groups = data['Ensembl_ID_x']
    return X, y, groups

X, y, groups = process_dataset(dataset)

# Define the model
model = RandomForestClassifier(random_state=42)

# Define the parameter grid
param_distributions = {
    'n_estimators': [50, 100, 200, 500],
    'max_features': ['auto', 'sqrt', 'log2'],
    'max_depth': [None, 10, 20, 30, 40, 50],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'bootstrap': [True, False]
}

# Define the cross-validation strategy
cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)

# Define the RandomizedSearchCV
random_search = RandomizedSearchCV(
    estimator=model,
    param_distributions=param_distributions,
    n_iter=100,  # Number of parameter settings that are sampled
    scoring=make_scorer(roc_auc_score),
    n_jobs=-1,  # Use all available cores
    cv=cv,
    verbose=2,
    random_state=42
)

# Fit the RandomizedSearchCV
random_search.fit(X, y, groups=groups)

# Print the best parameters and best score
print(f"Best Parameters: {random_search.best_params_}")
print(f"Best Score: {random_search.best_score_}")

# Get the best model
best_model = random_search.best_estimator_