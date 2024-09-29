import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
#from xgboost import XGBClassifier
#from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
import shap
import matplotlib.pyplot as plt
#from sklearn.model_selection import train_test_split
#from sklearn.metrics import roc_auc_score, average_precision_score


def process_dataset(data):
    # Separate features and labels
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
                              'Metabolism of cofactors'
                              ]
    '''
    X = data.iloc[:, 4:]
    #X = X.drop(columns=not_important_features)
    y = data['label']
    return X, y


def SHAP_analysis(X, y, max_depth):
    model = RandomForestClassifier(n_estimators= 50, min_samples_split= 2, min_samples_leaf= 2, max_features= 'sqrt', max_depth= max_depth, bootstrap= False)
    #model = RandomForestClassifier(n_estimators= 50, min_samples_split= 2, min_samples_leaf= 2, max_features= 'sqrt', max_depth= 25, bootstrap= False)  # Change the model for the data you're working on
    # Create SHAP explainer
    model.fit(X, y)
    explainer = shap.Explainer(model)

    # Calculate SHAP values
    shap_values = explainer.shap_values(X)
    # Summary plot
    plt.figure(figsize=(20, 20))
    shap.summary_plot(shap_values[:, :, 1], X, max_display=10)
    #save the plot
    #plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\SHAP Summary Plot.png")
    plt.close()
    # Bar plot
    plt.figure(figsize=(20, 20))
    #plt.text(1.05, 1.01, f'Max Depth: {max_depth}', horizontalalignment='right', verticalalignment='top',
             #transform=plt.gca().transAxes, fontsize=12)
    shap.summary_plot(shap_values[:, :, 1], X, plot_type='bar', max_display=10)
    #save the plot
    #plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\SHAP Features Importance Plot.png")
    plt.close()
# Get the top 10 most contributing features
    top_features = [feature for feature in X.columns if feature in ['mean proact score', 'Num_Elevated_Interactors',
                        'Age Coefficient', 'Tissue_Childhood', 'Preferential_Expression', 'diffnetmed',
                        'TS score', 'neutralmediancrisprscore', 'egene', 'Time_Neonatal']]
    # Loop through each pair of top features and create dependence plots
    for i in range(len(top_features)):
        for j in range(i + 1, len(top_features)):
            feature_x = top_features[i]
            feature_y = top_features[j]

            # Create SHAP dependence plot
            shap.dependence_plot(feature_x, shap_values[:, :, 1], X, interaction_index=feature_y, show=False)
            plt.savefig(
                rf"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\SHAP dependence\Dependence_{feature_x}_{feature_y}.png")



def main():
    # Load dataset
    dataset_path = r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 2 - Causality in Tissues\MOSTCO~3_imputed_dataset.csv"
    Dataset = pd.read_csv(dataset_path).iloc[:, 1:]

    # Process dataset
    X, y = process_dataset(Dataset)

    max_depth = 40
    SHAP_analysis(X, y, max_depth)
if __name__ == "__main__":
    main()
