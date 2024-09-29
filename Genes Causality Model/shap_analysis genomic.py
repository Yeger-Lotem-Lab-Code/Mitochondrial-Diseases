import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
import shap
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score

def process_dataset(data):
    # Separate features and labels
    X = data.iloc[:, 3:]
    #X = X.drop(columns=['Colon1_mean proact', 'Whole Colon_mean proact'])
    y = data.iloc[:,2]
    return X, y


def SHAP_analysis(X, y):
    model = XGBClassifier(use_label_encoder=False, eval_metric='logloss',colsample_bytree = np.float64(0.898027900926573), learning_rate= np.float64(0.02738831867639583), max_depth= 6, n_estimators= 175, subsample= np.float64(0.9315111447751979))  # Change the model for the data you're working on
    # Create SHAP explainer
    model.fit(X, y)
    explainer = shap.TreeExplainer(model)

    # Calculate SHAP values
    shap_values = explainer.shap_values(X)
    plt.figure(figsize=(10, 8))  # Increase figure size
    # Create the SHAP summary plot
    shap.summary_plot(shap_values, X, plot_size=(12, 8), plot_type="dot", show=False, max_display=10)
    # Customize the font sizes and other aesthetics
    plt.gca().tick_params(labelsize=12)  # Increase tick label size
    plt.xlabel('SHAP value (impact on model output)', fontsize=14)  # Increase x-axis label size
    plt.title('SHAP Summary Plot', fontsize=16)  # Add and size the title
    #plt.colorbar(label='Feature Value')  # Label for the color bar
    # Adjust layout for better spacing
    plt.tight_layout()
    #plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\SHAP\SHAP Summary Plot.png")  # Display the updated plot
    plt.show()


    # Bar plot
    plt.figure(figsize=(10, 8))  # Adjust the figure size as needed
    shap.summary_plot(shap_values, X, plot_size=(12, 8), plot_type='bar', show=False, max_display=10)
    # Customize the font sizes and other aesthetics
    plt.gca().tick_params(labelsize=12)  # Increase tick label size
    plt.xlabel('Mean |SHAP Value| (average impact on model output)', fontsize=14)  # Increase x-axis label size
    plt.title('Top 10 Feature Importances Based on SHAP Values', fontsize=16)  # Add and size the title
    # Adjust layout for better spacing
    plt.tight_layout()
    #plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Model 1 - Genes Causality\SHAP\SHAP Features Importance Plot.png")  # Display the updated plot
    # Display the updated plot
    plt.show()





def main():
    # Load dataset
    dataset_path = r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv"
    Dataset = pd.read_csv(dataset_path)

    # Process dataset
    X, y = process_dataset(Dataset)

    # Perform SHAP analysis
    SHAP_analysis(X, y)


if __name__ == "__main__":
    main()