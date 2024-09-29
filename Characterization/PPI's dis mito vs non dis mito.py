import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, ks_2samp
from statsmodels.stats.multitest import multipletests



#import data
genomic_ds = pd.read_csv(r"C:\Users\liorr\Dropbox\Lior\Mitochondrial Diseases - Genes Causality\OMIM Labeled\MitoCarta_OMIM\MitoCarta_Label_OMIM_subset_with_prefGorilla_imputed_dataset.csv")
#create df of mito genes and labeles from genomic_ds
mitocarta_genes_and_labels = genomic_ds[['ENSG ID', 'Label']]

trace_file = pd.read_csv(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\2_TRACE_Dataset_united_development_features.csv")
#filter trace to only mito genes
trace_file = trace_file[trace_file['Gene_ID'].isin(genomic_ds['ENSG ID'])]


#create df of mito genes and labeles
dismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 1]['ENSG ID']
non_dismito_genes = mitocarta_genes_and_labels[mitocarta_genes_and_labels['Label'] == 0]['ENSG ID']

# Filter for each PPI feature type across all tissues in the trace file
PPI_feature_types = ['num_interactors', 'num_elevated_interactors', 'num_specific_interactions']

# Initialize lists to store results for plotting


# Iterate over each PPI feature type
for PPI in PPI_feature_types:
    # Filter columns related to the current PPI feature type
    PPI_features = [col for col in trace_file.columns if PPI in col]
    p_values = []
    significance_markers = []
    # Initialize lists to store median values for each tissue
    dismito_medians = []
    nondismito_medians = []
    tissues = []

    # Iterate over each PPI feature column for the current feature type
    for feature in PPI_features:
        # Extract the tissue name from the column name (e.g., 'Brain_num_interactors' -> 'Brain')
        tissue = feature.split('_')[0]
        tissues.append(tissue)

        # Calculate median expression for each group (dismito and nondismito) for the current PPI feature
        dismito_values = trace_file[trace_file['Gene_ID'].isin(dismito_genes)][feature].dropna()
        nondismito_values = trace_file[trace_file['Gene_ID'].isin(non_dismito_genes)][feature].dropna()

        # Store median values for plotting
        dismito_medians.append(dismito_values.median())
        nondismito_medians.append(nondismito_values.median())

        # Perform the Mann-Whitney U test
        stat, p_value = mannwhitneyu(dismito_values, nondismito_values, alternative='two-sided')
        p_values.append(p_value)

        # Print the results of the test
        print(f"Tissue: {tissue}, PPI: {PPI}, Test: Mann-Whitney U, Group 1 (Dismito): n={len(dismito_values)}, "
              f"Group 2 (Non-Dismito): n={len(nondismito_values)}, p-value: {p_value}")

    # Adjust p-values for multiple testing correction using FDR
    _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

    # Determine significance markers after correction
    for p_value in corrected_p_values:
        if p_value < 0.05:
            significance_markers.append('*')
        else:
            significance_markers.append('')

    # Prepare data for plotting median expression of PPI features across tissues
    data = pd.DataFrame({
        'Tissue': tissues,
        'Dismito': dismito_medians,
        'Non-Dismito': nondismito_medians,
    })
    data_melted = data.melt(id_vars='Tissue', var_name='Gene Type', value_name='Median Value')

    # Plot median expression for current PPI type across tissues
    sns.set_palette("colorblind")
    plt.figure(figsize=(20, 10))
    ax = sns.barplot(x='Tissue', y='Median Value', hue='Gene Type', data=data_melted)

    # Set legend and title
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.title(f'Median Value of {PPI} Across Tissues for Dismito vs Non-Dismito Genes', pad=20)

    # Add asterisks above the bars
    #for i, marker in enumerate(significance_markers):
        #if marker:
            #ax.text(i, max(data_melted['Median Value']) + 1, marker, ha='center', va='bottom', color='black', fontsize=12)

    # Rotate x-axis labels and show plot
    plt.xticks(rotation=90)
    plt.xlabel('Tissue')
    plt.ylabel(f'Median {PPI} Level')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    #plt.savefig(
    #    fr"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\{PPI} across tissues dis mito vs non dis mito.png")

    plt.show()
    plt.close()