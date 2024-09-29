import csv
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Set seaborn color palette
sns.set_palette("colorblind")
palette = sns.color_palette("colorblind")


# Function to perform KS test and generate cumulative frequency graph
def glob_mito_PC_ks(arr1, arr2):
    '''Calculates KS test for 2 groups of genes and creates a cumulative graph'''

    # Sort arrays once
    arr1 = np.sort(arr1)
    arr2 = np.sort(arr2)

    # Combine and sort the joined array for the graph limits
    joined = np.sort(np.concatenate([arr1, arr2]))

    # Set colors from seaborn colorblind palette
    mito_color = palette[0]  # First color for mitochondrial genes
    pro_color = palette[1]   # Second color for protein coding genes

    # Controls the bins of the graph (increase for smoother graph)
    bins = np.append(np.linspace(joined.min(), joined.max(), 1000), [np.inf])

    plt.hist(arr1, bins=bins, density=1, histtype='step', cumulative=True, color=mito_color, linewidth=3)
    plt.hist(arr2, bins=bins, density=1, histtype='step', cumulative=True, color=pro_color, linewidth=3)

    # KS Test statistics
    print(scipy.stats.ks_2samp(arr1, arr2), 'mito vs prot')
    print('Number of mito = ', len(arr1))
    print('Number of protein-coding genes = ', len(arr2))

    # Set limits for the X-axis if needed
    plt.xlim(joined.min(), joined.max())
    plt.ylim(0, 1)

    # Save figure
    plt.legend(('mitochondrial genes', 'protein coding genes'), loc='upper left')
    plt.ylabel('Cumulative frequency')
    plt.xlabel('Impact on growth (CRISPR score)')
    plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Essentiality_mito_vs_pro_cod.png")
    plt.show()


# Load the generated file to compare mitochondrial vs protein-coding genes
depmap = pd.read_csv(r'DepMap_CRISPR_min_score_final_with_mito_12-2-23.csv')

# Group the data using vectorized operations
arr1 = depmap.loc[depmap['mito'] == "Yes", 'Min'].dropna().astype(float).values
arr2 = depmap.loc[depmap['mito'] == "No", 'Min'].dropna().astype(float).values

# Perform KS test and plot the cumulative graph
glob_mito_PC_ks(arr1, arr2)