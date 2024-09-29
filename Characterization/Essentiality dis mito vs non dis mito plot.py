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
def glob_dismito_PC_ks(arr1, arr2):
    '''Calculates KS test for 2 groups of genes and creates a cumulative graph'''

    # Sort arrays once
    arr1 = np.sort(arr1)
    arr2 = np.sort(arr2)

    # Combine and sort the joined array for the graph limits
    joined = np.sort(np.concatenate([arr1, arr2]))

    # Set colors from seaborn colorblind palette
    mito_color = palette[0]  # First color for disease mitochondrial genes
    pro_color = palette[1]   # Second color for non disease mitochondrial genes

    # Controls the bins of the graph (increase for smoother graph)
    bins = np.append(np.linspace(joined.min(), joined.max(), 1000), [np.inf])

    plt.hist(arr1, bins=bins, density=1, histtype='step', cumulative=True, color=mito_color, linewidth=3)
    plt.hist(arr2, bins=bins, density=1, histtype='step', cumulative=True, color=pro_color, linewidth=3)

    # KS Test statistics
    print(scipy.stats.ks_2samp(arr1, arr2), 'dis mito vs non dis mito')
    print('Number of dis mito = ', len(arr1))
    print('Number of non dis mito = ', len(arr2))

    # Set limits for the X-axis if needed
    plt.xlim(joined.min(), joined.max())
    plt.ylim(0, 1)

    # Save figure
    plt.legend(('disease mitochondrial genes', 'non disease mitochondrial genes'), loc='upper left')
    plt.ylabel('Cumulative frequency')
    plt.xlabel('Impact on growth (CRISPR score)')
    plt.savefig(r"C:\Users\liorr\Dropbox\Lior\Final Results\Characterization\Essentiality_dis mito_vs_non dis mito.png")
    plt.show()

# Load the generated file to compare mitochondrial vs protein-coding genes
depmap = pd.read_csv(r"C:\Users\liorr\Documents\pythonProject2\DepMap_CRISPR_min_score_final_with_Dismito_26-08-24.csv")

# Group the data using vectorized operations
arr1 = depmap.loc[depmap['Dis'] == "Yes", 'Min'].dropna().astype(float).values
arr2 = depmap.loc[depmap['Dis'] == "No", 'Min'].dropna().astype(float).values

# Perform KS test and plot the cumulative graph
glob_dismito_PC_ks(arr1, arr2)