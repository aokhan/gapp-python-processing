import pandas as pd
import numpy as np

# Author: Abdullah Khan (A.Khan.4@bham.ac.uk)
# THIS SCRIPT PRODUCES A LIST OF RARE GENES BASED ON KNOWN DISEASE CAUSING GENES.
# You will need to upload a .csv file with all of your potential rare or novel variants and a .csv with 'known' genes.


# Path to sum_rare.csv produced by 1_GAPP_Filter_Rarity.py
rare_genes_path = './data/sum_rare.csv'
# Path to file with with all of your potential rare or novel variants and a .csv with 'known' genes.
known_genes_path = './data/InheritedBleedingDisorder_ALL.csv'

# Read in csv files
rare_genes = pd.read_csv(rare_genes_path)
known_genes = pd.read_csv(known_genes_path)

# Check that you've uploaded your data frames correctly.
# print(rare_genes.head())
# print(known_genes.head()

rare_genes_list = rare_genes['gene']
known_genes_list = known_genes['gene'].to_list()

print(rare_genes_list.head())
print(known_genes_list)

# positive_known_genes = rare_genes[rare_genes['genes'] == known_genes['genes']]
positive_known_genes = rare_genes.isin(known_genes_list)
list_positives = rare_genes[positive_known_genes.gene.eq(True)]

print(positive_known_genes.head())
print(list_positives.head())

list_positives.to_csv('./data/list_positives.csv')

summary_list_positives = list_positives[['sample_ID','chr','start','end','ref_allele','alt_allele','exonic_or_splicing','gene','cosequence','details_code_change','SIFT','Novelty','1000g_freq','EVS_freq','no_hom_obv','no_het_obv']]
print(summary_list_positives)
summary_list_positives.to_csv('./data/summary_list_positives.csv')
