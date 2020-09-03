import pandas as pd
import numpy as np
import glob

# Author: Abdullah Khan (A.Khan.4@bham.ac.uk)

# Define path to collection of exomes:
exome_path = './data/sample_exomes/*.csv'

# Define column headers for exome files.
# Note '1000g_freq', 'EVS_freq', 'Novelty', 'cosequence' are essential for this script
column_names = ['sample_ID', 'chr', 'start', 'end', 'ref_allele', 'alt_allele', 'zygosity',
                'exonic_or_splicing', 'gene', 'cosequence', 'details_code_change', 'SIFT',
                'Novelty', 'dbSNP134', '1000g_freq', 'EVS_freq', 'no_hom_obv',
                'no_het_obv', 'QC', 'Phred', 'VCF_det', 'VCF_det2',
                'VCF_det3']

# empty dataframe filled within the loop
sum_rare = pd.DataFrame()

# for loop which reads csv in path, appends column titles, replaces '.' values.
# For fname (filename) in the above exome path...
for fname in glob.glob(exome_path):
    # name each column - adjust this depending on your exome output.
    exome_sample = pd.read_csv(fname, names=column_names)

    # name each column - adjust this depending on your exome output.
    exome_sample.replace('.', np.nan, inplace=True)
    exome_sample[exome_sample.isna().any(axis=1)]
    exome_sample["1000g_freq"] = pd.to_numeric(exome_sample["1000g_freq"])
    exome_sample["EVS_freq"] = pd.to_numeric(exome_sample["EVS_freq"])
    # filters data to rare according to a user defined threshold in relevant columns, including novel variants.
    exome_is_rare = exome_sample[(exome_sample['1000g_freq'] <= 0.0001) | (exome_sample['EVS_freq'] <= 0.0001) | (
            exome_sample['Novelty'] == 'NOVEL')]
    exome_is_rare_ns = exome_is_rare[exome_is_rare['cosequence'] != 'synonymous SNV']
    exome_rare_summary = exome_is_rare_ns[column_names]
    print(exome_rare_summary.info())
    # exome_rare_summary.to_csv(fname.replace('.csv',' rare.csv'))
    sum_rare = pd.concat([sum_rare, exome_rare_summary], axis=0)

# Exports the above into a single csv file which with rare and novel variants
sum_rare.to_csv('./data/sum_rare.csv')
