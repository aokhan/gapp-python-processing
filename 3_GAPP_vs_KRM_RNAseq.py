import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as shc

# Author: Abdullah Khan (A.Khan.4@bham.ac.uk)
# THIS SCRIPT will cross reference patient data with an RNA seq model.

# Read sum_rare.csv produced by 1_GAPP_Filter_Rarity.py
gapp_rare = pd.read_csv('./data/sum_rare.csv')
# Read csv file with known genes for filtering against- h_rount_MK_pro, h_proMK_pro
known_genes = pd.read_csv('./data/InheritedBleedingDisorder_ALL.csv')

# csv files pulled from KRM data with paired fold log change

proMK_pro = pd.read_csv('./data/KRM_pro_mega_vs_pro.csv')
roundMK_pro = pd.read_csv('./data/KRM_round_mega_vs_pro.csv')

# filter to remove where there is no human homologue

h_proMK_pro = proMK_pro.dropna(axis=0, subset=['human_homolog'])
h_roundMK_pro = roundMK_pro.dropna(axis=0, subset=['human_homolog'])

h_proMK_pro.sort_values('human_homolog', ascending=True)
h_roundMK_pro.sort_values('human_homolog', ascending=True)

h_roundMK_pro.set_index(['human_homolog'], inplace=True)
h_proMK_pro.set_index(['human_homolog'], inplace=True)

sum_round_comp_3 = h_proMK_pro.merge(h_roundMK_pro, left_index=True, right_index=True)
# print(sum_round_comp_3.info())
#  RIGHT DATAFRAME _x IS ROUND MK V PRO, LEFT IS PRO_MK_ROUND
sum_round_comp_3.to_csv('./data/CHECK THIS SUM_3.csv')

frame3 = sum_round_comp_3.plot.scatter(x='log2FoldChange_y',
                                       y='log2FoldChange_x',
                                       c='baseMean_x',
                                       s=(sum_round_comp_3['baseMean_y'] / 200),
                                       colormap='gist_heat', alpha=0.4)
plt.ylim(-4, +3)
plt.xlim(-8, +4)
# frame2.axes.get_xaxis().set_visible(False)
plt.title('Fold changes in expression of known genes on proplatelet formation and release')
plt.xlabel('log2Fold Proplatelet MKs v. Released')
plt.ylabel('log2FoldChange Round MKs v. Released')

known_genes_list = known_genes['gene'].to_list()  # made a list of known genes
# now need to compare index to known_genes
RNA_seq_gene_list1 = h_roundMK_pro[h_roundMK_pro.index.isin(known_genes_list)]  # remember human homolog is the index
# THIS GIVES A DATAFRAME OF KNOWN GENES EXPRESSED FOR h_roundMK_pro...
RNA_seq_gene_list2 = h_proMK_pro[h_proMK_pro.index.isin(known_genes_list)]  # remember human homolog is the index
# print(RNA_seq_gene_list2)
RNA_seq_known = RNA_seq_gene_list2.merge(RNA_seq_gene_list1, left_index=True, right_index=True)
frame4 = RNA_seq_known.plot.scatter(x='log2FoldChange_y',
                                    y='log2FoldChange_x',
                                    c='baseMean_x',
                                    s=(RNA_seq_known['baseMean_y'] / 200),
                                    colormap='gist_heat', alpha=0.4)
plt.ylim(-4, +3)
plt.xlim(-8, +4)
# frame2.axes.get_xaxis().set_visible(False)
plt.title('Fold changes in expression of known genes on proplatelet formation and release')
plt.xlabel('log2Fold Proplatelet MKs v. Released')
plt.ylabel('log2FoldChange Round MKs v. Released')
plt.annotate("GATA1", (1.318354084, 1.183399092))
plt.annotate("RUNX1", (0.556647709, 0.233321203))
plt.annotate("TUBB1", (1.251054466, 0.895252317))
plt.annotate("THPO", (0.923839261, 0.895252317))
plt.annotate("MYH9", (1.372742608, 0.874800241))

RNA_seq_known.to_csv('./data/List_of_known_genes_v_RNA_seq.csv')

# filtering against GAPP genes - h_rount_MK_pro, h_proMK_pro
gapp_rare_list = gapp_rare['gene'].to_list()
gapp_RNA_gene_list1 = h_roundMK_pro[h_roundMK_pro.index.isin(gapp_rare_list)]
gapp_RNA_gene_list2 = h_proMK_pro[h_proMK_pro.index.isin(gapp_rare_list)]
gapp_RNA_merge = gapp_RNA_gene_list2.merge(gapp_RNA_gene_list1, left_index=True, right_index=True)

# print(gapp_RNA_merge.info())
# print(RNA_seq_known.info())
# print(h_roundMK_pro.info())


frame5 = gapp_RNA_merge.plot.scatter(x='log2FoldChange_y',
                                     y='log2FoldChange_x',
                                     c='baseMean_x',
                                     s=(gapp_RNA_merge['baseMean_y'] / 200),
                                     colormap='gist_heat', alpha=0.4)
plt.ylim(-4, +3)
plt.xlim(-8, +4)
# frame2.axes.get_xaxis().set_visible(False)
plt.title('Fold changes in expression of GAPP genes on proplatelet formation and release')
plt.xlabel('log2Fold Proplatelet MKs v. Released')
plt.ylabel('log2FoldChange Round MKs v. Released')

# Begin clustering these badboys using https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html and https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.cluster.hierarchy.dendrogram.html
gapp_RNA_dendo = gapp_RNA_merge[['log2FoldChange_y', 'log2FoldChange_x']]
print(gapp_RNA_dendo)
plt.figure(figsize=(10, 7))
plt.title('Dendrograms')
dend = shc.dendrogram(shc.linkage(gapp_RNA_dendo, method='ward'))
# all this stuff is from this blog https://www.analyticsvidhya.com/blog/2019/05/beginners-guide-hierarchical-clustering/
cluster = AgglomerativeClustering(n_clusters=4, affinity='euclidean', linkage='ward')
cluster.fit_predict(gapp_RNA_dendo)

gapp_RNA_merge['cluster'] = np.array(cluster.labels_)  # cluster labels added here
print(gapp_RNA_merge.head())

# these are commented out but they are used to make graphs if you need to.
plt.figure()
plt.scatter(gapp_RNA_dendo['log2FoldChange_y'], gapp_RNA_dendo['log2FoldChange_x'], s=0.5, c=cluster.labels_, alpha=0.4)
# seaborn clustermap
sns.clustermap(gapp_RNA_dendo)

plt.show()
plt.clf()
