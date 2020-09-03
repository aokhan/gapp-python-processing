# GAPP Python Processing

Python scripts for analysing exome data, correlating to a list of clinically relevant bleeding genes, and cross indexing with RNA seq data (analysed using DESeq2). 

These scripts were used to analyse GAPP data from 129 exomes. The First script collates all rare and novel variants into a master sheet which is used for further analysis. 

The second script filters this list of variants and cross indexes them with a list of clinical bleeding genes. This highlights novel and rare variants in known bleeding genes and exports a .csv with this as a summary. 

The final script uploads RNA seq data and cross references it with rare variants in the patient cohort. This script includes a number of graphing options which will stratify the data and provide heat maps and dendograms, as well as scatter plots of genes which are expressed. 

If you need any help with this please contact a.khan.4@bham.ac.uk or j.a.pike@bham.ac.uk

## Dependencies

Scripts tested with Python 3.7, Numpy 1.19.1, Pandas 1.1.1, glob2 0.7, matplotlib 3.3.1, scikit-learn 0.23.2, seaborn 0.10.1, scipy 1.5.2

## Citation

If these scripts are useful please cite our publication: *under review*
