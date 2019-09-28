# SNP data frame treatmetn

print('################################################################################################')
print('\nImporting library\n')
print('################################################################################################')

# Import it as 'pd'
import pandas as pd

# Import numpy as np
import numpy as np

# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

# import OS
import os

#import tqdm
from tqdm import tqdm

# import scypi
from scipy.cluster.hierarchy import dendrogram, linkage

# import plotly
import plotly.plotly as py
import plotly.graph_objs as go

from pylab import figure, axes, pie, title, show

#####################################################################################

print('################################################################################################')
print('\nSetting WD\n')
print('################################################################################################')

print('\nYour working directory is :\n'+os.getcwd()) # Prints the working directory

# Set working directory
os.chdir('/home/benjamin/Desktop/professionnal/full_TCGA_data/TCGA_other_data/')
cwd = os.getcwd()

print('Now your working directory is :\n'+os.getcwd()) # Prints the working directory

#####################################################################################

# DL data and preparing data frame for further investigations

d_clin = pd.read_table('clinical/TCGA_clinique_45_reduced.tsv')

d_SNP = pd.read_table('results_files/SNP_45_gene_clean_BF.tsv')
d_SNP.head(n=5)
d_SNP.to_csv('results_files/SNP_45_gene_clean_BF.tsv', sep='\t', index=True)

d_mRNA = pd.read_table('results_files/TCGA_RNAseq_45_IfullG.tsv')
d_mRNA.head(n=5)

d_miRNA = pd.read_table('results_files/TCGA_miRNA_quant_cpm_45_threshold10.tsv')

d_CNV = pd.read_table('results_files/CNV_45_merged_clean_-1_1_row497zero.tsv') 
d_CNV.head(n=5)

d_SNP_mRNA = pd.merge(d_SNP,d_mRNA,on='barcode')
d_SNP_mRNA_miRNA = pd.merge(d_SNP_mRNA,d_miRNA,on='barcode')
d_SNP_mRNA_miRNA_CNV = pd.merge(d_SNP_mRNA_miRNA,d_CNV,on='barcode')
d_SNP_mRNA_miRNA_CNV_BCR = pd.merge(d_clin[['barcode','BCR_60']],d_SNP_mRNA_miRNA_CNV,on='barcode')

d_SNP_mRNA.to_csv('results_files/BF_SNP_mRNA.tsv', sep='\t', index=False)
d_SNP_mRNA_miRNA.to_csv('results_files/BF_SNP_mRNA_miRNA.tsv', sep='\t', index=False)
d_SNP_mRNA_miRNA_CNV.to_csv('results_files/BF_SNP_mRNA_miRNA_CNV.tsv', sep='\t', index=False)
d_SNP_mRNA_miRNA_CNV_BCR.to_csv('results_files/BF_SNP_mRNA_miRNA_CNV_BCR.tsv', sep='\t', index=False)

#####################################################################################

