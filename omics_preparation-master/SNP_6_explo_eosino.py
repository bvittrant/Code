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

print('############################################################################')
print('\nSetting WD\n')
print('############################################################################')

print('\nYour working directory is :\n'+os.getcwd()) # Prints the working directory

# Set working directory
os.chdir('/home/benjamin/Desktop/professionnal/full_TCGA_data/TCGA_other_data/')
cwd = os.getcwd()

print('Now your working directory is :\n'+os.getcwd()) # Prints the working directory

#####################################################################################

# DL data and preparing data frame for further investigations

d_MutByPatient = pd.read_table('results_files/SNP_count_mutation_by_patient.tsv')
d_clin = pd.read_table('clinical/TCGA_clinique_45_reduced.tsv')

SNP_45_gene_type_clean

d_SNP_45_gene_type_clean = pd.read_table('results_files/SNP_45_gene_type_clean.tsv')
d_SNP_mutation_gene = pd.read_table('results_files/SNP_mutation_per_gene.tsv', sep='\t')
d_SNP_mutation_gene_type = pd.read_table('results_files/SNP_mutation_per_gene_type.tsv', sep='\t')

d_eosino = pd.read_table('data/eosino_gene.txt', sep='\t')

#####################################################################################

# Eosino exploration

d_inter = pd.merge(d_SNP_mutation_gene_type['GENE'].to_frame(),
    d_eosino['Gene name'].to_frame(), left_on=['GENE'], right_on=['Gene name'])
    
# NULL
    
#####################################################################################










