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
d_MutByPatient = pd.read_table('results_files/SNP_count_mutation_by_patient_IGS_KIRDL.tsv')

# Add col BCR to DF
d_MutByPatient['BCR'] = d_clin['BCR_60']
d_MutByPatient['TNM'] = d_clin['pathological_grade']
d_MutByPatient['GS'] = d_clin['gleason']
d_MutByPatient['MTBCR'] = d_clin['Months_to_BCR']

###############################################################################

# Boxplot properties

boxprops = dict(linestyle='-', linewidth=4, color='k')
medianprops = dict(linestyle='-', linewidth=4, color='k')
flierprops = dict(marker='o', markerfacecolor='green', markersize=8,
                  linestyle='none')
whiskerprops = dict(linestyle='--', linewidth=3, color='k')
capprops = dict(linestyle='-', linewidth=2, color='red')

# Boxplot

plot_BCR = d_MutByPatient.boxplot(column='mutation_counts', by='BCR', 
    figsize=(10,8),boxprops=boxprops,medianprops=medianprops,showmeans=True,
    flierprops = flierprops,  whiskerprops = whiskerprops, capprops = capprops)
fig = plot_BCR.get_figure()
fig.savefig("pictures/SNP_BP_BCR_IGS_KIRDL.png")

d_MutByPatient.groupby('BCR').boxplot(column=['mutation_counts'],
    by='GS',figsize=(10,8),boxprops=boxprops,medianprops=medianprops,showmeans=True,
    flierprops = flierprops,  whiskerprops = whiskerprops, capprops = capprops)
plt.savefig("pictures/SNP_BP_GS_IGS_KIRDL.png")

d_MutByPatient.groupby('BCR').boxplot(column=['mutation_counts'],
    by='TNM',figsize=(10,8),boxprops=boxprops,medianprops=medianprops,showmeans=True,
    flierprops = flierprops,  whiskerprops = whiskerprops, capprops = capprops)
plt.savefig("pictures/SNP_BP_TNM_IGS_KIRDL.png")

# Scatter plot (better for continuous data)

color_dict = {1:'red',0:'blue'}

plot_MTBCR = d_MutByPatient.plot.scatter(x='mutation_counts', y='MTBCR', figsize=(10,8),
    color=[ color_dict[i] for i in d_MutByPatient['BCR'] ])
fig = plot_MTBCR.get_figure()
fig.savefig("pictures/SNP_BP_MTBCR_IGS_KIRDL.png")









