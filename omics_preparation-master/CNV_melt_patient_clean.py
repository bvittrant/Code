################ Preparing CNV set ##################

###############################################################################

print('######################################################################')
print('\nImporting library\n')
print('######################################################################')

# Import it as 'pd'
import pandas as pd

# Import numpy as np
import numpy as np

# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

# import OS
import os

# import tqdm
from tqdm import tqdm

print('######################################################################')
print('\nSetting WD\n')
print('######################################################################')

# Prints the working directory
print('\nYour working directory is :\n'+os.getcwd())

# Set working directory
os.chdir('/home/benjamin/Desktop/professionnal/full_TCGA_data/TCGA_other_data/')
cwd = os.getcwd()

# Prints the working directory
print('Now your working directory is :\n'+os.getcwd())

###############################################################################

d_cnv = pd.read_table(filepath_or_buffer='results_files/CNV_full_merged.tsv', sep='\t', 
	header=0)
print('\nYour CNV data frame has this structure :')
print(d_cnv.shape)

###############################################################################

# Remove line with 0

d_cnv_1 = d_cnv.loc[(d_cnv.ix[:,4:] != 0 ).any(axis=1),:]
d_cnv.to_csv(path_or_buf='results_files/CNV_full_merged_clean_0.tsv',
	sep ='\t', index = False, header=True)

# Remove line with values between -1 & 1
d_cnv_1 = d_cnv.loc[(d_cnv.ix[:,4:] > 1 ).any(axis=1),:]
print(d_cnv_1.shape)
d_cnv_2 = d_cnv.loc[(d_cnv.ix[:,4:] < -1 ).any(axis=1),:]
print(d_cnv_2.shape)

frames = [d_cnv_1,d_cnv_2]
d_cnv_3 = pd.concat(frames)
print(d_cnv_3.shape)

d_cnv_3.to_csv(path_or_buf='results_files/CNV_full_merged_clean_-1_1.tsv',
	sep ='\t', index = False, header=True)