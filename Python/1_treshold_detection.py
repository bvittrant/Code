#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
    Ceci est un test
"""
################################################################################

# Importing library

import os
import pandas as pd
import numpy as np

################################################################################

# Setting working directory

os.getcwd()
os.chdir('/home/benjamin/Desktop/test')

# Import data

df2 = pd.read_table('tmp3.tsv', sep='\t', header=None,  names=['Gene_list', 'HKG'])
df2 = df2.replace('', np.NaN)

# Loading matrice with results
dfTCGA_raw = pd.read_table('data/TCGA_results_raw.txt', sep='\t',index_col=0)
dfTCGA_raw.to_csv('dfTCGA_raw.tsv', sep='\t')
dfTCGA_raw_transpose = dfTCGA_raw.transpose()
dfTCGA_raw_transpose.to_csv('dfTCGA_raw_transpose.tsv', sep='\t')

dflong = pd.read_table('data/LONG_results_raw.txt', sep='\t', index_col=0)

#data_frame_raw_transpose[[]] # rows names
#data_frame_raw_transpose.columns.values # colnames

# Removing rows with only 0
dfTCGA_clean_0 = dfTCGA_raw.loc[(dfTCGA_raw!=0).any(axis=1)]
dfTCGA_clean_0.to_csv('dfTCGA_clean_0.tsv', sep='\t')

dfTCGA_clean_0_transpose = dfTCGA_clean_0.transpose()
dfTCGA_clean_0_transpose.to_csv('dfTCGA_clean_0_transpose.tsv', sep='\t')

# Removing rows with all values inferior or egal to 9

dfTCGA_rawTMP = dfTCGA_raw
dfTCGA_rawTMP[dfTCGA_rawTMP <= 9] = 0
dfTCGA_clean_lt9 = dfTCGA_rawTMP.loc[(dfTCGA_rawTMP!=0).any(axis=1)]
del dfTCGA_rawTMP

dflong_rawTMP = dflong
dflong_rawTMP[dflong_rawTMP <= 9] = 0
dflong_clean_lt9 = dflong_rawTMP.loc[(dflong_rawTMP!=0).any(axis=1)]
del dflong_rawTMP

dfTCGA_clean_lt9.to_csv('data/dfTCGA_clean_lt9.tsv', sep='\t')
dflong_clean_lt9.to_csv('data/dflong_clean_lt9.tsv', sep='\t')


dfTCGA_clean_lt9_transpose = dfTCGA_clean_lt9_transpose()
dfTCGA_clean_lt9_transpose.to_csv('dfTCGA_clean_lt9_transpose.tsv', sep='\t')

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Removing rows with all values inferior or egal to 30
dfTCGA_rawTMP = dfTCGA_raw
dfTCGA_rawTMP[dfTCGA_rawTMP <= 30] = 0
dfTCGA_clean_lt30 = dfTCGA_rawTMP.loc[(dfTCGA_rawTMP!=0).any(axis=1)]
del dfTCGA_rawTMP

dfTCGA_clean_lt30.to_csv('dfTCGA_clean_lt30.tsv', sep='\t')

dfTCGA_clean_lt30_transpose = dfTCGA_clean_lt30.transpose()
dfTCGA_clean_lt30_transpose.to_csv('dfTCGA_clean_lt30_transpose.tsv', sep='\t')

# Removing rows with all values inferior or egal to 1000
dfTCGA_rawTMP = dfTCGA_raw
dfTCGA_rawTMP[dfTCGA_rawTMP <= 1000] = 0
dfTCGA_clean_lt1000 = dfTCGA_rawTMP.loc[(dfTCGA_rawTMP!=0).any(axis=1)]
del dfTCGA_rawTMP

dfTCGA_clean_lt1000.to_csv('dfTCGA_clean_lt1000.tsv', sep='\t')

dfTCGA_clean_lt1000_transpose = dfTCGA_clean_lt1000.transpose()
dfTCGA_clean_lt1000_transpose.to_csv('dfTCGA_clean_lt1000_transpose.tsv', sep='\t')