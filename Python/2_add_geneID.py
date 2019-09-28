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
#import _thread

################################################################################

# Setting working directory

os.getcwd()
os.chdir('/home/benjamin/Desktop/test')

# Import data

dfTCGA = pd.read_table('data/dfTCGA_clean_lt9.tsv', sep='\t', index_col=0)
dflong = pd.read_table('data/dflong_clean_lt9.tsv', sep='\t', index_col=0)

# dftest = pd.read_table('dftest.tsv', sep='\t', index_col=0)

################################################################################

# Check colnames and rownames

dflong[[]] # rows names
dflong .columns.values # colnames
dfTCGA[[]] # rows names
dfTCGA.columns.values # colnames

################################################################################

# Go to Biomart and convert ID 
# http://useast.ensembl.org/biomart/martview/632ecb1fa82127a0c4b03db63e741b4d
# /home/benjamin/Downloads/ENSTtoXXX.txt

dfENSTtoXXX = pd.read_table('data/ENSTtoXXX.txt', sep='\t')
#dfENSTtoXXX.head(n=15)
dfENSTtoXXX.sort_values(by='Gene_ID', axis=0, ascending=True, inplace=True, kind='quicksort', na_position='last')
#dfENSTtoXXX.head(n=1)
#dfENSTtoXXX.tail(n=15)
#dfENSTtoXXX.columns.values # colnames
################################################################################

# DF copies

df1 = dfTCGA.copy(deep=True)
df2 = dfENSTtoXXX.copy(deep=True)



colname2 = 'Transcript_ID'
colname1 = 'Gene_ID'

df1[colname1] = ""

df1[[]] # rows names
df1.columns.values # colnames

for i in range(0,df1.shape[0]):
	# Check if the transcript ID in first col of df1 exists in the col
	# transcrit ID of df2 and save the line of df2 with this value in tmp
	tmp = df2[df2[colname2] == df1.index[i]].copy(deep=True)
	tmp = tmp.reset_index(drop=True)
	# Replace 1 from df1[colname1] = "1" by the gene ID if possible or
	# by NaN if not
	df1.ix[i, 497] = tmp.loc[0,colname1] if not tmp.empty else 'NaN'
	# Evaluating speed method with modulo
	if i%1000 == 0 :
		print(' ######### We are at %i'%(i))



df1.sort_values(by='Gene_ID', axis=0, ascending=True, inplace=True, kind='quicksort', na_position='last')
df1.to_csv('data/TCGAtranscritIDGeneID.tsv', sep='\t')


################################################################################

df3 = dflong.copy(deep=True)

df3[colname1] = ""

df3[[]] # rows names
df3.columns.values # colnames

for i in range(0,df3.shape[0]):
	# Check if the transcript ID in first col of df3 exists in the col
	# transcrit ID of df2 and save the line of df2 with this value in tmp
	tmp = df2[df2[colname2] == df3.index[i]].copy(deep=True)
	tmp = tmp.reset_index(drop=True)
	# Replace 1 from df3[colname1] = "1" by the gene ID if possible or
	# by NaN if not
	df3.ix[i, 106] = tmp.loc[0,colname1] if not tmp.empty else 'NaN'
	# Evaluating speed method with modulo
	if i%1000 == 0 :
		print(' ######### We are at %i'%(i))



df3.sort_values(by='Gene_ID', axis=0, ascending=True, inplace=True, kind='quicksort', na_position='last')
df3.to_csv('data/LONGtranscritIDGeneID.tsv', sep='\t')