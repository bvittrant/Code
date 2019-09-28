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

print("################################################################################")
print("Your working directory is " + os.getcwd())
print("Your data frame has %i lines and %i columns" % (df2.shape[0], df2.shape[1]) )
print("################################################################################")
print("The first column has %i unique ID, the second has %i unique ID" 
            % (len(df2.Gene_list.unique()), len(df2.HKG.unique())) )

# Display a specific value in df df2.iloc[0][0]

# Create column to receive the boolean if the gene ID exist

df2["exist_in_HKG"] = ""

################################################################################

# Check if a gene in C1 exist in C2 and add 1 in C3 if true

for i in range(0,df2.shape[0]):
    if df2.iloc[i][0] in df2.HKG.values:
        df2.iloc[i][2] = 1
        print(" %s exists in our reference and is an HKG gene" % (df2.iloc[i][0]) )
    else:
        df2.iloc[i][2] = 0


print("I found %i T-HKG in our reference" % (df2.loc[df2['exist_in_HKG'] == 1, 'exist_in_HKG'].sum()))

df3 = df2[df2['exist_in_HKG'] == 1]

################################################################################

#

# Test to see how ;any HKG we still have

compt_tmp = 0

for i in range(0,dfTCGA_clean_lt30.shape[0]):
    if dfTCGA_clean_lt30.index[i] in df2.HKG.values:
        compt_tmp  = compt_tmp + 1
        print(" %s exists in our reference and is an HKG gene" % (dfTCGA_clean_lt30.index[i]) )
        

print('%i' % (compt_tmp))
