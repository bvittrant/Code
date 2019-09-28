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

dfenst_ensg = pd.read_table('data/TCGAtranscritIDGeneID.tsv', sep='\t', index_col=0)
dfenst_ensg[[]] # rows names
dfenst_ensg.columns.values # colnames

# group and fuse row by gene name

TCGAgeneIDraw = dfenst_ensg.groupby('Gene_ID').sum()
TCGAgeneIDraw[[]]
TCGAgeneIDraw.columns.values

TCGAgeneIDraw.to_csv('data/dfTCGA_geneID_lt9.tsv', sep='\t')

################################################################################

# Import data

dfenst_ensg = pd.read_table('data/LONGtranscritIDGeneID.tsv', sep='\t', index_col=0)
dfenst_ensg[[]] # rows names
dfenst_ensg.columns.values # colnames

# group and fuse row by gene name

LONGgeneIDraw = dfenst_ensg.groupby('Gene_ID').sum()
LONGgeneIDraw[[]]
LONGgeneIDraw.columns.values

LONGgeneIDraw.to_csv('data/dfLONG_geneID_lt9.tsv', sep='\t')