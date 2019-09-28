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

dSNP_full = pd.read_table(filepath_or_buffer='results_files/SNP_full.tsv', sep='\t', header=0)
#dtype={"Chromosome": object}

dSNP_info = dSNP_full[['#CHROM', 'POS','REF','ALT','MUTATION_TYPE','GENE','condition']]
dSNP_info.shape
dSNP_info.head(n=5)
dTCGA_full = dSNP_full.drop(['#CHROM', 'POS','REF','ALT','MUTATION_TYPE','GENE','condition'], axis=1)
dTCGA_full.shape
dTCGA_full.head(n=5)

dclin = pd.read_table(filepath_or_buffer='clinical/TCGA_clinique_45_reduced.tsv', sep='\t', header=0)
dclin.shape
dclin.head(n=5)

# Selecting col from dTCGA_full which match patient selected in dclin
list_barcode = dclin['barcode'].tolist()
dTCGA_select = dTCGA_full[list_barcode]
dTCGA_select.shape
dTCGA_select.head(n=5)

# Re-merge data frame in one
dSNP_45 = pd.concat([dSNP_info,dTCGA_select], axis=1)
dSNP_45.to_csv(path_or_buf='results_files/SNP_45.tsv', sep ='\t', index = False, header=True)

# Remove with sum of value =< 1 (Only one 1 in row indeed)
dSNP_45_1 = dSNP_45.loc[(dSNP_45.ix[:,7:].sum(axis=1) > 1 ),:]
dSNP_45_1.ix[:,7:].sum(axis=1).describe()
dSNP_45_1.to_csv(path_or_buf='results_files/SNP_45_position_no_unique.tsv', sep ='\t', index = False, header=True)
mutation_load = dSNP_45.ix[:,7:].sum(axis=0)
mutation_load.to_csv(path='results_files/SNP_count_mutation_by_patient.tsv', sep ='\t', index = true, header=True)
print(dSNP_45.shape)
print(dSNP_45_1.shape)

# focus on gene and mutation type
dSNP_45_gene_type = dSNP_45.drop(['#CHROM', 'POS','REF','ALT','condition'], axis=1)
print(dSNP_45_gene_type.shape)
dSNP_45_gene = dSNP_45.drop(['#CHROM', 'POS','REF','ALT','MUTATION_TYPE','condition'], axis=1)

# Sum number of mutation by gene ( number of mutation by gene and patient at the end)
dSNP_45_gene_sum = dSNP_45_gene_type.reset_index().groupby('GENE').sum()
print(dSNP_45_gene_sum.shape)
dSNP_45_gene_sum.head(n=5)

# Select only row with values over 1
dSNP_45_gene_sum_1 = dSNP_45_gene_sum.loc[(dSNP_45_gene_sum.ix[:,1:].sum(axis=1) > 1 ),:]
print(dSNP_45_gene_sum_1.shape)
dSNP_45_gene_sum_1.head(n=5)
dSNP_45_gene_sum_1.drop(['index'])
dSNP_45_gene_sum_1.to_csv(path_or_buf='results_files/SNP_45_gene.tsv', sep ='\t', index = True, header=True)


# Sum number of mutation by gene and mutation type ( number of mutation by gene and type and patient at the end)
dSNP_45_gene_type_sum = dSNP_45_gene_type.reset_index().groupby(['GENE','MUTATION_TYPE']).sum()
print(dSNP_45_gene_type_sum.shape)
dSNP_45_gene_type_sum.head(n=5)

# Select only row with values over 1
dSNP_45_gene_type_sum_1 = dSNP_45_gene_type_sum.loc[(dSNP_45_gene_type_sum.ix[:,2:].sum(axis=1) > 1 ),:]
print(dSNP_45_gene_type_sum_1.shape)
dSNP_45_gene_type_sum_1.head(n=5)
dSNP_45_gene_type_sum_1.drop(['index'], axis=1)
dSNP_45_gene_type_sum_1.to_csv(path_or_buf='results_files/SNP_45_gene_type.tsv', sep ='\t', index = True, header=True)


# for gene alone
dSNP_45_gene_sum_1_clean = dSNP_45_gene_sum_1.ix[(dSNP_45_gene_sum_1 == 0).astype(int).sum(axis=1) < 44]
print(dSNP_45_gene_sum_1_clean.shape)
dSNP_45_gene_sum_1_clean.head(n=5)
dSNP_45_gene_sum_1_clean.to_csv(path_or_buf='results_files/SNP_45_gene_clean.tsv', sep ='\t', index = True, header=True)
dSNP_45_gene_sum_1_clean = dSNP_45_gene_sum_1_clean.reset_index()
dSNP_45_gene_sum_1_clean = dSNP_45_gene_sum_1_clean.drop(['index'], axis=1)

# for gene and mutation type
dSNP_45_gene_type_sum_1_clean = dSNP_45_gene_type_sum_1.ix[(dSNP_45_gene_type_sum_1 == 0).astype(int).sum(axis=1) < 44]
print(dSNP_45_gene_type_sum_1_clean.shape)
dSNP_45_gene_type_sum_1_clean.head(n=5)
dSNP_45_gene_type_sum_1_clean = dSNP_45_gene_type_sum_1_clean.reset_index()
dSNP_45_gene_type_sum_1_clean = dSNP_45_gene_type_sum_1_clean.drop(['index'], axis=1)
dSNP_45_gene_type_sum_1_clean.to_csv(path_or_buf='results_files/SNP_45_gene_type_clean.tsv', sep ='\t', index = False, header=True)

# For gene
list_number = [dSNP_45_gene_sum_1_clean.sum(axis=1)]
list_gene = [dSNP_45_gene_sum_1_clean['GENE']]

df_list_number =pd.DataFrame(list_number)
df_list_gene =pd.DataFrame(list_gene)
df_mutation_per_gene = pd.concat([df_list_number,df_list_gene])
df_mutation_per_gene.shape
df_mutation_per_gene.head(n=5)
df_mutation_per_gene = df_mutation_per_gene.transpose()
df_mutation_per_gene.columns = ['COUNT','GENE']
df_mutation_per_gene.to_csv(path_or_buf='results_files/SNP_mutation_per_gene.tsv', sep ='\t', index = True, header=True)

# For mutation
list_number = [dSNP_45_gene_type_sum_1_clean.sum(axis=1)]
list_gene = [dSNP_45_gene_type_sum_1_clean['GENE']]
list_mutation = [dSNP_45_gene_type_sum_1_clean['MUTATION_TYPE']]

df_list_number =pd.DataFrame(list_number)
df_list_gene =pd.DataFrame(list_gene)
df_list_mutation=pd.DataFrame(list_mutation)

df_mutation_per_gene_type = pd.concat([df_list_gene,df_list_mutation,df_list_number])
df_mutation_per_gene_type.shape
df_mutation_per_gene_type.head(n=5)
df_mutation_per_gene_type = df_mutation_per_gene_type.transpose()
df_mutation_per_gene_type.columns = ['GENE','MUTATION_TYPE','COUNT',]
df_mutation_per_gene_type.to_csv(path_or_buf='results_files/SNP_mutation_per_gene_type.tsv', sep ='\t', index = False, header=True)
