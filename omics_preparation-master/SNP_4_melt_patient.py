################ Preparing SNP set ##################

#####################################################################################

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

print('################################################################################################')
print('\nSetting WD\n')
print('################################################################################################')

print('\nYour working directory is :\n'+os.getcwd()) # Prints the working directory

# Set working directory
os.chdir('/home/benjamin/Desktop/professionnal/full_TCGA_data/TCGA_other_data/')
cwd = os.getcwd()

print('Now your working directory is :\n'+os.getcwd()) # Prints the working directory

#####################################################################################

# DL data

dSNP = pd.read_table(filepath_or_buffer='/home/benjamin/Desktop/professionnal/\
full_TCGA_data/TCGA_other_data/data/Simple_Nucleotide_Variation_clean/index.txt', sep='\t', header=None)
#dtype={"Chromosome": object}

print('\nYour data frame has this structure :')
print(dSNP.shape)
print('\nand looks like\n')
print(dSNP.head(n=5))
print('\nI rename columns\n')
dSNP.columns = ['#CHROM', 'POS','REF','ALT','MUTATION_TYPE','GENE']
print(dSNP.head(n=5))

#removed column 5 et 6

#dSNP.drop('Unnamed: 5', axis=1, inplace=True)
#dSNP.drop('Unnamed: 6', axis=1, inplace=True)



# Add a colum wich will be the fusion of the tree first (for simply futur test conditions)
# use apply to apply join on each value in row to fill the new columns condition

dSNP['condition'] = dSNP.apply(lambda x: '_'.join(x.dropna().astype(str).values), axis=1)
print('Now looks like\n')
print(dSNP.head(n=5))



#####################################################################################

# We're going to iterate on the SNP file and file dSNP with the right value

dir_w = cwd+'/data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned_uniq_merged_cleaned/'
count = 1
print('################################################################################################')
print('\nI initiate the count at '+str(count)+' and copy your workind data frame to protect it')
print('################################################################################################')
dSNP_copy = dSNP.copy(deep=True)
#dSNP_copy.head(n=5)
print('Iterating on your files')
#########################################################################################
for file in tqdm(os.listdir(dir_w)):
    filename = os.fsdecode(file)
    #########################################################################################
    if filename.endswith(".vcf"):
        # Reading file as a data frame
        d_tmp = pd.read_table(filepath_or_buffer=dir_w+file, sep ='\t', header=None)
        #d_tmp.head(n=5)
        d_tmp.columns = ['#CHROM', 'POS','REF','ALT','MUTATION_TYPE','GENE']
        d_tmp['condition'] = d_tmp[['#CHROM','POS','REF','ALT','MUTATION_TYPE','GENE']].apply(lambda x: '_'.join(x.dropna().astype(str).values), axis=1)
        #d_tmp.head(n=5)
        # Create empty data frame to receive the values of the file iterated
        d_tmp_e = pd.DataFrame(0,index=range(0,dSNP.shape[0]), columns=[str(filename).split(".")[0]])
        #d_tmp_e.head(n=5)
        #########################################################################################
        for index in range(d_tmp.shape[0]):
            if d_tmp.iloc[index]['condition'] in dSNP['condition'].values:
                #print('Value exists')
                idx = dSNP[dSNP['condition'] == d_tmp.iloc[index]['condition']].index.tolist()
                #print(idx)
                d_tmp_e.loc[idx,str(filename).split(".")[0]] = 1 # put 1 for presence of the mutation
        #########################################################################################
        #print(d_tmp_e.head(n=5))
        count += 1
        dSNP_copy = pd.concat([dSNP_copy,d_tmp_e], axis=1)
        #print(dSNP_copy.iloc[[count]])
        #Free memory to avoid fuite
        #print(dSNP_copy.head(n=5))
        #print(d_tmp_e.head(n=5))
        #del d_tmp
        #del d_tmp_e
dSNP_copy.to_csv(path_or_buf='results_files/SNP_full.tsv', sep ='\t', index = False, header=True)
#########################################################################################













