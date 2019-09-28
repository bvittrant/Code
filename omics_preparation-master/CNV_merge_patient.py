################ Preparing CNV set ##################

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

d_mean_CNV = pd.read_table(filepath_or_buffer='/home/benjamin/Desktop/professionnal/\
full_TCGA_data/TCGA_other_data/data/Copy_Number_Variation_clean/CNV_values_sorted.tsv', sep='\t', header=0\
, dtype={"Chromosome": object})

print('\nYour data frame has this structure :')
print(d_mean_CNV.shape)
print('\nand looks like\n')
print(d_mean_CNV.head(n=5))

# Add a colum wich will be the fusion of the tree first (for simply futur test conditions)
# use apply to apply join on each value in row to fill the new columns condition

d_mean_CNV['condition'] = d_mean_CNV.apply(lambda x: '_'.join(x.dropna().astype(str).values), axis=1)
print('Now looks like\n')
print(d_mean_CNV.head(n=5))



#####################################################################################

# We're going to iterate on the CNV file and file d_mean_CNV with the right value

dir_w = cwd+'/data/Copy_Number_Variation_clean/mean_segment/'
count = 1
print('################################################################################################')
print('\nI initiate the count at '+str(count)+' and copy your workind data frame to protect it')
print('################################################################################################')
d_mean_CNV_copy = d_mean_CNV.copy(deep=True)
#d_mean_CNV_copy.head(n=5)
#########################################################################################
for file in os.listdir(dir_w):
    filename = os.fsdecode(file)
    #########################################################################################
    if filename.endswith(".tsv"):
        print('#########################################################################################\n')
        print('Working on file :\n'+os.path.join(cwd, filename))
        print('This is the file number '+str(count))
        # Reading file as a data frame
        d_tmp = pd.read_table(filepath_or_buffer=dir_w+file, sep ='\t', header=0)
        #d_tmp.head(n=5)
        d_tmp['condition'] = d_tmp[['Chromosome','Start','End']].apply(lambda x: '_'.join(x.dropna().astype(str).values), axis=1)
        #d_tmp.head(n=5)
        # Create empty data frame to receive the values of the file iterated
        d_tmp_e = pd.DataFrame(0,index=range(0,d_mean_CNV.shape[0]), columns=[str(filename)])
        #d_tmp_e.head(n=5)
        #########################################################################################
        for index in range(d_tmp.shape[0]):
            if d_tmp.iloc[index]['condition'] in d_mean_CNV['condition'].values:
                #print('Value exists')
                idx = d_mean_CNV[d_mean_CNV['condition'] == d_tmp.iloc[index]['condition']].index.tolist()
                #print(idx)
                d_tmp_e.loc[idx,str(filename)] = d_tmp.iloc[index]['Segment_Mean']
        #########################################################################################
        #print(d_tmp_e.head(n=5))
        count += 1
        d_mean_CNV_copy = pd.concat([d_mean_CNV_copy,d_tmp_e], axis=1)
        #print(d_mean_CNV_copy.iloc[[count]])
        #Free memory to avoid fuite
        #print(d_mean_CNV_copy.head(n=5))
        #print(d_tmp_e.head(n=5))
        #del d_tmp
        #del d_tmp_e
        if count == 5:
            break
d_mean_CNV_copy.to_csv(path_or_buf='results_files/CNV_full.tsv', sep ='\t', index = False, header=True)
#########################################################################################













