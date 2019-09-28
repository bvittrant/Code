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

d_clin = pd.read_table(filepath_or_buffer='clinical/clin_CNV_2.csv', sep='\t', 
	header=0)

print('\nYour clinical data frame has this structure :')
print(d_clin.shape)
#print('\nand looks like\n')
#print(d_clin.head(n=5))

d_cnv = pd.read_table(filepath_or_buffer='results_files/CNV_reduced.tsv', sep='\t', 
	header=0)
print('\nYour CNV data frame has this structure :')
print(d_cnv.shape)

###############################################################################

# Keep only columns of interest (free memory and avoid to swap on swap)
d_clin = d_clin[['cases_0_submitter_id','file_name']]
print('\nYour reduced clinical data frame has this structure :')
print(d_clin.shape)
#print(d_clin.head(n=5))

# Sorting data frame by barcode
d_clin.sort_values(['cases_0_submitter_id'], inplace=True)
#print(d_clin.head(n=5))

# Re indexing data frame (If not index is sorted from 0 to n)
d_clin.reset_index(inplace=True, drop=True)
#print(d_clin.head(n=5))

###############################################################################

# Melt columns
# Iterate on barcode and add in a list all files corresponding on the same
# barcode. Then melt the column corresponding to this file by selecting them
# with the list df[[list]] and melt them in one

# Create empty matrix ready to be append by sum of columns linked by the
# same barcode
# ATTENTION Remember that Python does not slice inclusive of the ending index
d_cnv_barcode = d_cnv.ix[:,0:4]
#print(d_cnv_barcode.shape)
#print(d_cnv_barcode.head(n=5))

print('##############################################################')
print('Starting loop')
print('##############################################################')
print('Go take coffee and come back before progress bar be full')
print('\n')

for i in tqdm(range(d_clin.shape[0])):
	# Condition start
	if i == 0:
		#Initialize list
		list_prob = []
		# Stocking barcode in tmp
		tmp_barcode = d_clin['cases_0_submitter_id'][i]
		# Creating list to receive file names
		list_filename = []
		if d_clin['file_name'][i] in list(d_cnv):
			list_filename.append(d_clin['file_name'][i])
		else:
			list_prob.append([d_clin['file_name'][i]])
	# Condition end
	if i == d_clin.shape[0]-1:
		# Append list with the last and append DF
		if d_clin['file_name'][i] in list(d_cnv):
			list_filename.append(d_clin['file_name'][i])
			d_cnv_barcode[tmp_barcode] = d_cnv[list_filename].sum(axis=1)
		else:
			ist_prob.append([d_clin['file_name'][i]])
			d_cnv_barcode[tmp_barcode] = d_cnv[list_filename].sum(axis=1)
	# Condition in the middle
	else:
		# If barcode are the same on the iteration than the previous
		if tmp_barcode == d_clin['cases_0_submitter_id'][i]:
			# Then append the list with the new filename
			if d_clin['file_name'][i] in list(d_cnv):
				list_filename.append(d_clin['file_name'][i])
			else:
				list_prob.append(d_clin['file_name'][i])
		# If barcode are diffrents
		else:
			# Melt columns on the list in one columns with barcode name
			d_cnv_barcode[tmp_barcode] = d_cnv[list_filename].sum(axis=1)
			# Then reattribute tmp_barcode and create a new list
			tmp_barcode = d_clin['cases_0_submitter_id'][i]
			list_filename = []
			if d_clin['file_name'][i] in list(d_cnv):
				list_filename.append(d_clin['file_name'][i])
			else:
				list_prob.append(d_clin['file_name'][i])



print('\n')
print('##############################################################')
print('Ending loop')
print('##############################################################')

print('Your final DF has this shape')
print(d_cnv_barcode.shape)
#print(d_cnv_barcode.head(n=5))

d_cnv_barcode.to_csv(path_or_buf='results_files/CNV_full_merged.tsv',
	sep ='\t', index = False, header=True)

list_prob.to_csv(,sep ='\t')