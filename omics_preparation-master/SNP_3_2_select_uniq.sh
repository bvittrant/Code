#!/bin/bash -i

echo "###############################################################################"

START_TIME=$SECONDS

for file in data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned/*.vcf
	do

		IFS="/" read -a WAY <<< "$file"

		uniq $file > data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned_uniq/"${WAY[3]}"


	#if [ $compteur -eq 3 ]
		#then
			#break
	#fi

done;
