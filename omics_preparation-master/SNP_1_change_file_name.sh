#!/bin/bash -i

((i=1))

for file in data/Simple_Nucleotide_Variation_clean/vcf_filtered/*.vcf
	do
		tmp=$(grep -E -o "TCGA-.{2}-.{4}" $file | sort -u)
		#echo '######################################'
		#echo $i
		#echo $tmp
		#echo $file

		if [ -e data/Simple_Nucleotide_Variation_clean/vcf_filtered_2/$tmp.vcf ]
			then
				#echo 'exist'
				#echo $i
				#echo $tmp
				#echo ${tmp}_$i.vcf
				cp $file data/Simple_Nucleotide_Variation_clean/vcf_filtered_2/${tmp}_$i.vcf
			else
				#echo 'NO exist'
				cp $file data/Simple_Nucleotide_Variation_clean/vcf_filtered_2/$tmp.vcf
		fi
		((i++))
done;