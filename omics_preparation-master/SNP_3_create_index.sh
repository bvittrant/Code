#!/bin/bash -i

echo "###############################################################################"

START_TIME=$SECONDS

echo "Compteur initialization"
((compteur=1))

# Create empty file to receive values of CNV informations
touch data/Simple_Nucleotide_Variation_clean/index.txt
>data/Simple_Nucleotide_Variation_clean/index.txt

for file in data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut/*.vcf
	do
		echo "On file "$compteur

		((i=1))

		IFS="/" read -a WAY <<< "$file"

		touch data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned/"${WAY[3]}"
		>data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned/"${WAY[3]}"

		while read line; do
 	
		 	#echo '#################################################'
		 	#echo 'We are at line' $i

		 	# Echo the value at line i in column CHROM	POS	ID	REF	ALT	QUAL	INFO
		 	# -v j=$i link to FNR == j which is the line number of the iterated file lnked to the count i
		 	# -v OFS='\t' the output will be separated by tab
		 	# The input is read with tab as IFS
		 	tmp="$(awk  -v j=$i -v OFS='\t' -F '\t' 'FNR == j {print $1,$2,$3,$4}' $file)"

		 	INFO="$(awk  -v j=$i -v OFS='\t' -F '\t' 'FNR == j {print $5}' $file)"
 			# Parsing INFO field
 			echo $INFO | grep -E -o "[A-Z,a-z]*_variant\|MOD[A-Z]*\|[a-z,A-Z,0-9]*\|ENSG[0-9]{11}" | uniq | cut -d'|' -f1,3 | uniq > tmp.txt

		 	((k=1))

		 	# ON lit le fichier tmp qui comprends tous les possibles uniques de la categorie INFO de la ligne en lecture du fichier $file
		 	while read line_tmp ; do

		 		mutation_GeneName="$(awk  -v j=$k -v OFS='\t' -F '|' 'FNR == j {print $1,$2}' tmp.txt)"
		 		#echo $mutation_GeneName

		 		full=$tmp$'\t'$mutation_GeneName
		 		echo $full >> data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned/"${WAY[3]}"

			 	if grep -q "$full" data/Simple_Nucleotide_Variation_clean/index.txt
			 		then
			 			:
			 		else
			 			echo "$full" >> data/Simple_Nucleotide_Variation_clean/index.txt
			 	fi

			done <tmp.txt

		 	# Increment the count of the lines
	 		((i++))

		done <$file

	((compteur++))

	#if [ $compteur -eq 3 ]
		#then
			#break
	#fi

done;

# Check the that the lines in values.txt are unique

sed -r -i "s/\ /\t/g" data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned/*.vcf

echo 'NUmber of lines in index.txt:'
wc -l data/Simple_Nucleotide_Variation_clean/index.txt
echo 'Number of unique line in index.txt:'
sort data/Simple_Nucleotide_Variation_clean/index.txt | uniq -c | wc -l