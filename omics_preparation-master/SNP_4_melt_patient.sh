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

		touch tmp.tsv
		> tmp.tsv

		IFS="/" read -a WAY <<< "$file"
		IFS="_" read -a NAME <<< "${WAY[4]}"
		echo ${NAME[0]} >> tmp.tsv

		# Begin at line 2 we don4t care about line 1 header
		((i=2))

		while read line; do
 	
		 	#echo '#################################################'
		 	#echo 'We are at line' $i

		 	# Echo the value at line i in column CHROM	POS	ID	REF	ALT	QUAL	INFO
		 	# -v j=$i link to FNR == j which is the line number of the iterated file lnked to the count i
		 	# -v OFS='\t' the output will be separated by tab
		 	# The input is read with tab as IFS
		 	tmp="$(awk  -v j=$i -v OFS='\t' -F '\t' 'FNR == j {print $1,$2,$3,$4}' $file)"

		 	#echo $tmp
		 	# Check if the value chromosome$2 start$3 end$4 doesn't exist in the file value
		 	# if not write the line if yes nothing
		 	if grep -q "$tmp" $file
		 		then
		 			echo 1 >> tmp.tsv
		 		else
		 			echo 0 >> tmp.tsv
		 	fi

		 	# Increment the count of the lines
	 		((i++))

		done <data/Simple_Nucleotide_Variation_clean/index.txt

	((compteur++))

done;

# Check the that the lines in values.txt are unique

echo 'NUmber of lines in values.txt:'
wc -l data/Simple_Nucleotide_Variation_clean/index.txt
echo 'Number of unique line in values.txt:'
sort data/Simple_Nucleotide_Variation_clean/index.txt | uniq -c | wc -l