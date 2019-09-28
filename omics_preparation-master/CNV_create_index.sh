#!/bin/bash -i

echo "###############################################################################"

START_TIME=$SECONDS

echo "Compteur initialization"
((compteur=1))

# Create empty file to receive values of CNV informations
touch data/Copy_Number_Variation_clean/values.txt
>data/Copy_Number_Variation_clean/values.txt

for file in data/Copy_Number_Variation_clean/mean_segment/*.txt
	do
		echo "On file "$compteur

		((i=1))

		while read line; do
 	
		 	#echo '#################################################'
		 	#echo 'We are at line' $i

		 	# Echo the value at line i in column $2 == File UUID $5 == file Name $8 == Data category
		 	# -v j=$i link to FNR == j which is the line number of the iterated file lnked to the count i
		 	# -v OFS='\t' the output will be separated by tab
		 	# The input is read with tab as IFS
		 	tmp="$(awk  -v j=$i -v OFS='\t' -F '\t' 'FNR == j {print $2,$3,$4,$5}' $file)"

		 	#echo $tmp

		 	# Check if the value chromosome$2 start$3 end$4 doesn't exist in the file value
		 	# if not write the line if yes nothing
		 	if grep -q "$tmp" data/Copy_Number_Variation_clean/values.txt
		 		then
		 			:
		 		else
		 			echo "$tmp" >> data/Copy_Number_Variation_clean/values.txt
		 	fi

		 	# Increment the count of the lines
	 		((i++))

		done <$file

	((compteur++))

done;

# Check the that the lines in values.txt are unique

echo 'NUmber of lines in values.txt:'
wc -l data/Copy_Number_Variation_clean/values.txt
echo 'Number of unique line in values.txt:'
sort data/Copy_Number_Variation_clean/values.txt | uniq -c | wc -l