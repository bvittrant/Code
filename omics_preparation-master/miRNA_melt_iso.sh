#!/bin/bash -i

echo "###############################################################################"

START_TIME=$SECONDS

echo "Compteur initialization"
((compteur=1))

# Create empty file to receive values of CNV informations
touch data/Transcriptome_Profiling_clean/values_iso.txt
>data/Transcriptome_Profiling_clean/values_iso.txt

for file in data/Transcriptome_Profiling_clean/isoforms/*.tsv
	do
		echo "On file "$compteur

		((i=1))

		while read line; do
 	
		 	#echo '#################################################'
		 	#echo 'We are at line' $i

		 	# Echo the value at line i in column $1 == miRNA_ID $2 == isoform coords
		 	# -v j=$i link to FNR == j which is the line number of the iterated file lnked to the count i
		 	# -v OFS='\t' the output will be separated by tab
		 	# The input is read with tab as IFS
		 	tmp="$(awk  -v j=$i -v OFS='\t' -F '\t' 'FNR == j {print $1,$2}' $file)"

		 	#echo $tmp

		 	# Check if the value $1 miRNA $2 isoform coordsdoesn't exist in the file value
		 	# if not write the line if yes nothing
		 	if grep -q "$tmp" data/Transcriptome_Profiling_clean/values_iso.txt
		 		then
		 			#echo 'Got it'
		 			continue
		 		else
		 			echo "$tmp" >> data/Transcriptome_Profiling_clean/values_iso.txt
		 			#echo "$tmp" 
		 	fi

		 	# Increment the count of the lines
	 		((i++))

		done <$file

	((compteur++))

done;

# Check the that the lines in values.txt are unique

echo 'NUmber of lines in values.txt:'
wc -l data/Transcriptome_Profiling_clean/values_iso.txt
echo 'Number of unique line in values.txt:'
sort data/Transcriptome_Profiling_clean/values_iso.txt | uniq -c | wc -l