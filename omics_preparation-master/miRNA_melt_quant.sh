#!/bin/bash 

echo "###############################################################################"
echo "Compteur initialization"
((compteur=1))

mkdir results_files
	mkdir results_files/final_files
		mkdir results_files/final_files/Cut

for file in data/Transcriptome_Profiling_clean/quantification/*.tsv
	do

		IFS="/" read -a WAY <<< "$file"

		cut -f 3 $file > results_files/final_files/Cut/${WAY[3]}
		chmod 744 results_files/final_files/Cut/${WAY[3]}
		#sort -k 1 -o TCGA/Cut/${WAY[2]} TCGA/Cut/${WAY[2]}


		# Creating result file on first iteration
		if [ ${compteur} = '1' ]
			then
				#echo "Compteur = 1"
				cut -f 1,3 $file > results_files/results_${compteur}.tsv
				chmod 744 results_files/results_${compteur}.tsv
				#sort -k 1 -o $file TCGA_results.txt
		fi

		# Join file but not on the first iteration
		if [ ${compteur} -ge 2 ]
			then
			#echo "Compteur >= 2"
				num=$((compteur-1))
				paste results_files/results_${num}.tsv results_files/final_files/Cut/${WAY[3]} > results_files/results_${compteur}.tsv
				rm results_files/results_${num}.tsv
		fi


		((compteur++))
	done;

((compteur--))

mv results_files/results_${compteur}.tsv results_files/TCGA_miRNA_cpm.tsv

rm -r results_files/final_files/Cut
rm -r results_files/final_files