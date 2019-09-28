#!/bin/bash -i

echo "###############################################################################"

START_TIME=$SECONDS

tmp='obama'

echo "Compteur initialization"
((compteur=1))

for file in data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned_uniq/*
	do
		IFS="/" read -a WAY <<< "$file"
		IFS="." read -a NAME <<< "${WAY[3]}"
		IFS="_" read -a BARCODE <<< "${NAME[0]}"

		#echo $file
		#echo ${WAY[3]}

		if [ $tmp = "${BARCODE[0]}" ]
			then
				cat $file >> data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned_uniq_merged/"${BARCODE[0]}".vcf
				#echo $tmp
			else
				cat $file > data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned_uniq_merged/"${BARCODE[0]}".vcf
				tmp="${BARCODE[0]}"
				#echo $tmp
				#echo ${BARCODE[0]}
		fi

		#echo $tmp

		((compteur++))

		#if [ $compteur -eq 10 ]
		#then
		#	break
		#fi

	done;

for file in data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned_uniq_merged/*
	do
		IFS="/" read -a WAY <<< "$file"
		sort -u $file | uniq > data/Simple_Nucleotide_Variation_clean/vcf_filtered_headerOut_cut_cleaned_uniq_merged_cleaned/${WAY[3]}
	done;