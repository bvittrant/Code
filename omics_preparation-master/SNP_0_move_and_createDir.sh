#!/bin/bash -i

echo "###############################################################################"

echo "Setting path to bash to use module at source"
. /usr/share/Modules/init/bash

START_TIME=$SECONDS

echo "Compteur initialization"
((compteur=1))
((compteur2=1))
((compteur3=1))

mkdir data/Simple_Nucleotide_Variation_clean/annotations
mkdir data/Simple_Nucleotide_Variation_clean/vcf
mkdir data/Simple_Nucleotide_Variation_clean/parcel

for file in data/Simple_Nucleotide_Variation/*/*.vcf
	do

		IFS="/" read -a WAY <<< "$file"
		IFS="." read -a NAME <<< "${WAY[3]}"
		IFS="_" read -a NUM <<< "${NAME[0]}"

        # We need to recreate the variable name because we split on _

        read -a reconstruct <<< "${WAY[2]}"_"${NAME[0]}".vcf

        #echo $reconstruct

        if [ ${NAME[0]} == 'annotations' ] 
        	then
        	       cp $file data/Simple_Nucleotide_Variation_clean/annotations/$reconstruct
        	       ((compteur2++))
                else
                       cp $file data/Simple_Nucleotide_Variation_clean/vcf/$reconstruct
                        ((compteur3++)) 
        fi

done;

find data/Simple_Nucleotide_Variation/ -name *.parcel -exec cp {} data/Simple_Nucleotide_Variation_clean/parcel \;

touch data/Simple_Nucleotide_Variation_clean/log.txt
echo 'There are '$compteur2' annotations files' >> data/Simple_Nucleotide_Variation_clean/log.txt
echo 'There are '$compteur2' annotations files'
echo 'There are '$compteur3' mean_segment' >> data/Simple_Nucleotide_Variation_clean/log.txt
echo 'There are '$compteur3' mean_segment'

