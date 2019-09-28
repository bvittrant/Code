#!/bin/bash -i

echo "###############################################################################"

echo "Setting path to bash to use module at source"
. /usr/share/Modules/init/bash

START_TIME=$SECONDS

echo "Compteur initialization"
((compteur=1))
((compteur2=1))
((compteur3=1))

mkdir data/Copy_Number_Variation_clean/annotations
mkdir data/Copy_Number_Variation_clean/mean_segment
mkdir data/Copy_Number_Variation_clean/parcel

for file in data_raw/Copy_Number_Variation/*/*.txt
	do

		IFS="/" read -a WAY <<< "$file"
		IFS="." read -a NAME <<< "${WAY[3]}"
		IFS="_" read -a NUM <<< "${NAME[0]}"

        # We need to recreate the variable name because we split on _

        read -a reconstruct <<< "${NAME[0]}"_"${NAME[1]}".txt

        #echo $reconstruct

        if [ ${NAME[0]} == 'annotations' ] 
        	then
        	       cp $file data/Copy_Number_Variation_clean/annotations/$reconstruct
        	       ((compteur2++))
                else
                       cp $file data/Copy_Number_Variation_clean/mean_segment/$reconstruct
                        ((compteur3++)) 
        fi

done;

find data/Copy_Number_Variation/ -name *.parcel -exec cp {} data/Copy_Number_Variation_clean/parcel \;

touch data/Copy_Number_Variation_clean/log.txt
echo 'There are '$compteur2' annotations files' >> data/Copy_Number_Variation_clean/log.txt
echo 'There are '$compteur2' annotations files'
echo 'There are '$compteur3' mean_segment' >> data/Copy_Number_Variation_clean/log.txt
echo 'There are '$compteur3' mean_segment'

