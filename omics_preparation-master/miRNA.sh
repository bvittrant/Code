#!/bin/bash -i

echo "###############################################################################"

echo "Setting path to bash to use module at source"
. /usr/share/Modules/init/bash

START_TIME=$SECONDS

echo "Compteur initialization"
((compteur=1))
((compteur2=1))
((compteur3=1))

mkdir data/Transcriptome_Profiling_clean/quantification
mkdir data/Transcriptome_Profiling_clean/isoforms
mkdir data/Transcriptome_Profiling_clean/annotation

for file in data/Transcriptome_Profiling/*/*.txt
	do

		IFS="/" read -a WAY <<< "$file"
		IFS="." read -a NAME <<< "${WAY[3]}"
		IFS="_" read -a NUM <<< "${NAME[0]}"

        # We need to recreate the variable name because we split on _

        read -a reconstruct <<< "${WAY[2]}"_"${NAME[0]}"_"${NAME[1]}".tsv

        #echo $reconstruct

        if [ ${NAME[0]} == 'mirnas' ] 
        	then
        		cp $file data/Transcriptome_Profiling_clean/quantification/$reconstruct
        		((compteur++))
        fi
        if [ ${NAME[0]} == 'isoforms' ] 
        	then
        		cp $file data/Transcriptome_Profiling_clean/isoforms/$reconstruct
        		((compteur3++))
        fi
        if [ ${NAME[0]} == 'annotations' ] 
        	then
        		cp $file data/Transcriptome_Profiling_clean/annotation/$reconstruct
        		((compteur2++))
        fi

done;

touch data/Transcriptome_Profiling_clean/log.txt
echo 'There are '$compteur' mirnas quantification files' >> data/Transcriptome_Profiling_clean/log.txt
echo 'There are '$compteur' mirnas quantification files'
echo 'There are '$compteur2' annotations files' >> data/Transcriptome_Profiling_clean/log.txt
echo 'There are '$compteur2' annotations files'
echo 'There are '$compteur3' isoforms files' >> data/Transcriptome_Profiling_clean/log.txt
echo 'There are '$compteur3' isoforms files'

