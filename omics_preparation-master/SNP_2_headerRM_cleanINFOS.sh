# Remove line begining by ##
# run in vcf_filtered dir
for d in *.vcf ; do grep -v '^##' $d > ../vcf_filtered_headerOut/$d ; done

# Remove the part of infos we don4t want in two step
# CONSERVE ORDER OF THE SED COMMAND

# run in vcf_filtered_headerOut dir
#sed -r -i "s/\|ENSG.*R1\t/\t/g" *.vcf
#sed -r -i "s/\|ENSG.*\t.*/\t/g" *.vcf 
#sed -r -i "s/CSQ.*MOD[A-Z]{5}\|//g" *.vcf
#sed -r -i "s/CSQ.*LOW\|//g" *.vcf
#sed -r -i "s/CSQ.*HIGH\|//g" *.vcf
#sed -r -i "s/\|.*SNV\|*/intergenic_variant/g" *.vcf
#sed -r -i "s/\tintergenic_variant.*/\tintergenic_variant/g" *.vcf
#sed -r -i "s/intergenic_variantE_Multiple_observations/intergenic_variant/g" *.vcf
#sed -r -i "s/intergenic_variantT.*/intergenic_variant/g" *.vcf

#ACGT.*;
#\|[A-Z,a-z]*_variant\|MOD[A-Z]*\|[a-z,A-Z,0-9]*\|ENSG[0-9]{11}\|

#remove fisrt row, prob with header 
# 11 col and rest of the file 9 
#sed -i '1d' *.vcf

#keep col 1,2,3,4,5,6,8 of
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
# then
#CHROM	POS	ID	REF	ALT	QUAL	INFO
for d in *.vcf ; do cut -f1,2,4,5,8 $d > ../vcf_filtered_headerOut_cut/$d ; done
