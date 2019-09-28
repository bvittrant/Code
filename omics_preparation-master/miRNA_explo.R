# Exploration of the miRNA values

d1 = read.table(file = 'results_files/TCGA_miRNA_quant_cpm_52.tsv', header = T, check.names = F, sep = '\t')
row.names(d1) = d1$barcode
d1 = d1[,-c(1,2)]

# Remove columns with only 0 values
d2 = d1[, colSums(d1 != 0) > 0]
# Remove columns with values inferior to 1
d3 = d1[, colSums(d1 < 10 ) < 52]
# Remove columns with all values inferior to 10
d4 = d1[, colSums(d1 < 10 ) < 52]
# Remove columns with all values inferior to 20
d5 = d1[, colSums(d1 < 20 ) < 52]
# Remove columns with all values inferior to 100
d6 = d1[, colSums(d1 < 100 ) < 52]

write.table(d2, file = 'results_files/TCGA_miRNA_quant_cpm_52_threshold10.tsv', sep = '\t', row.names = T, quote = F)

summary(d1)
plot(d1)

d = read.table(file = 'results_files/TCGA_miRNA_iso_cpm.tsv', header = T,)
d1 = read.table(file = 'results_files/TCGA_miRNA_cpm_52.tsv', header = T)

join -e 'NA' --header -t $'\t' -21 -22 ffea2175-e304-4c56-b5bc-dde0cef3240a_isoforms_quantification.tsv ffec771e-e2d6-40d3-ba07-3932e8bccc0f_isoforms_quantification.tsv
join: incompatible join fields 1, 2

