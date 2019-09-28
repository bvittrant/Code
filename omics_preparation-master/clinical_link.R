d1 = read.table(file = 'results_files/TCGA_RNAseq_IGS.tsv', sep = '\t', header = T, check.names = F)
row.names(d1) = d1$file_id

d1$barcode = toupper(d1$barcode)
write.table(d1, file = 'results_files/TCGA_RNAseq_IGS.tsv', sep = '\t', row.names = F)

d7 = read.table(file = 'results_files/TCGA_RNAseq_IGS.tsv', sep = '\t', header = T, check.names = F )