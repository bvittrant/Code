
## Loading data

d_cnv = read.table(file = '../../data_final/CNV_45_merged_clean_-1_1_row497zero_SiteMeltByGene.tsv', sep = '\t', header = T, check.names = F)
d_cnv = as.data.frame(t(d_cnv))
d_cnv_IGS = read.table(file = '../../data_final/CNV_45_merged_clean_-1_1_row497zero_IGS_SiteMeltByGene.tsv', sep = '\t', header = T, check.names = F)
d_cnv_IGS = as.data.frame(t(d_cnv_IGS))
colnames(d_cnv)

d_rna = read.table(file = '../../data_final//TCGA_RNAseq_45_IfullG.tsv', sep = '\t', header = T, check.names = F)
d_rna_IGS = read.table(file = '../../data_final//TCGA_RNAseq_45_norm_IGS.tsv', sep = '\t', header = T, check.names = F)

d_mirna = read.table(file = '../../data_final//TCGA_miRNA_quant_cpm_45.tsv',sep = '\t', header = T, check.names = F)
#colnames(d_mirna)
#row.names(d_mirna)
row.names(d_mirna) = d_mirna$barcode
d_mirna = d_mirna[,-c(1,2)]

d_meth = read.table(file = '../../data_final/METH_45.tsv', sep = '\t', header = T, check.names = F)
row.names(d_meth) = d_meth$gene
d_meth = d_meth[,-1]
d_meth = as.data.frame(t(d_meth))
d_meth_IGS = read.table(file = '../../data_final/METH_45_IGS.tsv', sep = '\t', header = T, check.names = F)
row.names(d_meth) = d_meth_IGS$patient
d_meth_IGS = d_meth_IGS[,-1]

d_snp = read.table(file = '../../data_final/SNP_45_gene_clean.tsv', sep = '\t', header = T, check.names = F)
#row.names(d_snp)
#colnames(d_snp)
row.names(d_snp) = d_snp$GENE
d_snp = d_snp[,-1]
d_snp = as.data.frame(t(d_snp))
d_snp_IGS = read.table(file = '../../data_final/SNP_45_gene_clean_IGS.tsv', sep = '\t', header = T, check.names = F)

tmp = merge(d_cnv_IGS,d_rna_IGS, by=0, all=T)
tmp1 = merge(tmp, d_mirna,by=0, all=T)

tmp2 = merge(tmp1, d_mirna,by=0, all=T)
tmp3 = merge(tmp2, d_mirna,by=0, all=T)

duplicated(colnames(tmp3))
colnames(tmp3)
