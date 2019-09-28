###############################################################################

## Load data 

d_annotation = read.table(file = '../../data_final/CNV_45_annotation.tsv', header = T, check.names = F)
d_annotation_IGS = read.table(file = '../../data_final/CNV_45_annotation_IGS.tsv', header = T, check.names = F)

d_cnv = read.table(file = '../../data_final/CNV_45_merged_clean_-1_1_row497zero.tsv', header = T, check.names = F)
row.names(d_cnv) = d_cnv$barcode
d_cnv = d_cnv[,-1]

d_cnv_IGS = read.table(file = '../../data_final/CNV_45_merged_clean_-1_1_row497zero_IGS.tsv', header = T, check.names = F)
row.names(d_cnv_IGS) = d_cnv_IGS$barcode
d_cnv_IGS = d_cnv_IGS[,-1]

d_cnv = as.data.frame(t(d_cnv))
d_cnv_IGS = as.data.frame(t(d_cnv_IGS))

###############################################################################

## data preparation

d_cnv = d_cnv[order(row.names(d_cnv)),]
d_cnv_IGS = d_cnv_IGS[order(row.names(d_cnv_IGS)),]

d_annotation = d_annotation[order(d_annotation$tmp),]
d_annotation_IGS = d_annotation[order(d_annotation_IGS$tmp),]

d_cnv_2 = cbind(d_annotation$ensembl_gene_id, d_cnv)
row.names(d_cnv_2) = d_cnv_2$`d_annotation$ensembl_gene_id`

d_cnv_2_IGS = cbind(d_annotation_IGS$ensembl_gene_id, d_cnv_IGS)
row.names(d_cnv_2_IGS) = d_cnv_2_IGS$`d_annotation_IGS$ensembl_gene_id`

## Melt value of cnv by gene name (melt site in gene)
d_cnv_3 = rowsum(d_cnv_2[,2:ncol(d_cnv_2)] , d_cnv_2$`d_annotation$ensembl_gene_id`, reorder = TRUE)
d_cnv_3_IGS = rowsum(d_cnv_2_IGS[,2:ncol(d_cnv_2_IGS)] , d_cnv_2_IGS$`d_annotation_IGS$ensembl_gene_id`, reorder = TRUE)


tmp = c()
for(i in 1:dim(d_cnv_3)[1]){
  tmp = c(tmp, paste('CNV', row.names(d_cnv_3)[i], sep = '_'))
}
row.names(d_cnv_3) = tmp

tmp = c()
for(i in 1:dim(d_cnv_3_IGS)[1]){
  tmp = c(tmp, paste('CNV', row.names(d_cnv_3_IGS)[i], sep = '_'))
}
row.names(d_cnv_3_IGS) = tmp

write.table(d_cnv_3, file = '../../data_final/CNV_45_merged_clean_-1_1_row497zero_SiteMeltByGene.tsv', row.names = T, sep = '\t')
write.table(d_cnv_3_IGS, file = '../../data_final/CNV_45_merged_clean_-1_1_row497zero_IGS_SiteMeltByGene.tsv', row.names = T, sep = '\t')

###############################################################################

## On check que pour chaque patient les CNV par gene sont de memes signes
tmp_IGS = c()
for (i in 2:dim(d_cnv_2_IGS)[2]) {
  for (j in 2:dim(d_cnv_2_IGS)[1]) {
    a = 0
    b = 0
    if(d_cnv_2_IGS[j,i] == d_cnv_2_IGS[j-1,i]){
      a = d_cnv_2_IGS[j,i]
      b = d_cnv_2_IGS[j-1,i]
    }
    if((a < 0 && b > 0) || (a > 0 && b < 0)){tmp = c(tmp,as.character(d_cnv_2_IGS[j,1]))}
  }
}

tmp = c()
for (i in 2:dim(d_cnv_2)[2]) {
  for (j in 2:dim(d_cnv_2)[1]) {
    a = 0
    b = 0
    if(d_cnv_2[j,i] == d_cnv_2[j-1,i]){
      a = d_cnv_2[j,i]
      b = d_cnv_2[j-1,i]
    }
    if((a < 0 && b > 0) || (a > 0 && b < 0)){tmp = c(tmp,as.character(d_cnv_2[j,1]))}
  }
}

