########################## START ##############################################

# Tentative d'annotation

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
biocLite("clusterProfiler")
biocLite("org.Hs.eg.db")

require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(clusterProfiler)
library(org.Hs.eg.db)

# explo CNV #

d_IGS = read.table(file = '../../data/IGS.tsv', sep = '\t', header = T, check.names = F)

# Prepare data
d_cnv = read.table(file = '../../results_files/CNV_45_merged_clean_-1_1_row497zero.tsv', sep = '\t', header = T, check.names = F)
row.names(d_cnv) = d_cnv$barcode
d_cnv = d_cnv[,-1]
d_cnv = as.data.frame(t(d_cnv))
#colnames(d_cnv) = as.character(d_cnv['barcode',])
row.names(d_cnv)
d_cnv_2 = as.data.frame(matrix(NA, ncol = 17, nrow = length(row.names(d_cnv))))
colnames(d_cnv_2)[1:3] = c('CHR','START','END')

# Assign annotation
for(i in 1:length(row.names(d_cnv))){
  tmp = unlist(strsplit(row.names(d_cnv)[i], '_'))
  d_cnv_2[i,1] = paste('chr', tmp[1],sep='')
  d_cnv_2[i,2] = as.numeric(tmp[2])
  d_cnv_2[i,3] = as.numeric(tmp[3])
  
  #df_tmp = as.data.frame(peakAnno, row.names = NULL, optional = FALSE)
  
  peak <- GRanges(seqnames = d_cnv_2[i,1], ranges = IRanges(d_cnv_2[i,2], d_cnv_2[i,3]))
  peakAnno <- annotatePeak(peak,tssRegion=c(-3000,3000),TxDb = txdb, level = "transcript",
                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"))
  df_tmp = as.data.frame(peakAnno, row.names = NULL, optional = FALSE)
  d_cnv_2[i,4:17] = df_tmp[1,]
  
  if(i == 1){colnames(d_cnv_2)[4:17] = colnames(df_tmp)}
  
}

## Add col gene name et ENSG ID

# Add gene name
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library("biomaRt")

## List all data available
listMarts()
ensembl=useMart("ensembl")
#listDatasets(ensembl)

## Define the one we want

#ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

## Take ID references convertion
HGNC_names_final_full = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','entrezgene','geneID'),
                              filters = 'entrezgene', 
                              values = d_cnv_2$geneId, #d_final_bilan$gene
                              mart = ensembl)
## Change name, care here to do that in the right way in future

hgnc_symbol = as.data.frame(matrix(nrow = dim(d_cnv_2)[1], ncol = 3)) 
colnames(hgnc_symbol) = colnames(HGNC_names_final_full)
d_cnv_3 = cbind(d_cnv_2, hgnc_symbol)

for(i in 1:dim(d_cnv_3)[1]){
  
  tmp = HGNC_names_final_full[which(HGNC_names_final_full$entrezgene == d_cnv_3[i,15]),]
  
  if(!isEmpty(tmp)[1] && !isEmpty(tmp)[2] && !isEmpty(tmp)[3]){
    d_cnv_3[i,18:20] = HGNC_names_final_full[which(HGNC_names_final_full$entrezgene == d_cnv_3[i,15]),]
  }
}

#write.table(d_cnv_3, file = '../../results_files/CNV_45_annotation.tsv',sep = '\t', row.names = F)

length(unique(d_cnv_3$entrezgene))
length(unique(d_cnv_3$geneId))
d_cnv_4 = na.omit(d_cnv_3)

d_cnv_5 = d_cnv_4[which(d_cnv_4$ensembl_gene_id %in% d_IGS$Gene.ID),]

#write.table(d_cnv_5, file = '../../results_files/CNV_45_annotation_IGS.tsv',sep = '\t', row.names = F)

tmp = c()
for(i in 1:dim(d_cnv_3)[1]){
  
  if(d_cnv_3[i,10] == 23){tmp_0 = 'X'}
  if(d_cnv_3[i,10] == 24){tmp_0 = 'Y'}
  if(d_cnv_3[i,10] != 24 && d_cnv_3[i,10] != 23){tmp_0 = d_cnv_3[i,10]}
  
  tmp = c(tmp, paste(tmp_0 , d_cnv_3[i,2], d_cnv_3[i,3]  , sep = '_'))
  
}
tmp[1:5]

d_cnv_6 = cbind(tmp,d_cnv_3)
d_cnv_7 = d_cnv_6[which(d_cnv_6$ensembl_gene_id %in% d_IGS$Gene.ID),]

#write.table(d_cnv_6, file = '../../results_files/CNV_45_annotation.tsv',sep = '\t', row.names = T)
#write.table(d_cnv_7, file = '../../results_files/CNV_45_annotation_IGS.tsv',sep = '\t', row.names = T)

d_cnv_8 = d_cnv[which(row.names(d_cnv) %in% d_cnv_7$tmp),]

#write.table(d_cnv_8, file = '../../results_files/CNV_45_merged_clean_-1_1_row497zero_IGS.tsv',sep = '\t', row.names = T)

########################### END ###############################################
