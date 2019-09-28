# Checking of the prob size 

d = read.table(file = 'results_files/CNV_45_annotation.tsv' )

# Prepare data
d_cnv = read.table(file = 'results_files/CNV_45_merged_clean_-1_1_row497zero.tsv', sep = '\t', header = T, check.names = F)
d_cnv = as.data.frame(t(d_cnv))
colnames(d_cnv) = as.character(d_cnv['barcode',])
d_cnv = d_cnv[-1,]

d_cnv = as.data.frame(d_cnv)

test = rowSums(as.matrix(d_cnv))
typeof(d_cnv)

# ID conversion for ENSG

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
source("http://bioconductor.org/biocLite.R")
biocLite("IRanges")

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
head(listFilters(ensembl))
