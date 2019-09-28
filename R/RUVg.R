# RUV normalisation

# Set working directory
setwd("~/Desktop/test")
setwd("C:/Users/benjamin/Desktop/test")

# Loading packages

source("https://bioconductor.org/biocLite.R")
biocLite("RUVSeq")
library(RUVSeq)

# Load raw data 
TCGA_lt9_geneID = na.omit(read.table('data/dfTCGA_geneID_lt9.tsv',header=TRUE, sep="\t",check.names = F))
LONG_lt9_geneID = na.omit(read.table('data/dfLONG_geneID_lt9.tsv',header=TRUE, sep="\t",check.names = F))

row.names(TCGA_lt9_geneID) = TCGA_lt9_geneID$Gene_ID
TCGA_lt9_geneID = TCGA_lt9_geneID[,-1]
row.names(LONG_lt9_geneID) = LONG_lt9_geneID$Gene_ID
LONG_lt9_geneID = LONG_lt9_geneID[,-1]

# Round because of kallisto bootstrap

TCGA_lt9_geneID = round(TCGA_lt9_geneID)
LONG_lt9_geneID = round(LONG_lt9_geneID)

# Define spike in genes

# ENSG00000075624 	ACTB
# ENSG00000111640 	GAPDH
# ENSG00000169919 	GUSB
# ENSG00000196262 	PPIA

spikes = c('ENSG00000075624','ENSG00000111640','ENSG00000169919','ENSG00000196262')
intersect(spikes, rownames(TCGA_lt9_geneID))

# Normalize

#TCGA_raw_geneID_filtered = as.matrix(apply(TCGA_raw_geneID , 1, function(x) length(x[x>5])>=2))

RUVg_TCGA = RUVg(as.matrix(TCGA_lt9_geneID), intersect(spikes, rownames(TCGA_lt9_geneID)), k=1, isLog=F)
RUVg_LONG = RUVg(as.matrix(TCGA_lt9_geneID), intersect(spikes, rownames(TCGA_lt9_geneID)), k=1, isLog=F)

df_RUVg_TCGA = as.data.frame(RUVg_TCGA$normalizedCounts)
df_RUVg_LONG = as.data.frame(RUVg_LONG$normalizedCounts)

# Put all negativ values to 0
df_RUVg_TCGA[df_RUVg_TCGA<0]=0
df_RUVg_LONG[df_RUVg_LONG<0]=0

# Write df in a file and export them
write.table(df_RUVg_TCGA,'data/TCGA_norm_lt9.tsv', sep = '\t')
write.table(t(df_RUVg_TCGA),'data/TCGA_norm_lt9_transpose.tsv', sep = '\t')

write.table(df_RUVg_LONG,'data/LONG_norm_lt9.tsv', sep = '\t')
write.table(t(df_RUVg_LONG),'data/LONG_norm_lt9_transpose.tsv', sep = '\t')
