######### load data
# create data with A_cleanRNAseqFileToOnlyTumorProfiles.java before, to localize the housekeeping genes
?
PRADGenes <- 
  read.table("E:/cloud/Data/Long/GSE54460_FPKM-genes-TopHat2-106samples-12-4-13_gene_expression.TPM_geneInfos.JAVA.txt",
             header=TRUE, sep="\t", na.strings="", dec=".", strip.white=TRUE, check.names=FALSE)

######### normalize with log2
hist(unlist(PRADGenes))
a <- log2(PRADGenes[,3:ncol(PRADGenes)-1]+1)
hist(unlist(a))
PRADGeneslog2 <- cbind(PRADGenes[,1], a)
colnames(PRADGeneslog2)[1]<-"patient"
write.table(PRADGeneslog2, file="E:/OneDrive/ulaval/TCGA/gdac.broadinstitute.org.PRAD.illuminahiseq_rnaseqv2_Level_3_RSEM_genes_normalized_transposed_onlyPrimaryTumor.log2.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

######### normalize with DESeq
.libPaths("win-library")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
DESeq

######### normalize with TMM
biocLite("edgeR")

######### normalize with RUVSeq
biocLite("RUVSeq")
library(RUVSeq)

PRAD = setNames(data.frame(t(PRADGenes[,-1])), PRADGenes[,1])#transpose table to have genes in rows
filter <- apply(PRAD, 1, function(x) length(x[x>0])>=5)#remove non-expressed genes, by requiring more than 5 reads in at least two samples for each gene 
filtered <- PRAD[filter,]
genes <- rownames(filtered)[grep("gene_", rownames(filtered))] #set non-housekeeping genes
spikes <- rownames(filtered)[grep("hsk_", rownames(filtered))] #set housekeeping genes
patients <- as.factor(colnames(PRAD)) #patient list

# function to plot the distribution of the normalization
library(RColorBrewer)
plot = function(set,title){
  pdf(paste0("PRADexpression_",title,".pdf"), width=100, height = 20, onefile = T)
  plotRLE(set, outline=FALSE, ylim=c(-2.5, 2.5), main=title, xlab="Patients", las=2)
  dev.off()
}

set0 <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(patients, row.names=colnames(filtered)))
plot(set0,"No_normalization")
set_UQ <- betweenLaneNormalization(set0, which="upper")
plot(set_UQ, "Upper-quartile_normalization")
set_RUVg <- RUVg(set_UQ, spikes, k=1, isLog=F)
plot(set_RUVg, "RUVg_normalization")


#export data in log2
PRAD_ruvg=as.data.frame(normCounts(set_RUVg))
PRAD_ruvgT=as.data.frame(t(PRAD_ruvg))
hist(unlist(PRAD_ruvgT))
PRAD_log2_ruvgT <- log2(PRAD_ruvgT+1)
hist(unlist(PRAD_log2_ruvgT))
colnames(PRAD_log2_ruvgT)=gsub(".*_","",colnames(PRAD_log2_ruvgT))
PRAD_log2_ruvgT = cbind(rownames(PRAD_log2_ruvgT),PRAD_log2_ruvgT)
colnames(PRAD_log2_ruvgT)[1]="Patient"
write.table(PRAD_log2_ruvgT, file="E:/cloud/Data/TCGA_PRAD/gdac.broadinstitute.org.PRAD.illuminahiseq_rnaseqv2_Level_3_RSEM_genes_normalized_transposed_onlyPrimaryTumor.log2RUVg.txt",
            sep="\t", row.names=F, quote=FALSE)
write.table(PRAD_log2_ruvgT, file="E:/cloud/Data/Long/GSE54460_FPKM-genes-TopHat2-106samples-12-4-13_gene_expression.TPM_geneInfos.log2RUVg.txt",
            sep="\t", row.names=F, quote=FALSE)

######### normalize with Z-score
#convert first column as row names
rownames(PRADGenes)=PRADGenes[,1] #set first column as row names
PRADGenes=PRADGenes[,-1] #remove first column that is useless now

#Zscore
PRAD_zscore=as.data.frame(scale(PRADGenes))
PRAD_zscoreLog2=log2(PRAD_zscore+1)
write.table(PRAD_zscore, file="E:/cloud/Data/TCGA_PRAD/gdac.broadinstitute.org.PRAD.illuminahiseq_rnaseqv2_Level_3_RSEM_genes_normalized_transposed_onlyPrimaryTumor.Zscore.txt",
            sep="\t", row.names=T, quote=FALSE, na ="")
write.table(PRAD_zscore, file="E:/cloud/Data/Long/GSE54460_FPKM-genes-TopHat2-106samples-12-4-13_gene_expression.TPM.Zscore.txt",
            sep="\t", row.names=T, quote=FALSE, na ="")

######### normalize with log2 filtered table (removing non expressed data)
log2filtered = t(log2(filtered+1))
colnames(log2filtered)=gsub(".*_","",colnames(log2filtered))
log2filtered = cbind(rownames(log2filtered),log2filtered)
colnames(log2filtered)[1]="Patient"
write.table(log2filtered, file="E:/OneDrive/ulaval/TCGA/gdac.broadinstitute.org.PRAD.illuminahiseq_rnaseqv2_Level_3_RSEM_genes_normalized_transposed_onlyPrimaryTumor.filterNonExpressed_log2.txt",
            sep="\t", row.names=F, quote=FALSE)
