###############################################################################

# Use MO to predict binary BCR event 

###############################################################################

rm(list=ls())

# Installing and loading Mix Omic

#install.packages('mixOmics')
library('mixOmics')
library('RColorBrewer')

###############################################################################
###############################################################################

d_mo = read.table(file = 'results_files/CNV_45_merged_clean_-1_1_row497zero.tsv', sep = '\t', header = T,check.names = F)
d_mo = na.omit(d_mo)
d_mo = d_mo[order(row.names(d_mo)),]

d_mo = d_mo[-1,]

d_clin = read.table(file = 'clinical/TCGA_clinique_45_reduced.tsv', sep = '\t', header = T, check.names = F)
row.names(d_clin) = d_clin$barcode
d_clin = d_clin[order(row.names(d_clin)),]

###############################################################################

# Mix Omic

#create the combined data set X
X = data.matrix(d_mo[,])
dim(X)
typeof(X)
is.na(X)

# Creating the vector of factor value to class/predict
Y = as.factor(d_clin$BCR_60)
length(Y)
typeof(Y)

###############################################################################

#Prelimineray PCA
pca.gene = pca(X, ncomp = 10, center = TRUE, scale = FALSE)
pca.gene
png(file = 'pictures/pca_cnv.png')
plot(pca.gene, main = "PCA for CNV data")
dev.off()

pca.gene = pca(exp(X), ncomp = 10, center = TRUE, scale = FALSE)
pca.gene
png(file = 'pictures/pca_cnv_exp.png')
plot(pca.gene, main = "PCA for CNV data passed to exp")
dev.off()



# ncomp = 2

###############################################################################


# Choose k
pca.TCGA = pca(X, ncomp = 3, center = TRUE, scale = F)

plot(pca.TCGA, main = "Explained variance for the first 10 principal components")

plotIndiv(pca.TCGA, group = Y, ind.names = FALSE, 
          ellipse = FALSE, legend = TRUE, title = 'TCGA, PCA comp 1 - 2')


#this chunk takes ~ 5 min to run
set.seed(32) # for reproducibility of the outputs of this code that performs random cross-validation sampling. To be removed in proper analysis
TCGA.plsda.perf <- plsda(X, Y, ncomp = 10)
# to speed up computation in this example we choose 5 folds repeated 10 times:
perf.plsda <- perf(TCGA.plsda.perf, validation = 'Mfold', folds = 5,
                   progressBar = TRUE, nrepeat = 10)

head(perf.plsda$error.rate)
plot(perf.plsda, overlay = 'measure', sd=TRUE)

TCGA.plsda <- plsda(X, Y, ncomp = 3)

plotIndiv(TCGA.plsda , comp = c(1,2),
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'TCGA, PLSDA comp 1 - 2')

plotIndiv(TCGA.plsda , comp = c(1,3),
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'TCGA, PLSDA comp 1 - 3')

plotIndiv(TCGA.plsda , comp = c(2,3),
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'TCGA, PLSDA comp 2 - 3')

#this chunk takes ~ 6 min to run
set.seed(32) # for reproducibility of the outputs of this code that performs random cross-validation sampling. To be removed in proper analysis
# grid of possible keepX values that will be tested for comp 1 and comp 2
list.keepX <- c(1:10,  seq(20, 200, 10))
# to speed up computation in this example we choose 5 folds repeated 10 times:
tune.splsda.TCGA<- tune.splsda(X, Y, ncomp = 3, validation = 'Mfold', folds = 5, 
                               progressBar = FALSE, dist = 'mahalanobis.dist',
                               test.keepX = list.keepX, nrepeat = 10) #nrepeat 50-100
head(tune.splsda.TCGA$error.rate)
tune.splsda.TCGA$choice.keepX

###############################################################################
png(file = 'pictures/miRNA_TCGA_52_fullG_MO_k1.png')
plot(tune.splsda.TCGA, optimal = TRUE, sd = TRUE)
dev.off()

# optimal number of variables to select on 3 comps:
select.keepX = c(200,9,10) #from tuning step
splsda.TCGA <- splsda(X, Y, ncomp = 3, keepX = select.keepX) 
varnames1 = selectVar(splsda.TCGA, comp = 1)$name
write.table(varnames1,'data/list_miRNA_TCGA_k1.tsv', sep = '\t')
varnames2 = selectVar(splsda.TCGA, comp = 2)$name
write.table(varnames2,'data/list_miRNA_TCGA_k2.tsv', sep = '\t')
###############################################################################


plotLoadings(splsda.TCGA, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(splsda.TCGA, comp = 2, method = 'mean', contrib = 'max')
plotLoadings(splsda.TCGA, comp = 3, method = 'mean', contrib = 'max')

plotIndiv(splsda.TCGA, comp = c(1,2),
          group = d_mo$BCR_60, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'TCGA, sPLSDA comp 1 - 2')


set.seed(32)  
perf.splsda <- perf(splsda.TCGA, folds = 5, validation = "Mfold", 
                    dist = "max.dist", progressBar = FALSE, nrepeat = 10)
# perf.srbct  # lists the different outputs
perf.splsda$error.rate
