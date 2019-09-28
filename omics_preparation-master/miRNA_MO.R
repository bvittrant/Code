###############################################################################

# Prepare miRNA for Mix Omix for fun

d_mi = read.table(file = 'results_files/TCGA_miRNA_quant_cpm_52.tsv', sep = '\t', header = T)
row.names(d_mi) = d_mi$barcode
d_mi = d_mi[,-c(1,2)]

d_clin = read.table(file = 'clinical/TCGA_clinique_52.tsv', sep = '\t', header = T)
row.names(d_clin) = d_clin$barcode
# 1,6,7,10

d_mo = merge(d_clin[,c(1,6,7,10)], d_mi, by='row.names')
row.names(d_mo) = d_mo$Row.names
d_mo = d_mo[,-c(1,2)]

###############################################################################

# Use MO to predict binary BCR event 

###############################################################################

# Installing and loading Mix Omic

#install.packages('mixOmics')
library('mixOmics')
library('RColorBrewer')

###############################################################################

# Mix Omic

#create the combined data set X
X = data.matrix(d_mo[,4:1884])
dim(X)

# Creating the vector of factor value to class/predict
Y = as.factor(as.character(d_mo$BCR_60))
length(Y)

###############################################################################

# PCA
pca.TCGA = pca(X, ncomp = 10, center = TRUE, scale = F)
pca.TCGA

png(file = 'pictures/pca_miRNA.png')
plot(pca.TCGA, main = 'PCA for miRNA data')
dev.off()

# Choose k
pca.TCGA = pca(X, ncomp = 3, center = TRUE, scale = F)

plot(pca.TCGA, main = "Explained variance for the first 10 principal components")

plotIndiv(pca.TCGA, group = d_mo$BCR_60, ind.names = FALSE, 
          ellipse = FALSE, legend = TRUE, title = 'TCGA, PCA comp 1 - 2')


#this chunk takes ~ 5 min to run
set.seed(32) # for reproducibility of the outputs of this code that performs random cross-validation sampling. To be removed in proper analysis
TCGA.plsda.perf <- plsda(X, Y, ncomp = 3)
# to speed up computation in this example we choose 5 folds repeated 10 times:
perf.plsda <- perf(TCGA.plsda.perf, validation = 'Mfold', folds = 5,
                   progressBar = TRUE, nrepeat = 10)

head(perf.plsda$error.rate)
plot(perf.plsda, overlay = 'measure', sd=TRUE)

TCGA.plsda <- plsda(X, Y, ncomp = 3)

plotIndiv(TCGA.plsda , comp = c(1,2),
          group = d_mo$BCR_60, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'TCGA, PLSDA comp 1 - 2')

plotIndiv(TCGA.plsda , comp = c(1,3),
          group = d_mo$BCR_60, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'TCGA, PLSDA comp 1 - 3')

plotIndiv(TCGA.plsda , comp = c(2,3),
          group = d_mo$BCR_60, ind.names = FALSE, 
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
