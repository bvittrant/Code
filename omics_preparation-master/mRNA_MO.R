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

d_mo = read.table(file = 'data/MO_final.tsv', sep = '\t', header = T)
row.names(d_mo) = d_mo$patients
d_mo = d_mo[,-1]
d_mo = na.omit(d_mo)

###############################################################################

# Mix Omic

#create the combined data set X
X = data.matrix(d_mo[,3:908])
dim(X)
typeof(X)
is.na(X)

# Creating the vector of factor value to class/predict
Y = as.factor(d_mo$BCR)
length(Y)
typeof(Y)


# The vector indicating each independent study
study = as.factor(d_mo$study)
table(Y,study)
length(study)


###############################################################################

#Prelimineray PCA
pca.gene = pca(X, ncomp = 10, center = TRUE, scale = FALSE)
pca.gene
png(file = 'pictures/PCA_BCR_TCGA.png')
plot(pca.gene, main = "PCA for mRNA data")
dev.off()

# ncomp = 2

###############################################################################

# Chose best keepX
# A faire avec la version beta
tune.mint = tune(X = X, Y = Y, study = study, ncomp = 2, 
                 test.keepX = c(5, 10, 15, 20, 30, 100), method = 'mint.splsda', 
                 progressBar = T, scale = T, near.zero.var = F)

tune.mint
tune.mint$error.rate
plot(tune.mint)
tune.mint$choice.keepX.constraint
tune.mint$choice.keepX

png(file = 'pictures/MO_BER_BCR60_TCGA_LONG.png')
plot(tune.mint)
dev.off()

###############################################################################
###############################################################################

###############################################################################
###############################################################################

# Define colors method 1
jBrewColors1 = brewer.pal(n = length(unique(d_mo[1:141,1])), name = 'YlOrRd')
cond.row1 = jBrewColors1[as.factor(d_mo[1:141,2])]

jBrewColors2 = brewer.pal(n = length(unique(d_mo[1:141,2])), name = 'Blues')
cond.row2 = jBrewColors2[as.factor(d_mo[1:141,2])]

# Define color method 2
colfunc = colorRampPalette(c("red","yellow","springgreen","royalblue"))

###############################################################################
###############################################################################

# sSPLS-DA

mint.splsda.res = mint.splsda(X = X, Y = Y, study = study, ncomp = 2, keepX = c(10,30), scale = T)

varnames1 = selectVar(mint.splsda.res, comp = 1)$name
write.table(varnames1,'data/list_TCGA_LONG_k1.tsv', sep = '\t')
varnames2 = selectVar(mint.splsda.res, comp = 2)$name
write.table(varnames2,'data/list_TCGA_LONG_k2.tsv', sep = '\t')

###############################################################################


plotIndiv(mint.splsda.res, legend = TRUE, title = 'MINT PLS-DA', 
          subtitle = 'No variable selection')
plotIndiv(mint.splsda.res, study = 'global', legend = TRUE, title = 'MINT sPLS-DA', 
          subtitle = 'Global')
plotIndiv(mint.splsda.res, study = 'all.partial', legend = TRUE, 
          title = 'MINT sPLS-DA')
plotArrow(mint.splsda.res, legend = TRUE)
cim(mint.splsda.res)


png(file = 'pictures/HM_BCR60_TCGA_LONG_k2.png')
par(mar=c(9,5,6,1))
cim(mint.splsda.res,
    comp = 2,
    row.sideColors = cbind(cond.row1), # cbind(cond.row1,cond.row2),
    row.names = row.names(d_mo)) #row.names(d_MO_2)
#col.names = varnames2)
legend("topright", legend=sort(unique(d_MO_2$BCR_60)), 
       fill=jBrewColors1,
       title="BCR 60" , box.lty=0, cex = 0.70, ncol=1)
#legend("bottomright", legend=sort(unique(d_MO_2$BCR_50)), 
#fill=jBrewColors2,
#title="BCR 50" , box.lty=0, cex = 0.70, ncol=1)
#legend("bottomright", legend=sort(unique(d_MO_2$BCR_40)), 
#fill=jBrewColors3,
#title="BCR 40" , box.lty=0, cex = 0.70, ncol=1)
#legend(x = 180, y = 1, legend=sort(unique(d_MO_2$grade)), 
#fill=jBrewColors4,
#title="grade" , box.lty=0, cex = 0.70, ncol=1)
#legend(x = 180, y = -7, legend=sort(unique(d_MO_2$gleason)), 
#fill=jBrewColors5,
#title="gleason" , box.lty=0, cex = 0.70, ncol=1)
dev.off()



d1 = read.table(file='results_files/MO_TCGA.tsv', sep= '\t', header = T)
d2 = read.table(file='results_files/TCGA_RNAseq_52_IGS.tsv', sep= '\t', header = T)
d3 = read.table(file='results_files/TCGA_RNAseq_52_fullG.tsv', sep= '\t', header = T)
row.names(d2) = d2$barcode
row.names(d3) = d3$barcode

d2 = d2[intersect(row.names(d1), row.names(d2)),]
d3 = d3[intersect(row.names(d1), row.names(d3)),]

write.table(d2[,-c(2:9)], file = 'results_files/TCGA_RNAseq_45_IGS.tsv', row.names = F, sep = '\t')
write.table(d3[,-c(2)], file = 'results_files/TCGA_RNAseq_45_IfullG.tsv', row.names = F, sep = '\t')

d4 = read.table(file='clinical/TCGA_clinique_52_full.tsv', sep= '\t', header = T)
d5 = read.table(file='clinical/TCGA_clinique_52.tsv', sep= '\t', header = T)
d4$barcode = toupper(d4$barcode)
row.names(d4) = d4$barcode
row.names(d5) = d5$barcode


d4 = d4[intersect(row.names(d1), row.names(d4)),]
d5 = d5[intersect(row.names(d1), row.names(d5)),]

write.table(d4, file = 'clinical/TCGA_clinique_45_full.tsv', row.names = F, sep = '\t')
write.table(d5, file = 'clinical/TCGA_clinique_45_reduced.tsv', row.names = F, sep = '\t')

d3 = read.table(file = 'results_files/TCGA_RNAseq_45_IfullG.tsv', sep= '\t', header = T)
d4 = read.table(file = 'clinical/TCGA_clinique_45_full.tsv', sep= '\t', header = T)
tmp = cbind(d4[,c('age','BCR_sensor','gleason_score','pathologic_t','psa_preop','BCR_sensor_months')])
d6 = cbind(tmp,d3)

write.table(d6, file = 'data/TCGA_cox', row.names = F, sep = '\t' )

d_mi = read.table(file = 'results_files/TCGA_miRNA_quant_cpm_52_threshold10.tsv', sep= '\t', header = T)
row.names(d_mi) = d_mi$barcode
row.names(d4) = d4$barcode

d_mi = d_mi[intersect(row.names(d_mi),row.names(d4)),]

write.table(d_mi, file = 'results_files/TCGA_miRNA_quant_cpm_45_threshold10.tsv', row.names = F, sep = '\t' )
