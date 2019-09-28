# load data
setwd("E:/cloud/Code/R")
.libPaths("win-library")
library(ComplexHeatmap)
library(circlize)

# load data 
MyData_full = 
  read.table("E:/cloud/Data/TCGA_PRAD/datamining/TCGA_BCR_a.classification.data_to_train.csv",
             header=TRUE, sep="\t", na.strings="", dec=".", strip.white=TRUE, check.names=FALSE)
MyData = 
  read.table("E:/cloud/Data/TCGA_PRAD/datamining/TCGA_BCR_d.classification.model_features.csv",
             header=TRUE, sep=",", na.strings="", dec=".", strip.white=TRUE, check.names=FALSE)


#convert first column as row names
rownames(MyData)=MyData[,1] #set first column as row names
MyData=MyData[,-1] #remove first column that is useless now
#MyData=MyData[,1:ncol(MyData)-1] # remove last column (it's actually empty)

#genesOfInterest and measures we want in the heatmap
genesMat = as.matrix(MyData[,grep("\\|", colnames(MyData))])
genesMat=scale(genesMat)

#configure heatmap parts
haBottom = HeatmapAnnotation(boxplot = anno_boxplot(genesMat, axis = TRUE, axis_side = "left", size = unit(1, "mm"), pch = 16), 
                             name = "Gene expression (scaled)")


ha1 = rowAnnotation(df=data.frame(Class=MyData$gleason.CODE), 
                    annotation_legend_param = list(title = "Class"),width = unit(3, "mm"),
                    col = list(Class= c("g06" =  "limegreen", "g07_4_3" = "cyan", "g07_3_4" = "cyan",
                                             "g08" = "cyan4", "g09_10" = "blueviolet")))
ha2 = rowAnnotation(df=data.frame(pathological_t=MyData$pathologic_t), 
                    annotation_legend_param = list(title = "Pathological T"),width = unit(3, "mm"),
                    col = list(pathological_t = c("t2a" = "olivedrab1", "t2b" = "olivedrab3", "t2c" = "olivedrab", 
                                                  "t3a" = "orange1", "t3b" = "orange4", "t4" = "orangered", "?" = "white")))
ha3 = rowAnnotation(df=data.frame(error=MyData$`class`), 
                    annotation_legend_param = list(title = "BCR"),width = unit(3, "mm"),
                    col = list(error= c("FALSE" =  "grey", "TRUE" = "black")))

ht1 = Heatmap(genesMat, name = "Gene expression", bottom_annotation = haBottom, 
              top_annotation_height = unit(4, "mm"), 
              split= MyData$`class`,
              #km=4,
              bottom_annotation_height = unit (2, "cm"), show_row_names = F, show_column_names = TRUE, 
              column_names_side = "top",column_names_gp = gpar(fontsize = 4), cluster_columns = T) 

# draw heatmap
ht_list=ht1+ha1+ha2+ha3
draw(ht_list)

##PCA patients
library(FactoMineR)
MyData_genes=MyData[,grep("\\|", colnames(MyData))]
MyData_genes=cbind(MyData_genes)
pca = PCA(MyData_genes, scale.unit=TRUE, ncp=5, graph=T)
pca = PCA(MyData_genes, scale.unit=TRUE, ncp=5, quanti.sup=c(11: 12), quali.sup=13, graph=T)

plot.PCA(pca, axes=c(1, 2), choix="ind", label="none", habillage=MyData$class)



plotPCA(as.matrix(t(MyData_genes)), col=rep(1:2, each=2), cex=1.2, labels=F, col(t(MyData_genes$class)))  


MyData_genes=MyData_sign[,grep("\\|", colnames(MyData_sign))]
plotPCA(as.matrix(t(MyData_genes)), col=rep(1:2, each=2), cex=1.2, labels=F)  

