# explo CNV #

d_cnv = read.table(file = 'results_files/CNV_full_merged_clean_-1_1.tsv', sep = '\t', header = T, check.names = F)

# Remove duplicated rows (keep one of it !)
d_cnv = d_cnv[order(d_cnv$condition),]
d_cnv = d_cnv[!duplicated(d_cnv$condition),]

# Attribute row name
row.names(d_cnv) = d_cnv$condition
d_cnv = d_cnv[,-c(1:4)]


# Remove rows with only 1 values different of 0
vec_tmp = c()
for(i in 1:dim(d_cnv)[1]){
  count = 0
  for(j in 1:dim(d_cnv)[2]){if(d_cnv[i,j] == 0){count = count + 1}}
  
  if(count >= 497){vec_tmp = cbind(i,vec_tmp);print(i)}
}
d_cnv = d_cnv[-vec_tmp,]

#write.table(d_cnv, file = 'results_files/CNV_full_merged_clean_-1_1_row497zero.tsv', sep = '\t', row.names = T)

# Select the 52 patients

d_clin = read.table(file = 'clinical/TCGA_clinique_45_reduced.tsv', sep = '\t', header = T, check.names = F)
vec_name = as.character(d_clin$barcode)
d_cnv_clin = d_cnv[vec_name]
d_cnv_clin =d_cnv_clin[,order(colnames(d_cnv_clin))]

d_cnvt = as.data.frame(t(d_cnv_clin))

intersect(colnames(d_cnv), d_clin$barcode)

write.table(d_cnvt, file = 'results_files/CNV_45_merged_clean_-1_1_row497zero.tsv', sep = '\t', row.names = T)
