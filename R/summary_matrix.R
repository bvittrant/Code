d_tmp= matrix(ncol=6, nrow = 107)
colnames(d_tmp) = c('Min','1stQ','Median','Mean','3rdQ','Max')
for(i in 2:dim(d1)[2]){
  tmp = summary(d1[,i])
  d_tmp[i,1] = tmp[1]
  d_tmp[i,2] = tmp[2]
  d_tmp[i,3] = tmp[3]
  d_tmp[i,4] = tmp[4]
  d_tmp[i,5] = tmp[5]
  d_tmp[i,6] = tmp[6]
}
d1_2= as.data.frame(d_tmp[-1,])
