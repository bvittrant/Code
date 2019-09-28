###############################################################################

# test Cox

library("survival")
library("survminer")

###############################################################################

# Loading data 
df_TCGAGeneClin = read.table('data/TCGAcox.tsv', header = T, check.names = F)#, fill = T
row.names(df_TCGAGeneClin) = df_TCGAGeneClin$Row.names
df_TCGAGeneClin = df_TCGAGeneClin[,-1]


vec1 = c("BCR_sensor", "BCR_sensor_months","ENSG00000159189","ENSG00000274780", "ENSG00000074181", "ENSG00000151651", "ENSG00000106991")
d1 = df_TCGAGeneClin[,vec1]

###############################################################################

covar = c('ENSG00000159189','ENSG00000274780', 'ENSG00000074181', 'ENSG00000151651', 'ENSG00000106991')
df_logtest = c()
vec_name = c()

comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}


for(i in 1:length(covar)){
  for(j in 1:length(combn(covar,i, simplify = F))){
    
    # create p among n combination (tirage au sort sans ordre et sans remise)
    tmp = unlist(strsplit(unlist(combn(covar,i, simplify = F)[j]), "[ ] "))
    
    # Create formula combinaison of the covariates
    tmp = paste(tmp,collapse = '*')
    
    # create formula
    fl1= formula(paste('Surv(BCR_sensor_months, BCR_sensor)~',tmp,sep=''))
    
    # Analyse de survie
    res.cox <- coxph(fl1, data =  d1)
    sum.res.cox = summary(res.cox)
    
    # assign to vec_name the combinaisaons as futur row.name
    assign(paste(covar[i:j], collapse = '-'),c(sum.res.cox$logtest,sum.res.cox$waldtest,sum.res.cox$concordance,sum.res.cox$rsq,sum.res.cox$sctest))
    
    # concatenate several results of the cox analysis to the summarry DF
    df_logtest  = rbind(df_logtest,eval(as.name(paste(covar[i:j], collapse = '-'))))
    
    # prepare and assign vec name
    vec_name = cbind(vec_name,tmp)
    row.names(df_logtest) = vec_name
  }
}

# compute AIC 2k−2ln(L) k = df et L = likehood 
sum.res.cox$loglik

# Write results
write.table(df_logtest, file = 'data/df_summary.tsv', sep = '\t', row.names = T, quote = F)
