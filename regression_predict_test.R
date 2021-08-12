d = read.table(file = '/Users/rerleman/Dropbox/My Mac (CCI00BHV2JALT)/Documents/distributions2.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

obsexp = function(cluster,db){
  all = unique(c(cluster,db))
  obs = c()
  exp = c()
  for (el in all){
    obs = c(obs,length(cluster[cluster==el]))
    exp = c(exp,length(db[db==el]))
  }
  return(list(obs,exp))
}

for (row in 2:nrow(d)){
  cluster = strsplit(d[row,11], ',')[[1]]
  db = strsplit(d[row,12], ',')[[1]]
  conttable = obsexp(cluster,db)
  obs = conttable[[1]]
  exp = conttable[[2]]/length(db)
  pvalue = chisq.test(obs,p=exp,simulate.p.value = TRUE)$p.value
  d[row,13] = pvalue
}

d = d[-c(11,12)]
d[1,11] = 'Chi square bootstrap p-value'
d = d[-10]
write.table(d, file = '/Users/rerleman/Dropbox/My Mac (CCI00BHV2JALT)/Documents/distributions2.1.tsv', sep = '\t')


plot(d[,6],d[,10], xlab='Asymptotic Chi-square', ylab='Bootstrap Chi-square')
