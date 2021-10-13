library(extraDistr)
library(pracma)

#Calculate chi-square test statistic value
chisquare = function(cluster,db){
  all = unique(c(cluster,db))
  n = length(cluster)
  m = length(db)
  S = 0
  #Add 1/2 to both obs and exp to avoid division by 0
  for (el in all){
    obs = length(cluster[cluster==el]) + 1/2
    exp = length(db[db==el])*(n/m) + 1/2
    S = S + ((obs-exp)**2)/exp
  }
  return(S)
}

chi_pvalue = function(cluster,db,iter){
  n = length(cluster)
  m = length(db)
  statobs = chisquare(cluster,db)
  statdistr = c()
  for (kk in 1:iter){
    permutation = randperm(c(cluster,db), n+m)
    sample1 = permutation[1:n]
    sample2 = permutation[(n+1):(n+m)]
    statdistr = c(statdistr, chisquare(sample1,sample2))
  }
  pvalue = length(statdistr[statdistr>=statobs])/iter
  return(pvalue)
}

d = read.table(file = '/Users/rerleman/Documents/adap-kdp-gof-data/data1.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE, fill=FALSE, quote="")

for (row in 2:nrow(d)){
  cluster = strsplit(d[row,10], ',')[[1]]
  if (length(cluster)<10){
    d[row,12] = "2"
  }
  else{
    db = strsplit(d[row,11], ',')[[1]]
    pvalue = chi_pvalue(cluster,db,10000)
    d[row,12] = pvalue
  }
}

d[,10] = d[,12]
d = d[,-11]
d = d[,-11]

write.table(d, file = '/Users/rerleman/Documents/adap-kdp-gof-data/data2.tsv', sep = '\t')