library(extraDistr)
library(pracma)

#Significance level
sig = 0.05

#Label probabilities
clusterprobs = c(1/20, 1/20 , 1/20, 1/20, 1/10, 1/10, 1/10, 1/10, 1/5, 1/5)
dbprobs = c(1/5, 1/5, 1/10, 1/10, 1/10, 1/10, 1/20, 1/20, 1/20, 1/20)

#Sample sizes
n = 15
m = 20

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

#Number of Monte Carlo simulations for approximating the power
iter = 1000
#Number of Monte Carlo simulations for the permutation test
iter2 = 500

pow = c()
for (k in 1:iter){
  cluster = rcat(n, clusterprobs)
  db = rcat(m, dbprobs)
  statobs = chisquare(cluster,db)
  
  statdistr = c()
  for (kk in 1:iter2){
    permutation = randperm(c(cluster,db), n+m)
    sample1 = permutation[1:n]
    sample2 = permutation[(n+1):(n+m)]
    statdistr = c(statdistr, chisquare(sample1,sample2))
  }
  pvalue = length(statdistr[statdistr>=statobs])/iter2
  if (pvalue <= sig){
    pow = c(pow, 1)
  }
}

cat("Power:", sum(pow)/iter)
