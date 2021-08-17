library(extraDistr)
library(pracma)

#Significance level
sig = 0.05

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

#Label probabilities
clusterprobs = c(0, 0.01470588, 0.11764706, 0.45588235, 0, 0.04411765, 0.10294118, 0.05882353, 0, 0, 0.08823529, 0, 0.01470588, 0.02941176, 0, 0.01470588, 0.05882353, 0)
dbprobs = c(0, 0.012906448, 0.039190898, 0.041297935, 0.004214075, 0.015170670, 0.025705858, 0.044669195, 0.008849558, 0.432364096, 0.240202276, 0.006742520, 0.034134008, 0.034134008, 0.025284450, 0.013063633, 0.021070375, 0)

#Sample sizes
nvec = c(3,5,10,20,50)
m = 50
powvec = c()
for (n in nvec){

  #Number of Monte Carlo simulations for approximating the power
  iter = 1000
  #Number of Monte Carlo simulations for the permutation test
  iter2 = 1000
  
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
  powvec = c(powvec, sum(pow)/iter)
}