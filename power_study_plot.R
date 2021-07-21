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
    obs1 = length(cluster[cluster==el]) + 1/2
    obs2 = length(cluster[cluster==el]) + 1/2
    exp1 = n*(length(cluster[cluster==el])+length(db[db==el]))/(n+m)
    exp2 = m*(length(cluster[cluster==el])+length(db[db==el]))/(n+m)
    S = S + ((obs1-exp1)**2)/exp1 + ((obs2-exp2)**2)/exp2
  }
  return(S)
}


#Label probabilities
clusterprobs = rep(1/10,10)
dbprobs = rep(1/10,10)

#Sample sizes
nvec = c(3, 5, 10, 20, 50)
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

plot(x=nvec, y=powvec, ylim=c(0,1), xlab="", ylab='')
lines(x=nvec, y=powvec)
