#Calculate Jensen-Shannon divergence between two samples
JSdiv = function(cluster,db){
  all = unique(c(cluster,db))
  n = length(cluster)
  m = length(db)
  S = 0
  for (el in all){
    P = length(cluster[cluster==el])/n
    Q = length(db[db==el])/m
    M = 1/2*(P+Q)
    if (P==0){
      DPQ = 0
      DQP = Q*log(Q/M)
    }
    if (Q==0){
      DPQ = P*log(P/M)
      DQP = 0
    }
    if (P!=0 & Q!=0){
      DPQ = P*log(P/M)
      DQP = Q*log(Q/M)
    }
    S = S + 1/2*DPQ + 1/2*DQP
  }
  return(S)
}

JSmat = function(x){
  n = length(x)
  Dmat = matrix(nrow=n,ncol=n)
  for (ind1 in 1:n){
    for (ind2 in 1:n){
      sample1 = x[[ind1]]
      sample2 = x[[ind2]]
      Dmat[ind1,ind2]=JSdiv(sample1,sample2)
    }
  }
  return(Dmat)
}

d1 = c(1,2,3)
d2 = c(2,2,1)
d3 = c(1,2,4,4)
d4 = c(1,1,2,3)

testmat = as.dist(JSmat(list(d1,d2,d3,d4)))
clusters = hclust(testmat)
plot(clusters, ylab = "JS distance", xlab = "Samples")


