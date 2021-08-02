library(ape)

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

data = list()

#Set the working dir
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

d = read.table(file = 'distributions2species.txt', sep = '\t', header = FALSE, nrow=200, stringsAsFactors = FALSE)

i = 1
for (row in 2:nrow(d)){
  cluster = strsplit(d[row,11], ',')[[1]]
  data[[i]] = cluster
  i = i+1
}

testmat = as.dist(JSmat(data))
clusters = hclust(testmat)
#Dendrogram
#plot(clusters, ylab = "JS distance", xlab = "Cluster samples species")
res = pcoa(testmat)
plot(res$vectors[,1], res$vectors[,2],pch=20, col=cutree(clusters,3), xlab="1. principal coordinate", ylab="2. principal coordinate")
