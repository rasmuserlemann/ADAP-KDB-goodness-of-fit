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

d = read.table(file = 'distributions2species.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#d = d[sample(nrow(d), 1000), ]

i = 1
for (row in 2:nrow(d)){
  cluster = strsplit(d[row,11], ',')[[1]]
  if (length(cluster)>2){
    data[[i]] = cluster
    i = i+1
  }
}

testmat = as.dist(JSmat(data))
clusters = hclust(testmat, method = "average")
#Dendrogram
#plot(clusters, ylab = "JS distance", xlab = "Cluster samples species")
#res = pcoa(testmat)

mycl <- cutree(clusters, h=0.5)

#Variance explained by the PCoA
allvar = var(res$values[,1])+var(res$values[,2])+var(res$values[,3])+var(res$values[,4])+var(res$values[,5])+var(res$values[,6])
PC1var = round(var(res$values[,1])/allvar,2)
PC2var = round(var(res$values[,2])/allvar,2)
print(PC2var)

#plot(res$vectors[,1], res$vectors[,2], pch=20, col=cutree(clusters,3), xlab="1. PCo", ylab="2. PCo")
#plot(res$vectors[,1], res$vectors[,3],pch=20, col=cutree(clusters,3), xlab="1. principal coordinate", ylab="3. principal coordinate")
#plot(res$vectors[,2], res$vectors[,3],pch=20, col=cutree(clusters,3), xlab="2. principal coordinate", ylab="3. principal coordinate")
#clusterGroups<- cutree(clusters,k=3)
#cluster1ind = sample(which(clusterGroups==1),1)
#cluster2ind = sample(which(clusterGroups==2),1)
#cluster3ind = sample(which(clusterGroups==3),1)
#print(data[cluster1ind])
#print(data[cluster2ind])
#print(data[cluster3ind])


