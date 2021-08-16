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
plot(clusters, ylab = "JS distance", xlab = "Cluster samples species")

mycl <- cutree(clusters, h=0.4)

cl1 = list()
cl2 = list()
cl3 = list()
cl4 = list()
cl5 = list()
cl6 = list()


i=1
for (row in 2:nrow(d)){
  cluster = strsplit(d[row,11], ',')[[1]]
  if (length(cluster)>2){
    if (mycl[i]==1){
      cl1[[length(cl1)+1]] = cluster
    }
    if (mycl[i]==2){
      cl2[[length(cl2)+1]] = cluster
    }
    if (mycl[i]==3){
      cl3[[length(cl3)+1]] = cluster
    }
    if (mycl[i]==4){
      cl4[[length(cl4)+1]] = cluster
    }
    if (mycl[i]==5){
      cl5[[length(cl5)+1]] = cluster
    }
    if (mycl[i]==6){
      cl6[[length(cl6)+1]] = cluster
    }
    i = i+1
  }
}

barplot(table(unlist(cl5)))
