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

d = read.table(file = 'distributions_disease2.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#d = d[sample(nrow(d), 1000), ]

i = 1
for (row in 2:nrow(d)){
  cluster = strsplit(d[row,2], ',')[[1]]
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
  cluster = strsplit(d[row,2], ',')[[1]]
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

a = as.vector(table(unlist(cl1))/length(unlist(cl1)))
print(as.numeric(table(unlist(cl6))/length(unlist(cl6))))
print(table(unlist(cl1))/length(unlist(cl1)))
#ss
#a = 0.0398873768 0.0220553731 0.0168934772 0.0755513843 0.4026278742 0.0079774754 0.0337869545 0.0286250587 0.0483341154 0.0314406382 0.0004692633 0.0187705303 0.0009385265 0.0408259033 0.0215861098 0.0811825434, 0.0431722196 0.0070389489 0.0009385265 0.0201783200 0.0028155795 0.0192397935 0.0178320038 0.0178320038
#b = c(0.018382353, 0, 0.040441176, 0.352941176, 0.007352941, 0, 0.033088235, 0.033088235, 0.066176471, 0.088235294, 0, 0, 0, 0.077205882, 0.040441176, 0.018382353, 0.084558824, 0, 0, 0, 0, 0.014705882, 0, 0.125000000)
#c = c(0, 0.02857143, 0.22857143, 0, 0.02857143, 0, 0.20000000, 0, 0.25714286, 0, 0, 0, 0, 0.14285714, 0, 0.11428571, 0, 0, 0, 0, 0, 0, 0, 0)
#d = c(0.02941176, 0.02941176, 0.11764706, 0.11764706, 0, 0, 0.05882353, 0.05882353, 0.20588235, 0.11764706, 0, 0, 0, 0, 0.05882353, 0, 0.20588235, 0, 0, 0, 0, 0, 0, 0)
#e = c(0.033057851, 0.033057851, 0.008264463, 0.049586777, 0.090909091, 0, 0.016528926, 0, 0.082644628, 0.404958678, 0, 0, 0, 0.099173554, 0.024793388, 0.074380165, 0, 0, 0, 0.033057851, 0.008264463, 0.041322314, 0, 0)

#species
#a = c(0, 0.013, 0.039, 0.041, 0.004, 0, 0.015, 0.025, 0.044, 0.008, 0.432, 0, 0.240, 0.006, 0.034, 0.034, 0.025, 0.013, 0.021)
#b = c(0, 0, 0, 0, 0, 0, 0.031, 0, 0, 0, 0.058, 0, 0.75, 0.031, 0.022, 0.084, 0.008, 0.008, 0.004)
#c = c(0, 0, 0.032, 0.064516129, 0, 0, 0, 0.032258065, 0.361290323, 0, 0.135483871, 0, 0.006451613, 0, 0.012903226, 0.316129032, 0.012903226, 0.006451613, 0.019354839)
#d = c(0.205655527, 0, 0, 0, 0.015424165, 0, 0, 0, 0.023136247, 0.205655527, 0.205655527, 0.318766067, 0.007712082, 0, 0.017994859, 0, 0, 0)
#e = c(0, 0.01470588, 0.11764706, 0.45588235, 0, 0.04411765, 0.10294118, 0.05882353, 0, 0, 0, 0, 0.08823529, 0, 0.01470588, 0.02941176, 0, 0.01470588, 0.05882353)

plot(a, xlab = "", ylab="")
points(b, col='red', add=TRUE, xlab = "", ylab="")
points(c, col='blue', add=TRUE, xlab = "", ylab="")
points(d, col='green', add=TRUE, xlab = "", ylab="")
points(e, col='orange', add=TRUE, xlab = "", ylab="")


