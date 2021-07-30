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

d1 = c(1,2,3)
d2 = c(2,2,1)
d3 = c(1,2,4,4)

cat("JS-divergence between d1 and d2 is ", JSdiv(d1,d2))
cat("JS-divergence between d1 and d2 is ", JSdiv(d1,d3))
cat("JS-divergence between d1 and d2 is ", JSdiv(d2,d3))


