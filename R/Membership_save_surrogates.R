# install.packages("R.matlab")
# install.packages("igraph")
# library(R.matlab)
# library(igraph)
# Set path below
setwd("/Volumes/My Passport for Mac/code/substorm_community/Matlab_Functions/Communities_paper/R")
for(surr in seq(1,10,by=1)){
Substorms=readMat(paste("RSubstorms_surrogate_w128_t5_",surr,".mat",sep=""))
mysubs<-0
for(n in c(1,13,18,29,33)){
  mysubs<-Substorms[["Substorms.surrogate.w128.t5"]][[n]][[1]]
  m=dim(mysubs)[3]
  nstats=dim(mysubs)[2]
  
  modularity_eb<-matrix(nrow=m,ncol=1)
  cluster_list_eb<-matrix(nrow=m,ncol=nstats)
  
  for(i in seq(1,m,by=1)){
    g=graph_from_adjacency_matrix(mysubs[,,i],weighted=NULL,diag=FALSE)
    eb<-cluster_edge_betweenness(g,weights=NULL,directed=TRUE)
    
    count<-1
    q<-membership(eb)
    cluster_list_eb[i,]<-c(q)
    modularity_eb[i,]<-modularity(eb)
    
    print(i)
  }
writeMat(paste(n,"modularity_eb_surrogate_w128_t5_",surr,".mat",sep=""), modularity_eb=modularity_eb,fixNames = TRUE,onWrite=NULL)
writeMat(paste(n,"communities_eb_surrogate_w128_t5_",surr,".mat",sep=""), cluster_list_eb=cluster_list_eb,fixNames = TRUE,onWrite=NULL)}}