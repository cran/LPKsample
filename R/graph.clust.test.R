graph.clust.test <-
function(y,W,method='kmeans',LP.type='smooth',smooth.method='BIC',return.clust=F){
    out<-ls()
    k<-length(unique(y))
    WgtG<-igraph::graph.adjacency(W,weighted=TRUE,mode='undirected')
    L1<-diag(dim(W)[1])-as.matrix(igraph::graph.laplacian(WgtG,normalized=TRUE))
    Lap.svd<-svd(L1)

    U.Lap<-diag(sqrt(sum(W)/rowSums(W)))%*%(Lap.svd$u[,2:k])
    if(method=='mclust'){
      m1<-mclust::Mclust(U.Lap,G=k,verbose=F)
      y.c<-m1$classification
    }else if(method=='kmeans'){
      if(sum(1-is.na(U.Lap[cumsum(table(y)),]))==sum(1-is.na(unique(U.Lap[cumsum(table(y)),])))){
	   y.c<-kmeans(U.Lap,U.Lap[cumsum(table(y)),],iter.max=50)$cluster
      }else{
	   y.c<-kmeans(U.Lap,k,iter.max=20)$cluster
      }
    }
    test0<-LP.comean(y,y.c,type=LP.type,smooth.method)
    if(return.clust==F){
	   return(test0)
    }else if(return.clust==T){
	   out$test<-list()
	   out$clust<-y.c
	   out$test<-test0
	   return(out)
	}
	
}
