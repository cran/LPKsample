GLP <-
function(X,y,m.max=4,alpha=0.05,c=0.5,clust.alg='kmeans',LPtype='smooth',smooth.method='BIC',return.LPT=F,return.clust=F){

   t.res<-matrix(0,m.max,3)
   colnames(t.res)<-c('component','comp.GLP','pvalue')
   num.g<-length(unique(y))
   lp.sig<-c()
   LPTmat<-list()

   

   #check maximum order
   if(m.max<1 | m.max!=floor(m.max)){
      stop("maximum order m.max must be a positive integer")
   }else if(m.max>=1){
   #dual loop for centralize adjustment:
       W.overall<-matrix(0,nrow(X),nrow(X))
   #start with odd m:
       for(i in 1:ceiling(m.max/2)){
          m.i<-2*i-1
          w0<-W.Gen(X,m.i,c)
          res1<-graph.clust.test(y,w0$W,method=clust.alg,LP.type=LPtype,smooth.method)
          t.res[m.i,]<-c(m.i,res1$LPINFOR,res1$p.val)
          if(t.res[m.i,3]<=alpha){
              W.overall<-W.overall+w0$W;
              LPTmat[[m.i]]<-w0$LPT
          }
       }
   #check if median adjustment needed
       if(t.res[1,3]<alpha){X<-apply(X,2, function(x,Y){ x - rep(aggregate(x,list(Y),median)$x, 
          as.vector(table(Y)))}, Y=y)}
   #continue with even m:
       if(m.max>1){
          for(i in 1:floor(m.max/2)){
             m.i<-2*i
             w0<-W.Gen(X,m.i,c)
             res1<-graph.clust.test(y,w0$W,method=clust.alg,LP.type=LPtype,smooth.method)
             t.res[m.i,]<-c(m.i,res1$LPINFOR,res1$p.val)
             if(t.res[m.i,3]<=alpha){
                  W.overall<-W.overall+w0$W;
                  LPTmat[[m.i]]<-w0$LPT
             }
          }
       }
   }

   #record significant components
   lp.sig<-as.numeric(which(t.res[,3]<=alpha))
   #graph-based test:
   if(length(lp.sig)>0){  #when there are significant components
       res0<-graph.clust.test(y,W.overall,method=clust.alg,LP.type=LPtype,smooth.method,return.clust)
   }else{                 # otherwise, p-value forced to 1
       res0<-list()
       res0$LPINFOR<-0
       res0$p.val<-1
       res0$test<-list()
       res0$test$LPINFOR<-0
       res0$test$p.val<-1
       res0$clust<-NULL
   }

   #output
   out<-list()
   if(return.clust==F){
       result<-res0
   }else if(return.clust==T){
       result<-res0$test
   }

   
   out$GLP<-result$LPINFOR
   out$pval<-result$p.val
   out$table<-t.res
   out$component<-lp.sig
   if(return.LPT==T){
      out$LPT<-do.call(cbind,LPTmat)
   }
  
  
   if(return.clust==T){
      out$clust<-res0$clust
   }
  
  return(out)

}
