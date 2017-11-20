LP.comean <-
function(x,y,type='smooth',smooth.method='BIC'){
  n<-length(x)
  if(length(unique(y))<=1 |length(unique(x))<=1  ){
    LP.mat<-0
  }else{
    m=length(unique(x))-1
    out<-list()

    Tx<-LP.Poly(x,m)
    Ty<-LP.Poly(y,m)

    LP.mat<-cor(Tx,Ty)
  }
  matdim<-dim(LP.mat)

  df0<-nrow(LP.mat)*ncol(LP.mat)

  if(type=='smooth'){
     LP.mat<-LP.smooth(c(LP.mat),n,method=smooth.method) 
     df0<-sum(LP.mat!=0) 
  }
  
   pval<-pchisq(n*sum(LP.mat^2),df=df0,lower.tail=FALSE)

   out$LPINFOR<-sum(LP.mat^2)
   out$p.val<-pval
   out$LP.matrix<-matrix(LP.mat,matdim[1],matdim[2])
   return(out)
}
