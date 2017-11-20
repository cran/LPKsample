LP.smooth <-
function(CR,n,method){ ###--Method = "AIC" or "BIC"
  CR.s <- sort(CR^2,decreasing=TRUE,index=TRUE)$x
  aa <- rep(0,length(CR.s))
  if(method=="AIC"){ penalty <- 2}
  if(method=="BIC"){ penalty <- log(n)}
  aa[1] <- CR.s[1] - penalty/n
  if(length(aa)==1){
    CR<-CR*(aa[1]>0)
  }else{
    if(aa[1]< 0){ return(rep(0,length(CR))) }
    for(i in 2: length(CR.s)){
      aa[i] <- aa[(i-1)] + (CR.s[i] - penalty/n)
    }
  #plot(aa,type="b",ylab=method,cex.axis=1.2,cex.lab=1.2)
    CR[CR^2<CR.s[which(aa==max(aa))]] <- 0
  }
  return(CR)
}
