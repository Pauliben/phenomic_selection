snp_imputation=function(X,maf=0.05,ploidy=2,weighting=FALSE,ctr=FALSE,std=FALSE){
  
  #Note: X is the snp marker matrix
  S <- 0
  prop.MAF.j <- maf
  
  #Naive imputation 
  for(i in 1:ncol(X))
  {
    #meanXi <- mean(X[,i],na.rm=TRUE)
    meanXi <- round(mean(X[,i],na.rm=TRUE),0)
    X[,i] <- ifelse(is.na(X[,i]),meanXi, X[,i]) 
    S<-S+var(X[,i])
  }
  
  p <- colMeans(X,na.rm=T)/ploidy    
  p <- ifelse(p <= 0.5,p,1-p)
  index <- which( p >= prop.MAF.j ) 
  X <- X[,index]
  
  #Centering, standarizing and weighting.
  for(i in 1:ncol(X))
  {
    meanXi <- mean(X[,i],na.rm=TRUE)
    if(ctr){ X[,i]<-X[,i]-meanXi } #centering
    if(std){ X[,i]<-X[,i]/sd(X[,i]) } # standarizing
    if(weighting){ X[,i]<-X[,i]*weight[i] } # weighting
  }
  
  #Remove empty snps/redundant
  X2 = X[ , colSums(is.na(X))==0]
  
  return(X2)
}