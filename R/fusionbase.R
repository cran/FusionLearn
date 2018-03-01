fusionbase<-function(x,y,lambda,N,p,m,beta=0.1,thresh=0.05,maxiter=30,
                    methods="scad",Complete=TRUE){  
  diff<-10
  niter<-0
  ##Sample size
  if (length(N)==1) {
    N <- rep(N, m)
  } 
  M <- max (N)
  predictors <- array(dim=c(M,p,m))
  for (k in 1:m){
    predictors[1:N[k],,k] <- matrix ( as.numeric (unlist(x[k])), ncol = p,nrow = N[k] )
  }
  x <- array(0,dim=c(M,m,p))
  for (j in 1:p){
    for ( k in 1:m){
      x[1:N[k],k,j] <- predictors[1:N[k],j,k]
    }
  }
  response <- array(0, dim = c(M,m))
  for(k in 1:m) {
    response[1:N[k],k] <- as.numeric(unlist(y[k]))
  }
  y <- response
  
  ## initial beta
  if (length(beta)==1){
    beta = array(beta,dim = c(
      dim(x)[3],dim(x)[2]))
  }
  ## check x and y's dimensions
  if(dim(y)[1]!=dim(x)[1]){
    stop("dimensions of X and Y are not aligned")
  }
  if(dim(y)[2]!=dim(x)[2]){
    stop("dimensions of X and Y are not aligned")
  } 
  ## check input values
  if(max(N)!=dim(x)[1]){
    stop("dimension of X is not satisfied ")
  }
  if(p!=dim(x)[3]){
    stop("dimension of X is not satisfied ")
  }
  if(m!=dim(x)[2]){
    stop("dimension of X is not satisfied ")
  }
  
  ##Complete predictor
  K =matrix(0,p,m)
  if (Complete == FALSE) {
    x[is.na(x)]<-0
    for ( j in 1:p){
      for ( s in 1:m){
        K[j,s] = as.numeric(sum(abs(x[,s,j]))==0,1) 
      }
    }
  }
  k = rowSums(K)
  k = m/(m-k)
  
  ## Standardize x
  for (j in 1:p){
    for (s in 1:m){
      x[,s,j]<-x[,s,j]/sqrt(sum(x[,s,j]^2))*sqrt(N[s])
    }
  }
  x[is.nan(x)] <-0
  ## Main
  resi<-array(dim=c(max(N),m))
  z<-array(dim=c(p,m))
  while ((diff>thresh)&&(niter<maxiter)){
    for (s in 1:m){
      resi[,s]<-y[,s]-x[,s,]%*%beta[,s]
    }
    oldbeta<-beta
    for (j in 1:p){
      for (s in 1:m){
        z[j,s]<-(x[,s,j]%*%resi[,s])/N[s]+beta[j,s]
      }
      z[which(K==1)] = 0 
      norm<-sqrt(sum(z[j,]^2)*k[j])
      if ("lass" %in% methods){
        beta[j,]<- .lass(norm,lambda)*(z[j,])/norm
      }
      if ("scad" %in% methods) {
        beta[j,]<- .scad(norm,lambda)*(z[j,])/norm
      }
      for (s in 1:m){
        resi[,s]<-resi[,s]-x[,s,j]*(beta[j,s]-oldbeta[j,s])
      }
    }
    diff<-max(abs(beta-oldbeta))
    niter<-niter+1
  }
  if (Complete==FALSE){
    beta[ which(K==1)] = NA
  }
  mylist <- list("beta"=beta, "method"=methods,
                 "threshold"=diff,"iteration"=niter)
  return(mylist)
}
