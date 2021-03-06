fusionmixed<-function(x,y,lambda,N,p,m1,m2,beta=0.1,thresh=0.1,maxiter=100,
                     methods="scad",link="logit", Complete=TRUE){
  
  diff<-10
  niter<-0
  m=m1+m2
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
  
  ## check x and y's dimensions
  if(dim(y)[1]!=dim(x)[1]){
    stop("dimensions of X and Y are not aligned")
  }
  if(dim(y)[2]!=dim(x)[2]){
    stop("dimensions of X and Y are not aligned")
  } 
  
  ##Missing Values
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
  
  ## main
  resi<-array(dim=c(max(N),m))
  z<-array(dim=c(p,m))
  while ((diff>thresh)&&(niter<maxiter)){
    for (s in 1:m1){
      resi[,s]<-y[,s]-x[,s,]%*%beta[,s]
    }
    for (s in (m1+1):m){
      if (length(unique(y[,s]))!=2){
        warning("Binary data is not detected")
      }
      if ( "probit" %in% link){
        vu <- 2/pi
        linearp<-x[,s,]%*%beta[,s]
        predp<-pnorm(linearp,0,1)
        vec<-(abs(linearp)>3)*abs(linearp)/pnorm(abs(linearp),0,1)
        temp<-(abs(linearp)<=3)*dnorm(linearp,0,1)/(predp*(1-predp))
        temp[is.na(temp)]<-0
        vec<-vec+temp
        resi[,s]<-(y[,s]-predp)/vu*vec
      }
      if ("logit" %in% link){
        vu=1/4
        linearp<-x[,s,]%*%beta[,s]
        predp<-1-1/(1+exp(linearp))
        resi[,s]<-(y[,s]-predp)/vu  
      }
    }
    oldbeta<-beta
    for (j in 1:p){
      for (s in 1:m1){
        z[j,s]<-(x[,s,j]%*%resi[,s])/N[s]+beta[j,s]
      }
      for (s in ((m1+1):m)){
        z[j,s]<-((x[,s,j]%*%resi[,s])/N[s]+beta[j,s])*vu
      }
      z[which(K==1)] = 0 
      norm<-sqrt(sum((z[j,])^2)*k[j])
      if ("lass" %in% methods ){
        beta[j,]<- .lass(norm,lambda)*(z[j,])/(norm)
      }
      if ("scad" %in% methods){
        beta[j,]<- .scad(norm,lambda)*(z[j,])/(norm)
      }
      beta[j,(m1+1):m]<- beta[j,(m1+1):m]/vu
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
  mylist <- list("beta"=beta, "method"=methods, "link_fn"=link,
                 "threshold"=diff,"iteration"=niter)
  return(mylist)
}
