fusionmixed.fit<-function(x,y,lambda,N,p,m1,m2,beta=0.1,thresh=0.1,maxiter=100, 
                          methods="scad",link="logit", Complete=TRUE,depen ="IND",a=1){
  m = m1+m2
  ##Sample size
  if (length(N)==1) {
    N <- rep(N, m)
  } 
  M <- max (N)
  predictor <- array(dim=c(M,p,m))
  for (k in 1:m){
    predictor[1:N[k],,k] <- matrix ( as.numeric (unlist(x[k])), ncol = p,nrow = N[k] )
  }
  predictors <- array(0,dim=c(M,m,p))
  for (j in 1:p){
    for ( k in 1:m){
      predictors[1:N[k],k,j] <- predictor[1:N[k],j,k]
    }
  }
  response <- array(0, dim = c(M,m))
  for(k in 1:m) {
    response[1:N[k],k] <- as.numeric(unlist(y[k]))
  }
  
  n=length(lambda)
  output<-matrix(0,n,3)
  commen <- as.character(rep(".",n))
  
  for ( k in 1:n){
    
    result = fusionmixed(x,y,lambda[k],N,p,m1,m2,beta,thresh,maxiter,methods,
                         link,Complete)
    result = result$beta
    if ( Complete == FALSE){
      result[is.na(result)] <- 0
    }
    index <- result[,1]!=0
    p1 <- sum (index)
    if (p1==0){
      output[k,] = NA
      commen[k] = "p = 0, tuning parameter is too large"
      next
    }
    if (p1 >= min(N)){
      output[k,] = NA
      commen[k] = "p >= n, glm cannot converge at this time"
      next
    }
    xnam <- paste("X", (1:p)[index], sep="")
    fmla.s <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
    
    informat <- array(0,dim=c(p1*m,p1*m))
    score <- array(0, dim=c(max(N),p1*m))
    jest <- array(0,dim = c(p1*m,p1*m))
    m2log <- numeric (m)
    for (s in 1:m1){
      data = data.frame(cbind(response[1:N[s],s],predictors[1:N[s],s,]))
      data[is.na(data)]<-0
      colnames(data)[1]<-"y"
      colnames(data)[2:(p+1)]<-paste("X", 1:p, sep="")
      submodel <- glm(fmla.s,data,family=gaussian(link=identity))
      m2log[s] <- submodel$aic - 2*p1 #
      
      wmat <- diag(submodel$weights)
      designmat <- as.matrix(data[,2:(p+1)][index])
      informat[((s-1)*p1+1):(s*p1),((s-1)*p1+1):(s*p1)] <- t(designmat)%*%wmat%*%designmat #
      sv <- data$y - submodel$fitted.values
      score[1:N[s],((s-1)*p1+1):(s*p1)] <- designmat*c(sv) #
    }
    for (s in (m1+1):m){
      data = data.frame(cbind(response[1:N[s],s],predictors[1:N[s],s,]))
      data[is.na(data)]<-0
      colnames(data)[1]<-"y"
      colnames(data)[2:(p+1)]<-paste("X", 1:p, sep="")
      submodel <- glm(fmla.s,data,family=binomial(link=link))
      m2log[s] <- submodel$aic - 2*p1 #
      
      wmat <- diag(submodel$weights)
      designmat <- as.matrix(data[,2:(p+1)][index])
      informat[((s-1)*p1+1):(s*p1),((s-1)*p1+1):(s*p1)] <- t(designmat)%*%wmat%*%designmat #
      if ("logit" %in% link){
        sv <- data$y - submodel$fitted.values
      }
      if ("probit" %in% link){
        linearp<-designmat%*%c(submodel$coefficients[-1])
        predp<-pnorm(linearp,0,1)
        vec<-(abs(linearp)>3)*abs(linearp)/pnorm(abs(linearp),0,1)
        temp<-(abs(linearp)<=3)*dnorm(linearp,0,1)/(predp*(1-predp))
        temp[is.na(temp)]<-0
        vec<-vec+temp
        sv<-(data$y-predp)*vec
      }
      score[1:N[s],((s-1)*p1+1):(s*p1)]<-designmat*c(sv)
    }
    ## singularity 
    if(det(informat)==0){
      output[k,] = NA
      commen[k] = "information matrix computionally singular"
      next
    }
    if("IND" %in% depen){
      for(s in 1:m){
        jest[((s-1)*p1+1):(s*p1),((s-1)*p1+1):(s*p1)] = var(score[1:N[s],((s-1)*p1+1):(s*p1)] )*N[s]
      }
    }
    if("CORR" %in% depen){
      jest <- cov(score)*max(N)
    }
    output[k,2] = sum(m2log)
    output[k,3] = a*log(p)*sum(diag(jest%*%qr.solve(informat)))
  }
  numid <- which(!is.na(output[,3]))
  output[,1] = output[,2]+output[,3]
  if (length(lambda[numid])>1){
    bic.plot =.plot.bic(lambda,output,numid)
  }
  commen[which.min(output[,1])] = "The minimum BIC"
  output = data.frame(lambda,output,commen)
  colnames(output) <- c("lambda","BIC","-2Loglkh","Est_Df","Comment")
  return(output)
  return(bic.plot)
}


