.lass <-
function(z, lambda){
  thetahat<-sign(z)*max((abs(z)-lambda),0)
  return(thetahat)
}
.plot.bic <- function (x,y,numid){
  x = x[numid]
  y = y[numid,]
  ymax = max(y)
  ymin = min(y)
  a=plot(x,y[,1],type="l",xlab="Lambda",ylab="pseudo-BIC",
         main = "Model Selection",ylim = c(ymin,ymax))
  points(x,y[,2],col="blue",type = "l", lty = 2)
  points(x,y[,3],col="red",type = "l", lty = 2)
  legend("bottomleft", c("pseudo-BIC","-2Loglikelihood","Estimated Degree of Freedom"),
         cex = 0.5,lty =c(1,2,2),col = c("black","blue","red"),
         pch=c(NA,NA,NA))
  return(a)
}
.scad <-
function(z,lambda){
  gamma<-3.7
  if (abs(z)>(lambda*gamma)){
    thetahat<-z
  } else if (abs(z)<=(2*lambda)){
    thetahat<-.lass(z,lambda)
  } else{
    thetahat<-.lass(z,lambda*gamma/(gamma-1))/(1-1/(gamma-1))
  }
  return(thetahat)
}
