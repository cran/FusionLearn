\name{fusionbase}
\alias{fusionbase}
\alias{fusionbase.fit}
\title{
Fusion learning method for continuous responses
}
\description{
\code{fusionbase} conducts the group penalization to multiple linear models with a specified penalty value. \code{fusionbase.fit} can be used to search the best candidate model based on the pseudo Bayesian information criterion with a sequence of penalty values. 
}
\usage{
fusionbase(x, y, lambda, N, p, m, beta=0.1, thresh=0.05, 
           maxiter=30, methods="scad",Complete=TRUE)

fusionbase.fit(x, y, lambda, N, p, m, beta=0.1, thresh=0.05, 
               maxiter=30, methods="scad", Complete=TRUE, depen ="IND", a=1)
}

\arguments{
  \item{x}{
List. Listing matrices of the predictors from different platforms. 
}
  \item{y}{
List. A list of continuous responses vectors from different platforms following the same order as in \code{x}.
}
  \item{lambda}{
Numeric or vector. For \code{fusionbase}, lambda is a numeric value for the penalty; for \code{fusionbase.fit}, lambda is a vector with a list of penalty values. 
}
  \item{N}{
Numeric or vector. If only one numeric value is provided, equal sample size will be assumed for each data set. If a vector is provided, then the elements are the sample sizes for all the platforms.
}
  \item{p}{
Numeric. The number of predictors.
}
  \item{m}{
Numeric. The number of platforms.
}
  \item{beta}{
Numeric or Matrix. An initial value for the estimated parameters with dimensions nvars x nplatforms. The defaul value is 0.1.
}
  \item{thresh}{
Numeric. The stopping criteria. The default value is 0.05.  
}
  \item{maxiter}{
Numeric. Maximum number of iterations. The default value is 30.
}
  \item{methods}{
Character ("lass" or "scad"). \code{lass}: LASSO; \code{scad}: SCAD.
}
  \item{Complete}{
Logic input. If \code{Complete == TRUE}, the predictors \eqn{M_1},...,\eqn{M_p} are measured in all platforms. If \code{Compelte == FALSE}, in some platforms, not all of the predictors \eqn{\{M_1,M_2,...,M_p\}} are measured. The values of the corresponding estimated coefficients for the missing predictors will be \code{NA}.
}  
  \item{depen}{
Character. Input only for function \code{fusionbase.fit}. "IND" means the observations across different platforms are independent; "CORR" means the observations are correlated, and the sample sizes should be equal for different platforms. 
  }
  \item{a}{
Numeric. Input only for function \code{fusionbase.fit}. The free multiplicative constant used in \eqn{\gamma_n}. The default value is 1.
}
}

\details{
The basic fusion learning function to learn from multiple linear models with continuous responses. More details regarding the model assumptions and the algorithm can be found in \code{\link{FusionLearn}}.
}
\value{

\code{fusionbase} returns a list that has components:

  \item{beta}{
A matrix (nvars x nplatforms) containing estimated coefficients of each linear model. If some data sets do not have the complete set of predictors, the corresponding coefficients are output as \code{NA}.
}

  \item{method}{
Penalty function LASSO or SCAD.
}

  \item{threshold}{
The numeric value shows the difference in the estimates between the successive updates upon convergence. 
}

  \item{iteration}{
The numeric value shows the number of iterations upon convergence.
}


\code{fusionbase.fit} provides the results in a table:

  \item{lambda}{
The sequence of penalty values.
}

  \item{BIC}{
The pseudolikelihood Bayesian information criterion evaluated at the sequence of the penalty values.
}

  \item{-2Loglkh}{
Minus twice the pseudo loglikelihood of the chosen model.
}

  \item{Est_df}{
The estimated degrees of freedom quantifying the model complexity.
}

\code{fusionbase.fit} also returns a model selection plot showing the results above. 
}
\note{
The range of the penalty values should be carefully chosen. For some penalty values, the resulting models may have singular information matrix or the fitting of the glm cannot converge.
}
\references{
Gao, X and Carroll, R. J. (2017) Data integration with high dimensionality. Biometrika, 104, 2, pp. 251-272
}
\author{
Xin Gao, Yuan Zhong, and Raymond J. Carroll
}

\examples{
##analysis of the stock index data
#Responses contain indices "VIX","GSPC", and "DJI" 
y <- list(stockindexVIX[,1],stockindexGSPC[,1],stockindexDJI[,1]) 

#Predictors include 46 stocks
x <- list(stockindexVIX[,2:47],stockindexGSPC[,2:47],stockindexDJI[,2:47])  

##Implementing the model selection algorithm based on the psuedolikelihood 
##information criteria  
model <- fusionbase.fit(x,y,seq(0.03,5,length.out = 10),232,46,3,depen="CORR")
lambda <- model[which.min(model[,2]),1]
result <- fusionbase(x,y,lambda,232,46,3)

##Identify the significant predictors for the three indices
id <- which(result$beta[,1]!=0)+1
colnames(stockindexVIX)[id]
}




