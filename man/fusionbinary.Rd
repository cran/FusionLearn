\name{fusionbinary}
\alias{fusionbinary}
\alias{fusionbinary.fit}
\title{
Fusion learning algorithm for binary responses
}
\description{
\code{fusionbinary} conducts the group penalization with a specified penalty value learning from multiple generalized linear models with binary responses. \code{fusionbinary.fit} can be used to search the best candidate model based on the pseudo Bayesian information criterion with a sequence of penalty values.
}
\usage{
fusionbinary(x, y, lambda, N, p, m, beta=0.1, thresh=0.1, 
             maxiter=100, methods="scad", link="logit", Complete=TRUE)

fusionbinary.fit(x, y, lambda, N, p, m, beta=0.1, thresh=0.1, 
                 maxiter=100, methods="scad", link="logit", Complete=TRUE, 
                 depen ="IND", a=1)
}
\arguments{
  \item{x}{
List. Listing matrices of the predictors from different platforms. 
}
  \item{y}{
List. A list of binary responses vectors from different platforms following the same order as in \code{x}.
}
  \item{lambda}{
Numeric or vector. For \code{fusionbinary}, lambda is a numeric value for the penalty; for \code{fusionbinary.fit}, lambda is a vector with a list of penalty values.
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
Numeric. An initial value for the estimated parameters with dimensions nvars x nplatforms.
}
  \item{thresh}{
Numeric. The stopping criteria. The default value is 0.1.
}
  \item{maxiter}{
Numeric. Maximum number of iterations. The default value is 100.
}
  \item{methods}{
Character ("lass" or "scad"). \code{lass}: LASSO; \code{scad}: SCAD.
}
  \item{link}{
Character ("logit" or "probit"). Link functions: logistic or probit.
}
  \item{Complete}{
Logic input. If \code{Complete == TRUE}, the predictors \eqn{M_1},...,\eqn{M_p} are measured in all platforms. If \code{Compelte == FALSE}, in some platforms, not all of the predictors \eqn{\{M_1,M_2,...,M_p\}} are measured. The values of the corresponding estimated coefficients for the missing predictors will be \code{NA}.
}
  \item{depen}{
Character. Input only for function \code{fusionbinary.fit}. "IND" means the observations across different platforms are independent; "CORR" means the observations are correlated, and the sample sizes should be equal for different platforms. 
}
  \item{a}{
Numeric. Input only for function \code{fusionbinary.fit}. The free multiplicative constant used in \eqn{\gamma_n}. The default value is 1.
}
}
\details{
The generalized fusion learning function to learn from multiple models with binary responses. More details regarding the algorithm can be found in \code{\link{FusionLearn}}.
}
\value{

\code{fusionbinary} returns a list that has components:

  \item{beta}{
A matrix (nvars x nplatforms) containing estimated coefficients of each linear model. If some data sets do not have the complete set of predictors, the corresponding coefficients are output as \code{NA}.
}

  \item{method}{
Penalty function LASSO or SCAD.
}

  \item{link}{
The link function used in the estimation.
}

  \item{threshold}{
The numeric value shows the difference in the estimates between the successive updates upon convergence. 
}

  \item{iteration}{
The numeric value shows the number of iterations upon convergence.
}

\code{fusionbinary.fit} provides the results in a table:

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

\code{fusionbinary.fit} also returns a model selection plot showing the results above. 
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
##Analysis of the gene data 
y = list(mockgene1[,2],mockgene2[,2])           ## responses "status"
x = list(mockgene1[,3:502],mockgene2[,3:502])   ## 500 predictors 


##Implementing fusion learning algorithm 
result <- fusionbinary(x,y,0.3,N=c(98,286),500,2) 
id <- which(result$beta[,1]!=0)+2
genename <- colnames(mockgene1)[id]

}


