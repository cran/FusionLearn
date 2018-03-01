\name{FusionLearn-package}
\alias{FusionLearn}
\docType{package}
\title{
Fusion Learning 
}
\description{
\code{FusionLearn} package implements a new learning algorithm to integrate information from different experimental platforms. The algorithm applies the grouped penalization method in the pseudolikelihood setting. 
}
\details{
In the context of fusion learning, there are \eqn{k} different data sets from \eqn{k} different experimental platforms. The data from each platform can be modeled by a different generalized linear model. Assume the same set of predictors \eqn{\{M_1,M_2,...,M_j,...,M_p \}} are measured across \eqn{k} different experimental platforms.
\tabular{llcccccc}{
Platforms \tab  Formula \tab \eqn{M_1} \tab \eqn{M_2} \tab \eqn{\dots} \tab \eqn{M_j} \tab \eqn{\dots} \tab \eqn{M_p} \cr           
1 \tab \eqn{y_1: g_1(\mu_1) \sim} \tab \eqn{x_{11}\beta_{11}+} \tab \eqn{x_{12}\beta_{12}+} \tab \eqn{\dots}  \tab \eqn{x_{1j}\beta_{1j}+} \tab \eqn{\dots} \tab \eqn{x_{1p}\beta_{1p}} \cr
2 \tab \eqn{y_2: g_2(\mu_2) \sim} \tab \eqn{x_{21}\beta_{21}+} \tab \eqn{x_{22}\beta_{22}+} \tab \eqn{\dots} \tab \eqn{x_{2j}\beta_{2j}+} \tab \eqn{\dots}  \tab \eqn{x_{2p}\beta_{2p}}  \cr
\tab  ... \tab  \cr
k \tab \eqn{y_k: g_k(\mu_k) \sim} \tab \eqn{x_{k1}\beta_{k1}+} \tab \eqn{x_{k2}\beta_{k2}+} \tab \eqn{\dots} \tab \eqn{x_{kj}\beta_{kj}+} \tab \eqn{\dots} \tab \eqn{x_{kp}\beta_{kp}}  \cr
}
Here \eqn{x_{kj}} represents the observation of the predictor \eqn{M_j} on the \eqn{k}th platform, and \eqn{\beta^{(j)}} denotes the vector of regression coefficients for the predictor \eqn{M_j}. 
\tabular{llll}{
Platforms \tab \eqn{\bold{M_j}} \tab \tab \eqn{\bold{\beta^{(j)}}} \cr
1 \tab \eqn{x_{1j}} \tab \tab \eqn{\beta_{1j}} \cr
2 \tab \eqn{x_{2j}} \tab \tab \eqn{\beta_{2j}} \cr
 \tab ... \tab \tab ... \cr
k \tab \eqn{x_{kj}} \tab \tab \eqn{\beta_{kj}}
}

Consider the following examples.

\emph{Example 1.} Suppose \eqn{k} different types of experiments are conducted to study the genetic mechanism of a disease. The predictors in this research are different facets of individual genes, such as mRNA expression, protein expression, RNAseq expression and so on. The goal is to select the genes which affect the disease, while the genes are assessed in a number of ways through different measurement processes across \eqn{k} experimental platforms.

\emph{Example 2.} The predictive models for three different financial indices are simultaneously built from a panel of stock index predictors. In this case, the predictor values across different models are the same, but the regression coefficients are different.

In the conventional approach, the model for each of the \eqn{k} platforms is analyzed separately. \code{FusionLearn} algorithm selects significant predictors through learning from multiple models. The overall objective is to minimize the function: 
\deqn{Q(\beta)=l_I(\beta)- n \sum_{j=1}^{p} \Omega_{\lambda_n} ||\beta^{(j)}||,}
with \eqn{p} being the numbers of predictors, \eqn{\Omega_{\lambda_n}} being the penalty functions, and \eqn{||\beta^{(j)}|| = (\sum_{i=1}^{k}\beta_{ij}^2)^{1/2}} denoting the \eqn{L_2}-norm of the coefficients of the predictor \eqn{M_j}. 

The user can specify the penalty function \eqn{\Omega_{\lambda_n}} and the penalty values \eqn{\lambda_n}. This package also contains functions to provide the pseudolikelihood Bayesian information criterion:

\deqn{ pseu-BIC(s) = -2l_I(\hat{\beta}_I;Y) + d_s^{*} \gamma_n }

with \eqn{-2l_I(\hat{\beta}_I; Y)} denoting the pseudo loglikelihood, \eqn{d_s^{*}} measuring the model complexity and \eqn{\gamma_n} being the penalty on the model complexity. 

The basic function \code{\link{fusionbase}} deals with continuous responses. The function \code{\link{fusionbinary}} is applied to binary responses, and the function \code{\link{fusionmixed}} is applied to a mix of continuous and binary responses.
}

\note{
Here we provide two examples to illustrate the data structures. Assume \eqn{X_I} and \eqn{X_{II}} represent two sets of the predictors from 2 experimental platforms. 

\emph{Example 1}. If the observations from \eqn{X_I} and \eqn{X_{II}} are independent, the number of observations can be different. The order of the predictors \eqn{\{M_1, M_2, M_3, M_4\}} in \eqn{X_I} matches with the predictors in \eqn{X_{II}}. If \eqn{X_{II}} does not include the predictor \eqn{M_3}, then the \eqn{M_3} in \eqn{X_{II}} needs to be filled with \code{NA}.    
\tabular{lrrrrlrrrr}{
     \tab    \eqn{M_1} \tab \eqn{M_2} \tab \eqn{M_3} \tab \eqn{M_4} \tab
     \tab    \eqn{M_1} \tab \eqn{M_2} \tab \eqn{M_3} \tab \eqn{M_4} \cr 
\eqn{X_I = } \tab 0.1  \tab     0.3   \tab   0.5     \tab   20      \tab
\eqn{X_{II} = } \tab 100  \tab     8   \tab  \code{NA} \tab     100      \cr
             \tab 0.3  \tab     0.1   \tab   0.5     \tab   7       \tab
             \tab 30  \tab      1   \tab   \code{NA} \tab     2     \cr
             \tab 0.1  \tab     0.9   \tab     1     \tab   0       \tab
             \tab 43  \tab     19   \tab  \code{NA}   \tab   -3       \cr
             \tab -0.3  \tab    1.2   \tab     2     \tab     40  }
\emph{Example 2}. If the observations from \eqn{X_I} and \eqn{X_{II}} are correlated, the number of observations must be the same. The \eqn{i}th row in \eqn{X_I} is correlatd with the \eqn{i}th row in \eqn{X_{II}}. The predictors of \eqn{X_I} and \eqn{X_{II}} should be matched in order. The predictors which are not measured need to be filled with \code{NA}.
\tabular{lrrrrlrrrr}{
    \tab \eqn{M_1} \tab \eqn{M_2} \tab \eqn{M_3} \tab \eqn{M_4} \tab
     \tab \eqn{M_1} \tab \eqn{M_2} \tab \eqn{M_3} \tab \eqn{M_4} \cr
\eqn{X_I = } \tab 0.1  \tab     0.3   \tab   0.5     \tab   20      \tab
\eqn{X_{II} = } \tab 0.3  \tab     0.8   \tab  \code{NA} \tab     100      \cr
      \tab 0.3  \tab     0.1   \tab   0.5     \tab   70       \tab
      \tab 0.2  \tab      1   \tab   \code{NA} \tab    20     \cr
       \tab -0.1  \tab     0.9   \tab     1     \tab   0       \tab
             \tab  0.43  \tab   1.9   \tab  \code{NA}   \tab   -30  \cr
             \tab -0.3  \tab    1.2   \tab     2     \tab     40   \tab
             \tab  -0.4  \tab   -2   \tab    \code{NA} \tab   40}
In functions \code{fusionbase.fit}, \code{fusionbinary.fit}, and \code{fusionmixed.fit}, the option \code{depen} is used to specify whether observations from different platforms are correlated or independent. 
}
\author{
Xin Gao, Yuan Zhong and Raymond J Carroll

Maintainer: Yuan Zhong <aqua.zhong@gmail.com>
}
\references{
Gao, X and Carroll, R. J. (2017) Data integration with high dimensionality. Biometrika, 104, 2, pp. 251-272
}




