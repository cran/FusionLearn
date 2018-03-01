\name{mockgene}
\alias{mockgene1}
\alias{mockgene2}
\docType{data}
\title{
Mock Gene Data
}
\description{
This dataset is a mock version of two different microarray experiments on breast cancer cells.
}
\usage{data("mockgene1")
data("mockgene2")}
\format{
The first data "mockgene1" contains 98 subjects, and the second data "mockgene2" contains 286 subjects.

The first column for each data is ID number.

The second column is subjects' status. If the status is the estrogen-receptor-positive, \eqn{y=1}; if the status is estrogen-receptor-negative, \eqn{y=0}.
Other columns record the gene expression values.

}
\details{
This is an example to implement the \code{FusionLearn} algorithm for binary responses. In this case, the two experiments followed different protocols, and the two sets of gene expression profiles are different. The objective is to select a suitable subset gene predictors for the disease analysis based on both experiments. 
}
\source{
This data is a mock version of the original data. The original gene data contain over 20,000 profile expressions, and more details can be found on \url{https://www.ncbi.nlm.nih.gov/} with series numbers GSE2034 and GSE22093.
}

