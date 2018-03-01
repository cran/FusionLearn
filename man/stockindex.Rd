\name{stockindex}
\alias{stockindexVIX}
\alias{stockindexGSPC}
\alias{stockindexDJI}
\docType{data}
\title{
Finance Data
}
\description{
This is a dataset containing the log return on three financial market indices and 46 stocks between 2013 and 2015. The responses are the financial indices, "VIX", "SP500", and "DJI", and the predictors are 46 stocks from the market. The data are given in three-day gap from 700 trading days.
}
\usage{data("stockindexVIX")
data("stockindexGSPC")
data("stockindexDJI")}

\details{
This example is used to demonstrate the use of the functions \code{\link{fusionbase}} and \code{\link{fusionmixed}}. This dataset has correlated responses and the same predictors values for three models.
}

\source{
This data was obtained from \url{https://finance.yahoo.com}. 
}


