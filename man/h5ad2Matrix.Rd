% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h5ad_util.R
\name{h5ad2Matrix}
\alias{h5ad2Matrix}
\title{Parse matrix from h5ad file}
\usage{
h5ad2Matrix(filename, name)
}
\arguments{
\item{filename}{file name or H5IdComponent to read data.frame from}

\item{name}{name of group that contains matrix to be loaded}
}
\value{
R matrix, dense or sparse - in dependence on input
}
\description{
Parse matrix from h5ad file
}
\examples{
obs = h5ad2data.frame('adata.h5ad','X')
}
