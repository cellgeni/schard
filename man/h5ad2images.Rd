% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h5ad_util.R
\name{h5ad2images}
\alias{h5ad2images}
\title{Load visium images from h5ad}
\usage{
h5ad2images(filename)
}
\arguments{
\item{filename}{path 2 visium h5ad}
}
\value{
list of slides. Each slide is list with scale.fctors, hires and lowres images
}
\description{
Load visium images from h5ad
}
\examples{
imgs = h5ad2images('adata.h5ad')
}
