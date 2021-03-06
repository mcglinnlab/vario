% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vario.R
\name{estimateD}
\alias{estimateD}
\title{Estimates the fractal dimension D from a 2-dimensional grid}
\usage{
estimateD(x, coord = NULL)
}
\arguments{
\item{x}{is a matrix of real numbers}

\item{coord}{the spatial coordinates}
}
\description{
Estimates the fractal dimension D from a 2-dimensional grid
}
\examples{
x = matrix(runif(100), 10, 10)
estimateD(x)
}
