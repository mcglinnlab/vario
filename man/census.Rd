\name{census}
\alias{census}
\title{Extract census data from sim objects }
\description{
Returns a summary of census data from a sim object, optionally the raw
census, the total abundance of all individuals in each cell, or the
species richness in each cell. By default, returns the last snapshot.
}
\usage{
census(sim, snap = length(sim$snaps), type = "census")
}
\arguments{
  \item{sim}{An object of class sim }
  \item{snap}{Either an integer specifying the index of the snapshot to
    extract, or the character name of the snapshot. i.e., \dQuote{500}
    returns the snapshot corresponding to the 500th timestep (if it exists). }
  \item{type}{ One of \dQuote{census}, \dQuote{richness},
    \dQuote{abundance}, specifying which summary value is desired.}
}
\details{
Pretty much as described above.
}
\value{
  For type == \dQuote{census}, returns a matrix with one row for every
  grid cell, and one column for every species, with individual matrix
  element [i,j] containing the abundance of species j in cell i. The
  census order is the same as the order of the \sQuote{coords} element
  in the sim object, for use in rda or variation partitioning.

  For type == \dQuote{abundance}, returns a square matrix of dimension M
  x M, with each matrix element ij containing the total abundance of
  individuals in grid cell ij. Suitable for passing to image(),
  contour() or related functions for display.

  For type == \dQuote{richness}, returns a square matrix of dimension M
  x M, with each matrix element ij containing the species richness of
  grid cell ij. Suitable for passing to image(), contour() or related
  functions for display.

}
\author{Tyler Smith }
\seealso{ \code{\link{neut.simulate}}, \code{\link{sim.cut}}}
\examples{
## See examples in \code{\link{neut.simulate}}
}
