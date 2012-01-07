\name{estimateD}
\alias{estimateD}
\title{Provides an estimate of the Fractal Dimension, D }
\description{
Esimates the fractal dimension D from a 2-dimensional grid
}
\usage{
estimateD(z)
}
\arguments{
  \item{z}{ A matrix of real numbers }
}
\details{
Assumes that variance increases in a linear fashion in log-log space
}
\value{
A real number between 2 and 3
}
\references{
Palmer, M. W. 1988. Fractal geometry: a tool for describing spatial patterns of
plant communities. Vegetatio 74:91-102.
}
\author{
Dan McGlinn
}
\note{
Future versions of this function will allow the user to specify coordinates as
an input
}
\seealso{
\code{\link{vario}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z) 
{
    n = dim(z)[1]
    gcords = expand.grid(1:n, 1:n)
    v = vario(as.vector(z), gcords)$vario
    mod = lm(log(v$exp) ~ log(v$Dist))
    m = coef(mod)[2]
    D = (6 - m)/2
    return(D)
  }
}

