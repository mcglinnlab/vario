\name{GridCount}
\alias{GridCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GridCount(samps, grd, land = NULL, perc.land = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samps}{
%%     ~~Describe \code{samps} here~~
}
  \item{grd}{
%%     ~~Describe \code{grd} here~~
}
  \item{land}{
%%     ~~Describe \code{land} here~~
}
  \item{perc.land}{
%%     ~~Describe \code{perc.land} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (samps, grd, land = NULL, perc.land = NULL) 
{
    dat.x <- samps[, 1]
    dat.y <- samps[, 2]
    grd.x <- grd[, 1]
    grd.y <- grd[, 2]
    uni.x <- sort(unique(grd.x))
    uni.y <- sort(unique(grd.y))
    luni.x <- length(uni.x)
    luni.y <- length(uni.y)
    if (is.null(land)) {
        land <- matrix(0, nrow = luni.x - 1, ncol = luni.y - 
            1)
        perc.land <- 0
    }
    n <- matrix(0, nrow = luni.x - 1, ncol = luni.y - 1)
    for (i in 1:(luni.x - 1)) {
        for (j in 1:(luni.y - 1)) {
            if (land[i, j] >= perc.land) {
                n[i, j] <- sum(dat.x >= uni.x[i] & dat.x < uni.x[i + 
                  1] & dat.y >= uni.y[j] & dat.y < uni.y[j + 
                  1])
            }
        }
    }
    out <- list(x = uni.x, y = uni.y, n = n)
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
