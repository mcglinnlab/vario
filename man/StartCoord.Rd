\name{StartCoord}
\alias{StartCoord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
StartCoord(samps, grd, n, nsamp, psamp, land = NULL, perc.land = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samps}{
%%     ~~Describe \code{samps} here~~
}
  \item{grd}{
%%     ~~Describe \code{grd} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nsamp}{
%%     ~~Describe \code{nsamp} here~~
}
  \item{psamp}{
%%     ~~Describe \code{psamp} here~~
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
function (samps, grd, n, nsamp, psamp, land = NULL, perc.land = NULL) 
{
    if (is.null(land)) 
        cnt.map <- GridCount(samps, grd)
    else cnt.map <- GridCount(samps, grd, land, perc.land)
    startcoord <- NULL
    for (i in 1:(length(cnt.map$x) - n)) {
        for (j in 1:(length(cnt.map$y) - n)) {
            prop.fill <- sum(cnt.map$n[i:(i + n - 1), j:(j + 
                n - 1)] >= nsamp)/n^2
            if (prop.fill >= psamp) 
                startcoord <- rbind(startcoord, c(cnt.map$x[i], 
                  cnt.map$y[j], prop.fill))
        }
    }
    startcoord
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
