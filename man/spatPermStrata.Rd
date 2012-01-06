\name{spatPermStrata}
\alias{spatPermStrata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spatPermStrata(psp, shiftpos = NULL, rotate = NULL, meth = "shift", sp = FALSE, nstrata = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psp}{
%%     ~~Describe \code{psp} here~~
}
  \item{shiftpos}{
%%     ~~Describe \code{shiftpos} here~~
}
  \item{rotate}{
%%     ~~Describe \code{rotate} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{sp}{
%%     ~~Describe \code{sp} here~~
}
  \item{nstrata}{
%%     ~~Describe \code{nstrata} here~~
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
function (psp, shiftpos = NULL, rotate = NULL, meth = "shift", 
    sp = FALSE, nstrata = 1) 
{
    n = dim(psp)[2]
    if (length(dim(psp)) == 3) {
        S = dim(psp)[1]
        flag = FALSE
    }
    else {
        S = 1
        psp = array(psp, dim = c(S, n, n))
        flag = TRUE
    }
    strata.size = n/nstrata
    if (round(strata.size) != strata.size) {
        stop("Number of strata must be evenly divisable by the linear dimension\n          of the grid")
    }
    Rpsp = psp
    brks = seq(1, n, strata.size)
    for (i in 1:nstrata) {
        sub.i = brks[i]:(brks[i] + strata.size - 1)
        for (j in 1:nstrata) {
            sub.j = brks[j]:(brks[j] + strata.size - 1)
            Rpsp[, sub.i, sub.j] = spatPerm2D(psp[, sub.i, sub.j], 
                meth = meth, sp = sp, shiftpos = shiftpos, rotate = rotate)
        }
    }
    if (flag) 
        Rpsp = drop(Rpsp)
    return(Rpsp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
