\name{nullGen}
\alias{nullGen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nullGen(pop, vobject, coords, meth, sp, all = FALSE, RPargs = FALSE, median = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop}{
%%     ~~Describe \code{pop} here~~
}
  \item{vobject}{
%%     ~~Describe \code{vobject} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{sp}{
%%     ~~Describe \code{sp} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
}
  \item{RPargs}{
%%     ~~Describe \code{RPargs} here~~
}
  \item{median}{
%%     ~~Describe \code{median} here~~
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
function (pop, vobject, coords, meth, sp, all = FALSE, RPargs = FALSE, 
    median = FALSE) 
{
    S = dim(pop)[1]
    n = dim(pop)[2]
    n2 = n + 2
    grain = vobject$parms$grain
    hmax = vobject$parms$hmax
    pos.neg = vobject$parms$pos.neg
    median = vobject$parms$median
    if (class(vobject$parms$direction) == "factor") 
        direction = as.character(vobject$parms$direction)
    else direction = as.numeric(vobject$parms$direction)
    tolerance = vobject$parms$tolerance
    unit.angle = as.character(vobject$parms$unit.angle)
    if (all) {
        rmat = apply(pop, 1, as.vector)
        rv = vobject$vario
        if (pos.neg) {
            if (median) 
                r.vals = array(NA, dim = c(length(rv$Dist), 6, 
                  3))
            else r.vals = array(NA, dim = c(length(rv$Dist), 
                3, 3))
        }
        else r.vals = array(NA, dim = c(length(rv$Dist), 2, 3))
        for (j in 1:3) {
            if (j == 1) {
                rpop = spatPerm2D(pop, meth = "random", sp = TRUE)
                rpop = fixUnSampTrueBorder(pop, rpop)
            }
            if (j == 2) {
                if (RPargs[[1]]) {
                  out = randPatPar(psp = pop, nstrata = RPargs[[2]], 
                    mtrials1 = RPargs[[3]], mtrials2 = RPargs[[4]], 
                    alpha = RPargs[[5]], npar = RPargs[[6]])
                  rpop = array(0, dim = c(S, n2, n2))
                  for (i in 1:S) {
                    rpop[i, , ] = array(out[-(1:5), i], dim = c(n2, 
                      n2))
                  }
                  rpop = rpop[, -c(1, n2), -c(1, n2)]
                }
                else rpop = spatPerm2D(pop, meth = "both", sp = TRUE)
            }
            if (j == 3) {
                rpop = spatPerm2D(pop, meth = "random", sp = FALSE)
                rpop = fixUnSampTrueBorder(pop, rpop)
            }
            rmat = apply(rpop, 1, as.vector)
            rv = vario(x = rmat, coord = coords, grain = grain, 
                hmax = hmax, pos.neg = pos.neg, median = median, 
                direction = direction, tolerance = tolerance, 
                unit.angle = unit.angle)$vario
            if (pos.neg) {
                if (median) 
                  r.vals[, , j] = as.matrix(rv[, c(5, 7:11)])
                else r.vals[, , j] = as.matrix(rv[, c(5, 7:8)])
            }
            else r.vals[, , j] = as.matrix(rv[, 4:5])
        }
    }
    else {
        if (RPargs[[1]]) {
            out = randPatPar(psp = pop, nstrata = RPargs[[2]], 
                mtrials1 = RPargs[[3]], mtrials2 = RPargs[[4]], 
                alpha = RPargs[[5]], npar = RPargs[[6]])
            rpop = array(0, dim = c(S, n2, n2))
            for (i in 1:S) {
                rpop[i, , ] = array(out[-(1:5), i], dim = c(n + 
                  2, n + 2))
            }
            rpop = rpop[, -c(1, n + 2), -c(1, n + 2)]
        }
        else {
            rpop = spatPerm2D(pop, meth = meth, sp = sp)
            rpop = fixUnSampTrueBorder(pop, rpop)
        }
        rmat = apply(rpop, 1, as.vector)
        rv = vario(x = rmat, coord = coords, grain = grain, hmax = hmax, 
            pos.neg = pos.neg, median = median, direction = direction, 
            tolerance = tolerance, unit.angle = unit.angle)$vario
        if (pos.neg) {
            if (median) 
                r.vals = as.matrix(rv[, c(5, 7:11)])
            else r.vals = as.matrix(rv[, c(5, 7:8)])
        }
        else r.vals = as.matrix(rv[, 4:5])
    }
    retrun(r.vals)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
