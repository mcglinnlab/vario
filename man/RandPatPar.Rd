\name{randPatPar}
\alias{randPatPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
randPatPar(psp, nstrata, mtrials1 = 1000, mtrials2 = 1e+06, alpha = 0.01, npar = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psp}{
%%     ~~Describe \code{psp} here~~
}
  \item{nstrata}{
%%     ~~Describe \code{nstrata} here~~
}
  \item{mtrials1}{
%%     ~~Describe \code{mtrials1} here~~
}
  \item{mtrials2}{
%%     ~~Describe \code{mtrials2} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{npar}{
%%     ~~Describe \code{npar} here~~
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
function (psp, nstrata, mtrials1 = 1000, mtrials2 = 1e+06, alpha = 0.01, 
    npar = 1) 
{
    n = dim(psp)[2]
    if (length(dim(psp)) == 3) 
        S = dim(psp)[1]
    else {
        S = 1
        psp = array(psp, dim = c(S, n, n))
    }
    n2 = n + 2
    psp.temp = array(0, dim = c(S, n2, n2))
    psp.temp[, -c(1, n2), -c(1, n2)] = psp
    psp.temp[, 1, -c(1, n2)] = -999
    psp.temp[, n2, -c(1, n2)] = -999
    psp.temp[, -c(1, n2), 1] = -999
    psp.temp[, -c(1, n2), n2] = -999
    psp.temp[, 1, 1] = -999
    psp.temp[, 1, n2] = -999
    psp.temp[, n2, 1] = -999
    psp.temp[, n2, n2] = -999
    psp = psp.temp
    pl = 1:n2^2
    c1 = NA
    c2 = NA
    skip = pl[-999 == as.vector(psp[1, , ])]
    pl = pl[is.na(match(pl, skip))]
    rpsp = psp
    rpsp[, -c(1, n2), -c(1, n2)] = spatPermStrata(psp[, -c(1, 
        n2), -c(1, n2)], meth = "reflect", sp = TRUE, nstrata = nstrata)
    nplaces = length(pl)
    nswaps = (nplaces * (nplaces - 1))/2
    if (npar > 1) {
        require(snowfall)
        sfInit(parallel = TRUE, cpus = npar, type = "SOCK")
        sfClusterSetupRNG()
        sfExport("randPat", "fixUnSampFalseBorder", "psp", "rpsp", 
            "n", "nstrata", "pl", "mtrials1", "mtrials2", "alpha")
        sfClusterEval(dyn.load("randompatternspar.dll"))
        out = unlist(sfSapply(1:S, function(i) {
            randPat(i = i, psp = psp, rpsp = rpsp, n = n, nstrata = nstrata, 
                pl = pl, mtrials1 = mtrials1, mtrials2 = mtrials2, 
                alpha = alpha)
        }))
        sfStop()
    }
    else {
        out = NULL
        for (i in 1:S) out = cbind(out, randPat(i = i, psp = psp, 
            rpsp = rpsp, n = n, nstrata = nstrata, pl = pl, mtrials1 = mtrials1, 
            mtrials2 = mtrials2, alpha = alpha))
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
