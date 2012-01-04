\name{RandPat}
\alias{RandPat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RandPat(i, psp, rpsp, n, nstrata, pl, mtrials1 = 1000, mtrials2 = 1e+06, alpha = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{psp}{
%%     ~~Describe \code{psp} here~~
}
  \item{rpsp}{
%%     ~~Describe \code{rpsp} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nstrata}{
%%     ~~Describe \code{nstrata} here~~
}
  \item{pl}{
%%     ~~Describe \code{pl} here~~
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
function (i, psp, rpsp, n, nstrata, pl, mtrials1 = 1000, mtrials2 = 1e+06, 
    alpha = 0.01) 
{
    psp <- psp[i, , ]
    rpsp <- rpsp[i, , ]
    n2 <- n + 2
    rpsp.tmp <- rpsp[-c(1, n2), -c(1, n2)]
    coords <- cbind(rep(1:nstrata, nstrata), rep(1:nstrata, each = nstrata))
    rcoords <- coords[sample(nstrata^2), ]
    for (j in 1:nstrata^2) {
        rows <- ((coords[j, 1] - 1) * n/nstrata + 1):((coords[j, 
            1] - 1) * n/nstrata + n/nstrata)
        cols <- ((coords[j, 2] - 1) * n/nstrata + 1):((coords[j, 
            2] - 1) * n/nstrata + n/nstrata)
        rrows <- ((rcoords[j, 1] - 1) * n/nstrata + 1):((rcoords[j, 
            1] - 1) * n/nstrata + n/nstrata)
        rcols <- ((rcoords[j, 2] - 1) * n/nstrata + 1):((rcoords[j, 
            2] - 1) * n/nstrata + n/nstrata)
        rpsp[-c(1, n2), -c(1, n2)][rows, cols] <- rpsp.tmp[rrows, 
            rcols]
    }
    rpsp <- FixUnSamp(psp, rpsp)
    ostat <- .C("spatstat", as.double(as.vector(psp)), as.integer(n), 
        as.double(rep(0, 4)))[[3]]
    nstat <- .C("spatstat", as.double(as.vector(rpsp)), as.integer(n), 
        as.double(rep(0, 4)))[[3]]
    phi <- .C("calcphi", as.double(nstat), as.double(ostat), 
        as.double(0))[[3]]
    ntrials <- 0
    gtrials <- 0
    rpsp.tmp1 <- rpsp
    while (phi > alpha & ntrials < mtrials1) {
        rpsp.tmp2 <- rpsp.tmp1[-c(1, n2), -c(1, n2)]
        rcoords <- coords[sample(nstrata^2, 2), ]
        startrows <- ((rcoords[1, 1] - 1) * n/nstrata + 1):((rcoords[1, 
            1] - 1) * n/nstrata + n/nstrata)
        startcols <- ((rcoords[1, 2] - 1) * n/nstrata + 1):((rcoords[1, 
            2] - 1) * n/nstrata + n/nstrata)
        endrows <- ((rcoords[2, 1] - 1) * n/nstrata + 1):((rcoords[2, 
            1] - 1) * n/nstrata + n/nstrata)
        endcols <- ((rcoords[2, 2] - 1) * n/nstrata + 1):((rcoords[2, 
            2] - 1) * n/nstrata + n/nstrata)
        rpsp.tmp1[-c(1, n2), -c(1, n2)][startrows, startcols] <- rpsp.tmp2[endrows, 
            endcols]
        rpsp.tmp1[-c(1, n2), -c(1, n2)][endrows, endcols] <- rpsp.tmp2[startrows, 
            startcols]
        rpsp.tmp1 <- FixUnSamp(psp, rpsp.tmp1)
        nstat <- .C("spatstat", as.double(as.vector(rpsp.tmp1)), 
            as.integer(n), as.double(rep(0, 4)))[[3]]
        phiTemp <- .C("calcphi", as.double(nstat), as.double(ostat), 
            as.double(0))[[3]]
        if (phiTemp < phi) {
            phi <- phiTemp
            gtrials <- gtrials + 1
            rpsp <- rpsp.tmp1
        }
        else {
            rpsp.tmp1 <- rpsp
        }
        ntrials <- ntrials + 1
    }
    if (phi > alpha & mtrials2 > 0) {
        psp <- as.vector(psp)
        rpsp <- as.vector(rpsp)
        tmp <- .C("randpatpar", psp = as.double(psp), rpsp = as.double(rpsp), 
            n = as.integer(n), ostat = as.double(rep(0, 4)), 
            nstat = as.double(rep(0, 4)), phi = as.double(0), 
            phiTemp = as.double(0), alpha = as.double(alpha), 
            pl = as.integer(pl - 1), nplaces = as.integer(length(pl) - 
                1), ntrials = as.double(0), gtrials = as.double(0), 
            mtrials = as.double(mtrials2))
        out <- c(i, phi, gtrials, tmp$phi, tmp$gtrials, tmp$rpsp)
    }
    else out <- c(i, phi, gtrials, NA, NA, as.vector(rpsp))
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
