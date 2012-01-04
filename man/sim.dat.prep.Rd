\name{sim.dat.prep}
\alias{sim.dat.prep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.dat.prep(stcrd, psp, envir, perc.land, res, n, K, K.cont = FALSE, export = FALSE, DESTDIR = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stcrd}{
%%     ~~Describe \code{stcrd} here~~
}
  \item{psp}{
%%     ~~Describe \code{psp} here~~
}
  \item{envir}{
%%     ~~Describe \code{envir} here~~
}
  \item{perc.land}{
%%     ~~Describe \code{perc.land} here~~
}
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{K.cont}{
%%     ~~Describe \code{K.cont} here~~
}
  \item{export}{
%%     ~~Describe \code{export} here~~
}
  \item{DESTDIR}{
%%     ~~Describe \code{DESTDIR} here~~
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
function (stcrd, psp, envir, perc.land, res, n, K, K.cont = FALSE, 
    export = FALSE, DESTDIR = getwd()) 
{
    i <- (stcrd[1]/res) + 1
    j <- (stcrd[2]/res) + 1
    envir.grd <- envir[1, i:(i + n - 1), j:(j + n - 1)]
    psp.grd <- psp[, i:(i + n - 1), j:(j + n - 1)]
    psp.grd <- psp.grd[apply(psp.grd, 1, sum, na.rm = TRUE) > 
        0, , ]
    S.grd <- dim(psp.grd)[1]
    if (K.cont) 
        K <- rep(as.vector(K * envir[2, i:(i + n - 1), j:(j + 
            n - 1)]), each = S.grd)
    else K <- rep(as.vector(ifelse(envir[2, i:(i + n - 1), j:(j + 
        n - 1)] >= perc.land, K, 0)), each = S.grd)
    no.bbs <- ifelse(apply(psp.grd, c(2, 3), sum) == 0, TRUE, 
        FALSE)
    sim.input <- list(M = n, S = S.grd, K = K, no.bbs = no.bbs, 
        envir = envir.grd)
    if (export) 
        save(sim.input, file = paste(DESTDIR, "/SimEmpirInputs.M", 
            n, ".x", stcrd[1], ".y", stcrd[2], ".Rdata", sep = ""))
    else sim.input
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
