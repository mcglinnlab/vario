\name{Empir}
\alias{Empir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Empir(stcrds, xmn, ymn, len, n, pos.neg = TRUE, meth = "randpat", nstrata = 3, rperm = 1000, nperm = 299, npar = 8, all.dat = psp.bbs, linux = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stcrds}{
%%     ~~Describe \code{stcrds} here~~
}
  \item{xmn}{
%%     ~~Describe \code{xmn} here~~
}
  \item{ymn}{
%%     ~~Describe \code{ymn} here~~
}
  \item{len}{
%%     ~~Describe \code{len} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{pos.neg}{
%%     ~~Describe \code{pos.neg} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{nstrata}{
%%     ~~Describe \code{nstrata} here~~
}
  \item{rperm}{
%%     ~~Describe \code{rperm} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{npar}{
%%     ~~Describe \code{npar} here~~
}
  \item{all.dat}{
%%     ~~Describe \code{all.dat} here~~
}
  \item{linux}{
%%     ~~Describe \code{linux} here~~
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
function (stcrds, xmn, ymn, len, n, pos.neg = TRUE, meth = "randpat", 
    nstrata = 3, rperm = 1000, nperm = 299, npar = 8, all.dat = psp.bbs, 
    linux = TRUE) 
{
    for (icord in 1:length(stcrds[, 1])) {
        xmin <- stcrds[icord, 1]
        xmax <- xmin + n * len
        ymin <- stcrds[icord, 2]
        ymax <- ymin + n * len
        xs <- seq(xmin, xmax, len)
        ys <- seq(ymin, ymax, len)
        psp <- all.dat[, ((xmin - xmn + len)/len):((xmax - xmn)/len), 
            ((ymin - ymn + len)/len):((ymax - ymn)/len)]
        psp <- psp[apply(psp, 1, sum) > 0, , ]
        S <- dim(psp)[1]
        mat <- matrix(as.vector(psp), ncol = S, nrow = n^2, byrow = TRUE)
        coords <- expand.grid(xs[-(n + 1)], ys[-(n + 1)])
        v <- vario(mat, coords, pos.neg = pos.neg, grain = len)
        if (meth == "randpat") 
            rpargs <- list(TRUE, nstrata, rperm, 0, 0.01, 1)
        else rpargs <- list(FALSE, nstrata, rperm, 0, 0.01, 1)
        out <- null.perms(x = psp, vobject = v, coords = coords, 
            nperm = nperm, meth = meth, npar = npar, RPargs = rpargs, 
            linux = linux)
        save(out, file = paste("empir.n", n, ".", meth, ".x", 
            xmin, ".y", ymin, ".Rdata", sep = ""))
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
