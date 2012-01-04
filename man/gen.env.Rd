\name{gen.env}
\alias{gen.env}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gen.env(M, env = "fractal", D = NA, err = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{env}{
%%     ~~Describe \code{env} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{err}{
%%     ~~Describe \code{err} here~~
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
function (M, env = "fractal", D = NA, err = 0.1) 
{
    if (D == 2) 
        env = "grad"
    if (D == 3) 
        env = "random"
    if (env == "fractal") {
        require(FieldSim)
        if (is.na(D)) 
            stop("must supply fractal dimension (D) when enviornment is a fractal")
        H = 3 - D
        res <- midpoint(H = H, nblevel = ceiling(log2(M)))
        res$Z <- res$Z[-nrow(res$Z), ]
        res$Z <- res$Z[, -ncol(res$Z)]
        ndrop <- 2^(ceiling(log2(M))) - M
        if (ndrop > 0) {
            res$Z <- res$Z[-(1:ndrop), ]
            res$Z <- res$Z[, -(1:ndrop)]
        }
        z <- scale(res$Z)
        z <- uni.trans(z, M)
        D.sim <- 3 - quadvar(z)
        icount <- 1
        while (abs(D.sim - D) >= err) {
            if (is.na(D)) 
                stop("must supply fractal dimension (D) when enviornment is a fractal")
            H = 3 - D
            res <- midpoint(H = H, nblevel = ceiling(log2(M)))
            res$Z <- res$Z[-nrow(res$Z), ]
            res$Z <- res$Z[, -ncol(res$Z)]
            ndrop <- 2^(ceiling(log2(M))) - M
            if (ndrop > 0) {
                res$Z <- res$Z[-(1:ndrop), ]
                res$Z <- res$Z[, -(1:ndrop)]
            }
            z <- scale(res$Z)
            z <- uni.trans(z, M)
            D.sim <- 3 - quadvar(z)
        }
    }
    if (env == "grad") {
        z <- matrix(rep(seq(0, M, length.out = M), times = M), 
            ncol = M, nrow = M)
    }
    if (env == "peaks") {
        x <- rep(seq(0, M, length.out = M), each = M)
        y <- rep(seq(0, M, length.out = M), times = M)
        z <- peaks(x, y, M * c(0.25, 0.75), M * c(0.25, 0.75), 
            M * 0.15)
        z <- uni.trans(z, M)
    }
    if (env == "random") {
        z <- matrix(runif(M^2), M, M)
        z <- uni.trans(z, M)
    }
    z
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
