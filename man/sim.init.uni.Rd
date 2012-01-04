\name{sim.init.uni}
\alias{sim.init.uni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.init.uni(M = 7, K = 500, S = 50, pop = round(rep(K/S, M * M * S)), s = 0.1, u = 0.01, d = 0.5, b = 0.505, m = 0.5 * u, habitat = 1, fitness = 1, ef.strength = 0, recruit = 1, fec = 1, grad.len = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{pop}{
%%     ~~Describe \code{pop} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{u}{
%%     ~~Describe \code{u} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{habitat}{
%%     ~~Describe \code{habitat} here~~
}
  \item{fitness}{
%%     ~~Describe \code{fitness} here~~
}
  \item{ef.strength}{
%%     ~~Describe \code{ef.strength} here~~
}
  \item{recruit}{
%%     ~~Describe \code{recruit} here~~
}
  \item{fec}{
%%     ~~Describe \code{fec} here~~
}
  \item{grad.len}{
%%     ~~Describe \code{grad.len} here~~
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
function (M = 7, K = 500, S = 50, pop = round(rep(K/S, M * M * 
    S)), s = 0.1, u = 0.01, d = 0.5, b = 0.505, m = 0.5 * u, 
    habitat = 1, fitness = 1, ef.strength = 0, recruit = 1, fec = 1, 
    grad.len = NULL) 
{
    u <- u * M
    if (is.null(grad.len)) 
        s <- s * M
    else s <- s * grad.len
    if (length(habitat) > 1) {
        fd <- array(dim = S * M^2)
        dim(fd) <- c(S, M, M)
        for (i in 1:S) for (j in 1:M) for (k in 1:M) fd[i, j, 
            k] <- (habitat[j, k] - fitness[i])
        fd <- as.vector(fd)
    }
    else {
        fd <- rep(0, S * M * M)
    }
    bv <- gauss.niche.diff(b, fd, s)
    dv <- 1 - gauss.niche.diff(1 - d, fd, s)
    ef <- recruit - fd * ef.strength
    params <- list(K = K, M = M, S = S, s.rel = s, u.rel = u, 
        dv = dv, bv = bv, m = m, habitat = habitat, fitness = fitness, 
        ef = ef, fec = fec)
    coords <- data.frame(X = rep(1:M, each = M), Y = rep(1:M, 
        M))
    ret <- list(p = params, coords = coords, snaps = list(`0` = pop))
    class(ret) <- "sim"
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
