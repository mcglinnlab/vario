\name{neut.sim.uni}
\alias{neut.sim.uni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
neut.sim.uni(M = 7, K = 500, S = 50, pop = round(rep(max(K)/S, M * M * S)), s = 0.1, u = 0.01, d = 0.5, b = 0.505, m = 0.5 * u, habitat = 1, fitness = 1, ef.strength = 0, recruit = 1, fec = 1, time = 100, cycles = 20, mig = 0, grad.len = NULL)
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
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{cycles}{
%%     ~~Describe \code{cycles} here~~
}
  \item{mig}{
%%     ~~Describe \code{mig} here~~
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
function (M = 7, K = 500, S = 50, pop = round(rep(max(K)/S, M * 
    M * S)), s = 0.1, u = 0.01, d = 0.5, b = 0.505, m = 0.5 * 
    u, habitat = 1, fitness = 1, ef.strength = 0, recruit = 1, 
    fec = 1, time = 100, cycles = 20, mig = 0, grad.len = NULL) 
{
    ptm <- proc.time()
    sim <- sim.init.uni(M = M, K = K, S = S, pop = pop, s = s, 
        u = u, d = d, b = b, m = m, habitat = habitat, fitness = fitness, 
        ef.strength = ef.strength, recruit = recruit, fec = fec, 
        grad.len = grad.len)
    if (mig == 0) 
        sim <- sim.cycle(sim = sim, time = time, cycles = cycles)
    if (mig == 1) 
        sim <- sim.cycle.uni(sim = sim, time = time, cycles = cycles)
    print("system.time")
    print(proc.time() - ptm)
    return(sim)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
