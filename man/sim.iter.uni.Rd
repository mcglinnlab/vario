\name{sim.iter.uni}
\alias{sim.iter.uni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.iter.uni(sim, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{
%%     ~~Describe \code{sim} here~~
}
  \item{time}{
%%     ~~Describe \code{time} here~~
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
function (sim, time) 
{
    nb = rep(0, sim$p$M * sim$p$M * sim$p$S)
    K = sim$p$K
    if (length(K) == 1) 
        K = rep(K, sim$p$M * sim$p$M * sim$p$S)
    prelength <- length(sim$snaps)
    pretime <- as.numeric(names(sim$snaps)[prelength])
    tmp <- .C("smith08c", M = as.integer(sim$p$M), K = as.integer(K), 
        S = as.integer(sim$p$S), pop = as.integer(sim$snaps[[prelength]]), 
        nb = as.integer(nb), u = as.double(sim$p$u), dv = as.double(sim$p$dv), 
        bv = as.double(sim$p$bv), m = as.double(sim$p$m), t = as.integer(time), 
        ef = as.double(sim$p$ef), fec = as.integer(sim$p$fec))
    sim$snaps[[prelength + 1]] <- tmp$pop
    names(sim$snaps)[prelength + 1] <- as.character(pretime + 
        time)
    return(sim)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
