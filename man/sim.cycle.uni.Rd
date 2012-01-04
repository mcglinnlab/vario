\name{sim.cycle.uni}
\alias{sim.cycle.uni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.cycle.uni(sim, time = 500, cycles = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{
%%     ~~Describe \code{sim} here~~
}
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{cycles}{
%%     ~~Describe \code{cycles} here~~
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
function (sim, time = 500, cycles = 20) 
{
    for (i in 1:cycles) {
        cat(paste("Cycle", i, ";", sep = " "))
        sim <- sim.iter.uni(sim, time)
    }
    return(sim)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
