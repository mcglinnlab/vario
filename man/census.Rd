\name{census}
\alias{census}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
census(sim, snap = length(sim$snaps), type = "census")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{
%%     ~~Describe \code{sim} here~~
}
  \item{snap}{
%%     ~~Describe \code{snap} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (sim, snap = length(sim$snaps), type = "census") 
{
    pop <- sim$snaps[[snap]]
    dim(pop) <- c(sim$p$S, sim$p$M, sim$p$M)
    if (type == "census") {
        output <- t(pop[, , 1])
        for (i in 2:sim$p$M) output <- rbind(output, t(pop[, 
            , i]))
    }
    else if (type == "richness") 
        output <- apply(pop, MARGIN = 2:3, FUN = function(x) sum(as.logical(x)))
    else if (type == "abundance") 
        output <- apply(pop, MARGIN = 2:3, FUN = sum)
    else stop("Invalid type")
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
