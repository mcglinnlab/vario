\name{grid.SAR}
\alias{grid.SAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grid.SAR(psp, grains, single = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psp}{
%%     ~~Describe \code{psp} here~~
}
  \item{grains}{
%%     ~~Describe \code{grains} here~~
}
  \item{single}{
%%     ~~Describe \code{single} here~~
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
function (psp, grains, single = FALSE) 
{
    grains <- grains[sqrt(grains) == round(sqrt(grains))]
    lens <- sqrt(grains)
    ns <- dim(psp)[2]
    sr <- rep(0, length(lens))
    cs <- rep(0, length(lens))
    for (l in 1:length(lens)) {
        if (lens[l] == 1) {
            sr[l] <- sum(psp)/n^2
            cs[l] <- 1
        }
        else {
            for (col in 1:(ns - lens[l] + 1)) {
                for (row in 1:(ns - lens[l] + 1)) {
                  if (single) {
                    sr[l] <- sr[l] + sum(sum(psp[row:(row + (lens[l] - 
                      1)), col:(col + (lens[l] - 1))]) > 0)
                    cs[l] <- cs[l] + 1
                  }
                  else {
                    sr[l] <- sr[l] + sum(apply(psp[, row:(row + 
                      (lens[l] - 1)), col:(col + (lens[l] - 1))], 
                      1, sum) > 0)
                    cs[l] <- cs[l] + 1
                  }
                }
            }
        }
    }
    out <- cbind(grains, sr/cs)
    colnames(out) <- c("grains", "richness")
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
