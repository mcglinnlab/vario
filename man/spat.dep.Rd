\name{spat.dep}
\alias{spat.dep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spat.dep(x, g, sp.dep = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{sp.dep}{
%%     ~~Describe \code{sp.dep} here~~
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
function (x, g, sp.dep = FALSE) 
{
    y <- x
    p <- apply(y, 2, sum)/dim(y)[1]
    for (k in 1:g) {
        for (j in 1:nrow(y)) {
            if (j == 1) {
                y[j, ] <- y[nrow(y), ] + y[1, ] + y[2, ]
            }
            else {
                if (j == nrow(y)) {
                  y[j, ] <- y[j - 1, ] + y[j, ] + y[1, ]
                }
                else {
                  y[j, ] <- y[j - 1, ] + y[j, ] + y[j + 1, ]
                }
            }
        }
    }
    y <- sapply(1:ncol(y), function(x) ifelse(y[, x] >= quantile(y[, 
        x], 1 - p[x]), 1, 0))
    if (sp.dep) {
        for (i in 1:ncol(y)) {
            if (runif(1) >= 0.5) {
                y[, i] <- y[N:1, i]
            }
        }
    }
    else {
        for (i in 1:ncol(y)) {
            ind <- sample(nrow(y), 1):nrow(y)
            if (length(ind) < nrow(y)) 
                ind <- c(ind, 1:(nrow(y) - length(ind)))
            y[, i] <- y[ind, i]
        }
    }
    y
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
