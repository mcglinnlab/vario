\name{Empir.fit}
\alias{Empir.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Empir.fit(stcrds, n, sim.array, meth = "avg", on = "desktop")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stcrds}{
%%     ~~Describe \code{stcrds} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{sim.array}{
%%     ~~Describe \code{sim.array} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{on}{
%%     ~~Describe \code{on} here~~
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
function (stcrds, n, sim.array, meth = "avg", on = "desktop") 
{
    output <- NULL
    for (icord in 1:length(stcrds[, 1])) {
        if (on == "desktop") {
            try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
                n, "/empir.n", n, ".x", stcrds[icord, 1], ".y", 
                stcrds[icord, 2], ".Rdata", sep = "")), silent = TRUE)
        }
        if (on == "laptop") {
            try(load(paste("C:/Users/Daniel McGlinn/Desktop/Lab data/abundance surface/BBS/empir_study_n", 
                n, "/empir.n", n, ".x", stcrds[icord, 1], ".y", 
                stcrds[icord, 2], ".Rdata", sep = "")), silent = TRUE)
        }
        if (!is.null(out)) {
            null.means <- apply(out$vario, c(1, 2), mean)
            null.sds <- apply(out$vario, c(1, 2), sd)
            d.pos <- (out$vario[, 2, 1] - null.means[, 2])/null.sds[, 
                2]
            d.neg <- (out$vario[, 3, 1] - null.means[, 3])/null.sds[, 
                3]
            sim.array <- sim.array[, 1:length(d.pos), , , ]
            obs.array <- array(rbind(d.neg, d.pos), dim = c(2, 
                length(d.pos)))
            obs.array <- array(rep(obs.array, prod(dim(sim.array[1, 
                1, , , ]))), dim = dim(sim.array))
            resid <- obs.array - sim.array
            tot.resid <- apply(abs(resid), 3:5, sum)
            dat <- cbind(rep(unique(D), 25), rep(rep(unique(s), 
                5), each = 3), rep(unique(u), each = 15), tot.resid)
            colnames(dat) <- c("D", "s", "u", "tot.resid")
            best <- dat[order(dat[, 4])[1:5], ]
            worst <- dat[order(dat[, 4], decreasing = TRUE)[1:5], 
                ]
            x <- stcrds[icord, 1] + (n * 100)/2
            y <- stcrds[icord, 2] + (n * 100)/2
            if (meth == "avg") 
                output <- rbind(output, c(1, x, y, apply(best, 
                  2, mean), apply(worst, 2, mean)))
            if (meth == "top") 
                output <- rbind(output, c(1, x, y, best[1, ], 
                  worst[1, ]))
            if (meth == "all") 
                output <- rbind(output, cbind(1:5, rep(x, 5), 
                  rep(y, 5), best, worst))
        }
    }
    colnames(output) <- c("rank", "x", "y", "Dbest", "sbest", 
        "ubest", "resid.best", "Dwor", "swor", "uwor", "resid.wor")
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
