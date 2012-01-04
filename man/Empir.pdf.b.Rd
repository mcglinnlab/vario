\name{Empir.pdf.b}
\alias{Empir.pdf.b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Empir.pdf.b(stcrds, n, nsamp = 2, both = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stcrds}{
%%     ~~Describe \code{stcrds} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nsamp}{
%%     ~~Describe \code{nsamp} here~~
}
  \item{both}{
%%     ~~Describe \code{both} here~~
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
function (stcrds, n, nsamp = 2, both = FALSE) 
{
    if (both) 
        pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
            n, "_random/empir_n", n, "_random.pdf", sep = ""), 
            width = 7 * 3, height = 7 * 2)
    else pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
        n, "_random/empir_n", n, "_random.pdf", sep = ""), width = 7 * 
        2.5, height = 7)
    icount <- 0
    for (icord in 1:length(stcrds[, 1])) {
        for (i in 1:2) {
            out <- NULL
            if (i == 1) {
                try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
                  n, "_random/empir.n", n, ".random.x", stcrds[icord, 
                    1], ".y", stcrds[icord, 2], ".Rdata", sep = "")), 
                  silent = TRUE)
                if (both) 
                  par(mfrow = c(2, 3))
                else par(mfrow = c(1, 3))
                brks <- seq(nsamp, 85, 5)
                cls <- terrain.colors(length(brks) - 1)
                image(Bn100$x, Bn100$y, Bn100$n, breaks = brks, 
                  col = cls, ylim = c(-100, 3500), xlim = c(-100, 
                    5500), main = paste("min samp", nsamp), xlab = "", 
                  ylab = "")
                points(stcrds[icord, 1], stcrds[icord, 2], col = "purple", 
                  pch = 19, cex = 3)
                polygon(c(rep(stcrds[icord, 1], 2), rep(stcrds[icord, 
                  1] + n * len1, 2)), c(stcrds[icord, 2], rep(stcrds[icord, 
                  2] + n * len1, 2), stcrds[icord, 2]))
                lines(wrld.pts, col = "red", lwd = 2)
                v.graph.all(out, plot = FALSE, draw = "exp", 
                  main = paste("geo coords: x=", stcrds[icord, 
                    1], " y=", stcrds[icord, 2], sep = ""))
            }
            if (i == 2) {
                try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
                  n, "/empir.n", n, ".x", stcrds[icord, 1], ".y", 
                  stcrds[icord, 2], ".Rdata", sep = "")), silent = TRUE)
                v.graph.all(out, plot = FALSE, draw = "pos-neg", 
                  main = paste("geo coords: x=", stcrds[icord, 
                    1], " y=", stcrds[icord, 2], sep = ""))
            }
        }
    }
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
