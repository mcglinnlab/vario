\name{Empir.pdf}
\alias{Empir.pdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Empir.pdf(stcrds, n, nsamp = 2, both = TRUE)
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
function (stcrds, n, nsamp = 2, both = TRUE) 
{
    if (both) 
        pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
            n, "/empir_n", n, ".pdf", sep = ""), width = 7 * 
            3, height = 7 * 2)
    else pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
        n, "/empir_n", n, ".pdf", sep = ""), width = 7 * 2.5, 
        height = 7)
    icount <- 0
    for (icord in 1:length(stcrds[, 1])) {
        out <- NULL
        try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n", 
            n, "/empir.n", n, ".x", stcrds[icord, 1], ".y", stcrds[icord, 
                2], ".Rdata", sep = "")), silent = TRUE)
        if (!is.null(out)) {
            if (icount == 0) {
                quants <- apply(out$vario, 1:2, quantile, c(0.025, 
                  0.975))
                ylims <- range(list(quants, out$vario[, , 1]))
            }
            if (both) 
                par(mfrow = c(2, 3))
            else par(mfrow = c(1, 3))
            brks <- seq(nsamp, 85, 5)
            cls <- terrain.colors(length(brks) - 1)
            image(Bn100$x, Bn100$y, Bn100$n, breaks = brks, col = cls, 
                ylim = c(-100, 3500), xlim = c(-100, 5500), main = paste("min samp", 
                  nsamp), xlab = "", ylab = "")
            points(stcrds[icord, 1], stcrds[icord, 2], col = "purple", 
                pch = 19, cex = 3)
            polygon(c(rep(stcrds[icord, 1], 2), rep(stcrds[icord, 
                1] + n * len1, 2)), c(stcrds[icord, 2], rep(stcrds[icord, 
                2] + n * len1, 2), stcrds[icord, 2]))
            lines(wrld.pts, col = "red", lwd = 2)
            v.graph.all(out, plot = FALSE, ylims = ylims, main = paste("geo coords: x=", 
                stcrds[icord, 1], " y=", stcrds[icord, 2], sep = ""))
            if (both) {
                plot(1:10, 1:10, type = "n", axes = F, xlab = "", 
                  ylab = "", frame.plot = F)
                diff <- array(rep(out$vario[, , 1], dim(out$vario)[3]), 
                  dim = dim(out$vario)) - out$vario
                DIFF <- apply(diff, 1:2, mean)
                sds <- apply(diff, 1:2, sd)
                DELTA <- DIFF/sds
                if (icount == 0) 
                  ylims2 <- c(0, max(DELTA))
                plot(out$vdist, DELTA[, 1], ylim = ylims2, type = "o", 
                  main = paste("geo coords: x=", stcrds[icord, 
                    1], " y=", stcrds[icord, 2], sep = ""))
                plot(out$vdist, DELTA[, 2], ylim = ylims2, type = "o", 
                  col = "red", main = paste("geo coords: x=", 
                    stcrds[icord, 1], " y=", stcrds[icord, 2], 
                    sep = ""))
                lines(out$vdist, DELTA[, 3], type = "o", col = "blue")
                legend("bottomleft", c("pos", "neg"), col = c("red", 
                  "blue"), pch = 19, lty = 1, lwd = 2, bty = "n")
            }
            icount <- 1
        }
    }
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
