\name{v.graph}
\alias{v.graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
v.graph(vobject, optim = NA, exp.only = FALSE, flip.neg = FALSE, ylim = NULL, xlab = NULL, ylab = NULL, cls = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vobject}{
%%     ~~Describe \code{vobject} here~~
}
  \item{optim}{
%%     ~~Describe \code{optim} here~~
}
  \item{exp.only}{
%%     ~~Describe \code{exp.only} here~~
}
  \item{flip.neg}{
%%     ~~Describe \code{flip.neg} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{cls}{
%%     ~~Describe \code{cls} here~~
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
function (vobject, optim = NA, exp.only = FALSE, flip.neg = FALSE, 
    ylim = NULL, xlab = NULL, ylab = NULL, cls = NULL) 
{
    n <- sqrt(vobject$parms$N)
    N <- n^2
    S <- vobject$parms$S
    coords <- cbind(rep(1:n, each = n), rep(1:n, times = n))
    v <- vobject$vario
    p <- vobject$p
    if (is.null(cls)) {
        cls <- c(rgb(127, 201, 127, alpha = 255 * 0.5, maxColorValue = 255), 
            rgb(190, 174, 212, alpha = 255 * 0.5, maxColorValue = 255), 
            rgb(253, 192, 134, alpha = 255 * 0.5, maxColorValue = 255))
    }
    if (is.null(xlab)) 
        xlab <- "Spatial Lag"
    if (is.null(ylab)) 
        ylab <- "Variance"
    if (!is.na(optim[1])) {
        y <- 1:S
        plot(sort(optim), y, axes = F, frame.plot = F, xlab = "grad", 
            ylab = "", xlim = c(0, max(optim)))
        axis(side = 1)
        axis(side = 2, at = c(1, S))
        arrows(sort(optim) - vobject$parms$niche.wid.rel/2, y, 
            sort(optim) + vobject$parms$niche.wid.rel/2, y, angle = 90, 
            len = 0.01, code = 3)
    }
    if (vobject$perm) {
        d <- vobject$vdists
        quants <- apply(v, c(1, 2), quantile, c(0.025, 0.975))
        if (is.null(ylim)) 
            ylim = range(list(quants, v[, , 1]))
        if (vobject$parms$pos.neg) {
            plot(d, v[, 1, 1], col = 1, pch = 19, type = "o", 
                ylab = ylab, xlab = xlab, ylim = ylim)
            points(d, v[, 1, 1] - v[, 2, 1] - v[, 3, 1], col = "green3", 
                type = "o", pch = 19)
            lines(d, v[, 2, 1], col = "purple", lwd = 2)
            lines(d, v[, 3, 1], col = "orange", lwd = 2)
        }
        else {
            plot(d, v[, 1, 1], pch = 19, type = "o", col = "green3", 
                ylab = ylab, xlab = xlab, ylim = ylim)
            lines(d, v[, 2, 1], col = 1, lwd = 2, type = "o", 
                pch = 19)
        }
        lines(d, rep(sum(p * (1 - p)), length(d)), col = "blue", 
            lwd = 2)
        for (i in 1:dim(quants)[3]) polygon(c(d, d[length(d):1]), 
            c(quants[1, , i], quants[2, length(d):1, i]), col = cls[i], 
            border = NA)
    }
    else {
        if (is.null(ylim)) {
            if (vobject$parms$pos.neg) {
                if (flip.neg) 
                  v[, 8] <- v[, 8] * -1
                ylim = range(list(v[, c(4:5, 7:8)], sum(p * (1 - 
                  p))))
            }
            else {
                if (exp.only) 
                  ylim = range(list(v[, 4], sum(p * (1 - p))))
                else ylim = range(list(v[, 4:5], sum(p * (1 - 
                  p))))
            }
        }
        plot(v$Dist, v$obs.var, col = 1, pch = 19, type = "n", 
            ylab = ylab, xlab = xlab, ylim = ylim)
        points(v$Dist, v$exp.var, col = "green3", type = "o", 
            pch = 19)
        if (!exp.only) 
            points(v$Dist, v$obs.var, col = 1, pch = 19, type = "o")
        if (vobject$parms$pos.neg) {
            lines(v$Dist, v$pos, col = "purple", lwd = 2)
            lines(v$Dist, v$neg, col = "orange", lwd = 2)
        }
        lines(v$Dist, rep(sum(p * (1 - p)), length(v$Dist)), 
            col = "blue", lwd = 2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
