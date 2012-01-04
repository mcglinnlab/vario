\name{v.graph.all2}
\alias{v.graph.all2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
v.graph.all2(vrand = NULL, vspat = NULL, obs.var = FALSE, flip.neg = FALSE, ylim1 = NA, ylim2 = NA, xlim = NA, ylab = "variance", xlab = "lag", cls = NA, lwd = 1, cex.axis = 1, box = TRUE, plot.new = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vrand}{
%%     ~~Describe \code{vrand} here~~
}
  \item{vspat}{
%%     ~~Describe \code{vspat} here~~
}
  \item{obs.var}{
%%     ~~Describe \code{obs.var} here~~
}
  \item{flip.neg}{
%%     ~~Describe \code{flip.neg} here~~
}
  \item{ylim1}{
%%     ~~Describe \code{ylim1} here~~
}
  \item{ylim2}{
%%     ~~Describe \code{ylim2} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{cls}{
%%     ~~Describe \code{cls} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{cex.axis}{
%%     ~~Describe \code{cex.axis} here~~
}
  \item{box}{
%%     ~~Describe \code{box} here~~
}
  \item{plot.new}{
%%     ~~Describe \code{plot.new} here~~
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
function (vrand = NULL, vspat = NULL, obs.var = FALSE, flip.neg = FALSE, 
    ylim1 = NA, ylim2 = NA, xlim = NA, ylab = "variance", xlab = "lag", 
    cls = NA, lwd = 1, cex.axis = 1, box = TRUE, plot.new = TRUE) 
{
    rflag <- !is.null(vrand)
    sflag <- !is.null(vspat)
    dr <- NULL
    ds <- NULL
    if (rflag) {
        dr <- vrand$vdists
        vr <- vrand$vario[, 1, ]
    }
    if (sflag) {
        ds <- vspat$vdists
        vs <- vspat$vario
        if (flip.neg) 
            vs[, 3, ] <- vs[, 3, ] * -1
    }
    if (is.na(xlim[1])) 
        xlim <- range(dr, ds)
    if (is.na(cls)) {
        cls <- c(rgb(190, 174, 212, alpha = 255 * 0.5, maxColorValue = 255), 
            "#99CCFF", rgb(127, 201, 127, alpha = 255 * 0.5, 
                maxColorValue = 255))
    }
    if (rflag & sflag & plot.new) {
        if (obs.var) 
            par(mfrow = c(1, 3))
        else par(mfrow = c(1, 2))
    }
    if (rflag) {
        q.rand <- apply(vr, 1, quantile, c(0.025, 0.975))
        if (is.na(ylim1[1])) 
            ylim1 = range(q.rand, vr[, 1])
        plot(dr, vr[, 1], type = "n", ylab = ylab, xlab = xlab, 
            ylim = ylim1, xlim = xlim, main = "Within-species Agg.", 
            axes = F)
        axis(side = 1, cex.axis = cex.axis)
        axis(side = 2, cex.axis = cex.axis)
        polygon(c(dr, dr[length(dr):1]), c(q.rand[1, ], q.rand[2, 
            length(dr):1]), border = NA, col = cls[1])
        lines(dr, vr[, 1], col = 1, lwd = lwd)
        if (box) 
            box()
    }
    if (sflag) {
        if (obs.var) {
            if (flip.neg) 
                ovar <- vs[, 2, ] - vs[, 3, ]
            else ovar <- vs[, 2, ] + vs[, 3, ]
            q.obs <- apply(ovar, 1, quantile, c(0.025, 0.975))
            if (is.na(ylim1[1])) 
                ylim1 = range(q.obs, ovar)
            plot(ds, ovar[, 1], type = "n", ylab = ylab, xlab = xlab, 
                ylim = ylim1, xlim = xlim, main = "Total Between-species Agg.", 
                axes = F)
            axis(side = 1, cex.axis = cex.axis)
            axis(side = 2, cex.axis = cex.axis)
            polygon(c(ds, ds[length(ds):1]), c(q.obs[1, ], q.obs[2, 
                length(ds):1]), border = NA, col = cls[1])
            lines(ds, ovar[, 1], col = 1, lwd = lwd)
            if (box) 
                box()
        }
        q.spat <- apply(vs, 1:2, quantile, c(0.025, 0.975))
        if (is.na(ylim2[1])) 
            ylim2 = range(q.spat, vs[, , 1])
        plot(ds, vs[, 2, 1], type = "n", ylab = ylab, xlab = xlab, 
            ylim = ylim2, xlim = xlim, main = "Pos/Neg Between-species Agg.", 
            axes = F)
        axis(side = 1, cex.axis = cex.axis)
        axis(side = 2, cex.axis = cex.axis)
        polygon(c(ds, ds[length(ds):1]), c(q.spat[1, , 2], q.spat[2, 
            length(ds):1, 2]), border = NA, col = cls[2])
        polygon(c(ds, ds[length(ds):1]), c(q.spat[1, , 3], q.spat[2, 
            length(ds):1, 3]), border = NA, col = cls[3])
        lines(ds, vs[, 2, 1], col = "dodgerblue", lwd = lwd)
        lines(ds, vs[, 3, 1], col = "green3", lwd = lwd)
        if (box) 
            box()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
