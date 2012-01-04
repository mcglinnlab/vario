\name{v.graph.all}
\alias{v.graph.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
v.graph.all(vrand = NULL, vspat = NULL, obs.var = FALSE, ylims = NA, xlims = NA, ylab = "variance", xlab = "lag", cls = NA, lwd = 1, plot.new = TRUE)
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
  \item{ylims}{
%%     ~~Describe \code{ylims} here~~
}
  \item{xlims}{
%%     ~~Describe \code{xlims} here~~
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
function (vrand = NULL, vspat = NULL, obs.var = FALSE, ylims = NA, 
    xlims = NA, ylab = "variance", xlab = "lag", cls = NA, lwd = 1, 
    plot.new = TRUE) 
{
    rflag <- !is.null(vrand)
    sflag <- !is.null(vspat)
    if (rflag) {
        dr <- vrand$vdists
        vr <- vrand$vario[, 1, ]
    }
    if (sflag) {
        ds <- vspat$vdists
        vs <- vspat$vario
    }
    if (is.na(xlims[1])) 
        xlims <- range(dr, ds)
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
        ylims = range(q.rand, vr[, 1])
        plot(dr, vr[, 1], type = "n", ylab = ylab, xlab = xlab, 
            ylim = ylims, xlim = xlims, main = "Within-species Agg.")
        polygon(c(dr, dr[length(dr):1]), c(q.rand[1, ], q.rand[2, 
            length(dr):1]), border = NA, col = cls[1])
        lines(dr, vr[, 1], col = 1, lwd = lwd)
    }
    if (sflag) {
        if (obs.var) {
            ovar <- vs[, 2, ] + vs[, 3, ]
            q.obs <- apply(ovar, 1, quantile, c(0.025, 0.975))
            ylims = range(q.obs, ovar)
            plot(ds, ovar[, 1], type = "n", ylab = ylab, xlab = xlab, 
                ylim = ylims, xlim = xlims, main = "Total Between-species Agg.")
            polygon(c(ds, ds[length(ds):1]), c(q.obs[1, ], q.obs[2, 
                length(ds):1]), border = NA, col = cls[1])
            lines(ds, ovar[, 1], col = 1, lwd = lwd)
        }
        q.spat <- apply(vs, 1:2, quantile, c(0.025, 0.975))
        ylims = range(q.spat, vs[, , 1])
        plot(ds, vs[, 2, 1], type = "n", ylab = ylab, xlab = xlab, 
            ylim = ylims, xlim = xlims, main = "Pos/Neg Between-species Agg.")
        polygon(c(ds, ds[length(ds):1]), c(q.spat[1, , 2], q.spat[2, 
            length(ds):1, 2]), border = NA, col = cls[2])
        polygon(c(ds, ds[length(ds):1]), c(q.spat[1, , 3], q.spat[2, 
            length(ds):1, 3]), border = NA, col = cls[3])
        lines(ds, vs[, 2, 1], col = "dodgerblue", lwd = lwd)
        lines(ds, vs[, 3, 1], col = "green3", lwd = lwd)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
