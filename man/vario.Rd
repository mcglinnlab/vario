\name{vario}
\alias{vario}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vario(x, coord, grain = 1, breaks = NA, hmax = NA, round.up = FALSE, pos.neg = FALSE, binary = TRUE, snap = NA, median = FALSE, direction = "omnidirectional", tolerance = pi/8, unit.angle = c("radians", "degrees"), distance.metric = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{grain}{
%%     ~~Describe \code{grain} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
}
  \item{hmax}{
%%     ~~Describe \code{hmax} here~~
}
  \item{round.up}{
%%     ~~Describe \code{round.up} here~~
}
  \item{pos.neg}{
%%     ~~Describe \code{pos.neg} here~~
}
  \item{binary}{
%%     ~~Describe \code{binary} here~~
}
  \item{snap}{
%%     ~~Describe \code{snap} here~~
}
  \item{median}{
%%     ~~Describe \code{median} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
}
  \item{tolerance}{
%%     ~~Describe \code{tolerance} here~~
}
  \item{unit.angle}{
%%     ~~Describe \code{unit.angle} here~~
}
  \item{distance.metric}{
%%     ~~Describe \code{distance.metric} here~~
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
function (x, coord, grain = 1, breaks = NA, hmax = NA, round.up = FALSE, 
    pos.neg = FALSE, binary = TRUE, snap = NA, median = FALSE, 
    direction = "omnidirectional", tolerance = pi/8, unit.angle = c("radians", 
        "degrees"), distance.metric = "euclidean") 
{
    if (distance.metric != "euclidean") {
        if (pos.neg) 
            stop("cannot commpute pos-neg covariance using a turnover metric")
    }
    else require(vegan)
    unit.angle = match.arg(unit.angle)
    if (mode(direction) == "numeric") {
        if (length(direction) > 1) 
            stop("only one direction is allowed")
        if (length(tolerance) > 1) 
            stop("only one tolerance value is allowed")
        if (unit.angle == "degrees") {
            ang.deg = direction
            ang.rad = (ang.deg * pi)/180
            tol.deg = tolerance
            tol.rad = (tol.deg * pi)/180
        }
        else {
            ang.rad = direction
            ang.deg = (ang.rad * 180)/pi
            tol.rad = tolerance
            tol.deg = (tol.rad * 180)/pi
        }
        if (ang.rad > pi | ang.rad < 0) 
            stop("direction must be an angle in the interval [0,pi[ radians")
        if (tol.rad > pi/2 | tol.rad < 0) 
            stop("tolerance must be an angle in the interval [0,pi/2] radians")
        if (tol.deg >= 90) {
            direction = "omnidirectional"
        }
    }
    if (class(x) == "sim") {
        coord = x$coords
        if (is.na(snap)) 
            snap = length(sim$snaps)
    }
    else x = ifelse(x == -999, NA, x)
    Dist = dist(coord)
    if (is.na(breaks[1])) {
        if (round.up) 
            H = ceiling(Dist/grain) * grain
        else H = round(Dist/grain) * grain
    }
    else {
        H = Dist
        for (i in 1:(length(breaks) - 1)) H[H >= breaks[i] & 
            H < breaks[i + 1]] = breaks[i]
    }
    if (is.na(hmax)) 
        hmax = round((max(Dist)/2)/grain) * grain
    vobject = list()
    if (class(x) == "sim") {
        if (binary) 
            resp = apply(census(x, snap = snap), c(1, 2), as.logical) * 
                1
        else resp = census(x, snap = snap)
        vobject$parms = data.frame(grain, hmax, S = ncol(resp), 
            N = nrow(resp), pos.neg, median, niche.wid.rel = x$p$s.rel, 
            disp.wid.rel = x$p$u.rel, direction, tolerance, unit.angle, 
            distance.metric)
    }
    else {
        resp = as.matrix(x)
        vobject$parms = data.frame(grain, hmax, S = ncol(resp), 
            N = nrow(resp), pos.neg, median, direction, tolerance, 
            unit.angle, distance.metric)
    }
    H[H > hmax] = NA
    H = as.vector(H)
    if (direction != "omnidirectional") {
        u.ang = .C("tgangle", as.double(as.vector(coord[, 1])), 
            as.double(as.vector(coord[, 2])), as.integer(dim(coord)[1]), 
            res = as.double(rep(0, length(as.vector(Dist)))))$res
        if (any(is.na(u.ang))) 
            stop("NA returned in angle calculations maybe due to co-located data")
        u.ang = atan(u.ang)
        u.ang[u.ang < 0] = u.ang[u.ang < 0] + pi
        ang.lower = ang.rad - tol.rad
        ang.upper = ang.rad + tol.rad
        if (ang.lower >= 0 & ang.upper < pi) 
            ang.ind = (!is.na(u.ang) & ((u.ang >= ang.lower) & 
                (u.ang <= ang.upper)))
        if (ang.lower < 0) 
            ang.ind = (!is.na(u.ang) & ((u.ang < ang.upper) | 
                (u.ang > (pi + ang.lower))))
        if (ang.upper >= pi) 
            ang.ind = (!is.na(u.ang) & ((u.ang > ang.lower) | 
                (u.ang < (ang.upper - pi))))
        Dist[!ang.ind] = NA
        H[!ang.ind] = NA
    }
    Dist = sapply(split(Dist, H), mean, na.rm = TRUE)
    vobject$vario = data.frame(H = names(table(H)), Dist = Dist, 
        n = as.numeric(table(H)))
    if (distance.metric == "euclidean") 
        exp.split = split(dist(resp)^2 * 0.5, H)
    else exp.split = split(vegdist(resp, method = distance.metric), 
        H)
    exp.gamma = sapply(exp.split, mean, na.rm = TRUE)
    vobject$vario = cbind(vobject$vario, exp.var = exp.gamma)
    if (median) 
        exp.med = sapply(exp.split, median, na.rm = TRUE)
    if (ncol(resp) > 1) {
        rich = apply(resp, 1, sum)
        obs.gamma = sapply(split(dist(rich)^2 * 0.5, H), mean, 
            na.rm = TRUE)
        if (pos.neg) {
            cov.mat = getCovFractions(resp)
            pos.split = split(cov.mat$pos, H)
            neg.split = split(cov.mat$neg, H)
            pos = sapply(pos.split, mean)
            neg = sapply(neg.split, mean)
            if (median) {
                pos.med = sapply(pos.split, median)
                neg.med = sapply(neg.split, median)
                vobject$vario = cbind(vobject$vario, obs.var = obs.gamma, 
                  ratio = obs.gamma/exp.gamma, pos = pos, neg = neg, 
                  exp.med = exp.med, pos.med = pos.med, neg.med = neg.med)
            }
            else vobject$vario = cbind(vobject$vario, obs.var = obs.gamma, 
                ratio = obs.gamma/exp.gamma, pos = pos, neg = neg)
        }
        else {
            if (median) {
                obs.gamma.med = sapply(split(dist(rich)^2 * 0.5, 
                  H), median, na.rm = TRUE)
                vobject$vario = cbind(vobject$vario, obs.var = obs.gamma, 
                  ratio = obs.gamma/exp.gamma, exp.med = exp.med, 
                  obs.med = obs.gamma.med)
            }
            else vobject$vario = cbind(vobject$vario, obs.var = obs.gamma, 
                ratio = obs.gamma/exp.gamma)
        }
    }
    vobject$p = apply(resp, 2, sum, na.rm = TRUE)/nrow(resp)
    vobject$perm = FALSE
    return(vobject)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
