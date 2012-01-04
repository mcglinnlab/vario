\name{null.perms}
\alias{null.perms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
null.perms(x, vobject, nperm, coords = NULL, meth = "both", sp = TRUE, all = FALSE, snap = NULL, npar = 1, linux = FALSE, RPargs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{vobject}{
%%     ~~Describe \code{vobject} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{sp}{
%%     ~~Describe \code{sp} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
}
  \item{snap}{
%%     ~~Describe \code{snap} here~~
}
  \item{npar}{
%%     ~~Describe \code{npar} here~~
}
  \item{linux}{
%%     ~~Describe \code{linux} here~~
}
  \item{RPargs}{
%%     ~~Describe \code{RPargs} here~~
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
function (x, vobject, nperm, coords = NULL, meth = "both", sp = TRUE, 
    all = FALSE, snap = NULL, npar = 1, linux = FALSE, RPargs = FALSE) 
{
    dists <- vobject$vario$Dist
    grain = vobject$parms$grain
    hmax = vobject$parms$hmax
    pos.neg = vobject$parms$pos.neg
    median = vobject$parms$median
    if (class(vobject$parms$direction) == "factor") 
        direction = as.character(vobject$parms$direction)
    else direction = as.numeric(vobject$parms$direction)
    tolerance = vobject$parms$tolerance
    unit.angle = as.character(vobject$parms$unit.angle)
    hmax <- vobject$parms$hmax
    if (class(x) == "sim") {
        coords <- x$coords
        if (is.null(snap)) 
            snap <- length(x$snaps)
    }
    else {
        if (is.null(coords)) 
            stop("need to supply spatial coordinates if not a simulation product")
    }
    r.vals <- list()
    r.vals$parms <- vobject$parms
    r.vals$p <- vobject$p
    if (median & !pos.neg) 
        stop("if computing medians must also compute pos.neg fractions, set pos.neg=TRUE")
    if (pos.neg) {
        if (all) {
            if (median) {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  6, 3, nperm + 1))
                r.vals$vario[, , , 1] <- as.matrix(vobject$vario[, 
                  c(5, 7:11)])
            }
            else {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  3, 3, nperm + 1))
                r.vals$vario[, , , 1] <- as.matrix(vobject$vario[, 
                  c(5, 7:8)])
            }
        }
        else {
            if (median) {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  6, nperm + 1))
                r.vals$vario[, , 1] <- as.matrix(vobject$vario[, 
                  c(5, 7:11)])
            }
            else {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  3, nperm + 1))
                r.vals$vario[, , 1] <- as.matrix(vobject$vario[, 
                  c(5, 7:8)])
            }
        }
    }
    else {
        if (all) {
            if (median) {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  4, 3, nperm + 1))
                r.vals$vario[, , , 1] <- as.matrix(vobject$v[, 
                  c(4:5, 7:8)])
            }
            else {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  2, 3, nperm + 1))
                r.vals$vario[, , , 1] <- as.matrix(vobject$v[, 
                  4:5])
            }
        }
        else {
            if (median) {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  4, nperm + 1))
                r.vals$vario[, , 1] <- as.matrix(vobject$v[, 
                  c(4:5, 7:8)])
            }
            else {
                r.vals$vario <- array(0, dim = c(length(dists), 
                  2, nperm + 1))
                r.vals$vario[, , 1] <- as.matrix(vobject$v[, 
                  4:5])
            }
        }
    }
    if (class(x) == "sim") {
        pop <- as.logical(x$snaps[[snap]])
        dim(pop) <- c(x$p$S, x$p$M, x$p$M)
    }
    else {
        pop <- array(x, dim = c(sqrt(nrow(x)), sqrt(nrow(x)), 
            ncol(x)))
        pop <- aperm(pop, c(3, 1, 2))
    }
    if (RPargs[[1]] & npar == 1) {
        r.vals$p.conv1 <- 0
        r.vals$p.conv2 <- 0
    }
    if (npar == 1) {
        for (i in 1:nperm) {
            if (RPargs[[1]]) {
                out <- RandPatPar(psp = pop, nstrata = RPargs[[2]], 
                  mtrials1 = RPargs[[3]], mtrials2 = RPargs[[4]], 
                  alpha = RPargs[[5]], npar = RPargs[[6]])
                S <- dim(pop)[1]
                n2 <- dim(pop)[2] + 2
                r.vals$p.conv1 <- r.vals$p.conv1 + sum(out[2, 
                  ] <= RPargs[[5]])/S/nperm
                r.vals$p.conv2 <- r.vals$p.conv2 + sum(out[4, 
                  ] <= RPargs[[5]], na.rm = TRUE)/S/nperm
                rpop <- array(0, dim = c(S, n2, n2))
                for (k in 1:S) {
                  rpop[k, , ] <- array(out[-(1:5), k], dim = c(n2, 
                    n2))
                }
                rpop <- rpop[, -c(1, n2), -c(1, n2)]
            }
            else {
                rpop <- SpatPerm2D(pop, meth = meth, sp = sp)
                rpop <- FixUnSamp2(pop, rpop)
            }
            rmat <- apply(rpop, 1, as.vector)
            rv <- vario(x = rmat, coord = coords, grain = grain, 
                hmax = hmax, pos.neg = pos.neg, median = median, 
                direction = direction, tolerance = tolerance, 
                unit.angle = unit.angle)$vario
            if (pos.neg) {
                if (all) {
                  if (median) 
                    r.vals$vario[, , , i + 1] <- as.matrix(rv[, 
                      , c(5, 7:11)])
                  else r.vals$vario[, , , i + 1] <- as.matrix(rv[, 
                    , c(5, 7:8)])
                }
                else {
                  if (median) 
                    r.vals$vario[, , i + 1] <- as.matrix(rv[, 
                      c(5, 7:11)])
                  else r.vals$vario[, , i + 1] <- as.matrix(rv[, 
                    c(5, 7:8)])
                }
            }
            else {
                if (all) {
                  if (median) 
                    r.vals$vario[, , , i + 1] <- as.matrix(rv[, 
                      , c(4:5, 7:8)])
                  else r.vals$vario[, , , i + 1] <- as.matrix(rv[, 
                    , 4:5])
                }
                else {
                  if (median) 
                    r.vals$vario[, , i + 1] <- as.matrix(rv[, 
                      c(4:5, 7:8)])
                  else r.vals$vario[, , i + 1] <- as.matrix(rv[, 
                    4:5])
                }
            }
            print(i)
        }
    }
    else {
        require(snowfall)
        sfInit(parallel = TRUE, cpus = npar, type = "SOCK")
        sfClusterSetupRNG()
        sfExport("pop", "vobject", "coords", "meth", "all", "sp", 
            "RPargs", "RandPatPar", "RandPat", "FixUnSamp", "FixUnSamp2", 
            "SpatPerm2D", "SpatPerm2D.str", "vario", "dist.cross.real", 
            "null.gen")
        if (linux) 
            sfClusterEval(dyn.load("danspkg.so"))
        else {
            sfLibrary(danspkg)
        }
        out <- unlist(sfLapply(1:nperm, function(...) null.gen(pop = pop, 
            vobject = vobject, coords = coords, meth = meth, 
            sp = sp, all = all, RPargs = RPargs)))
        sfStop()
        if (pos.neg) {
            if (all) {
                if (median) {
                  dim(out) <- c(length(dists), 6, 3, nperm)
                  r.vals$vario[, , , -1] <- out
                }
                else {
                  dim(out) <- c(length(dists), 3, 3, nperm)
                  r.vals$vario[, , , -1] <- out
                }
            }
            else {
                if (median) {
                  dim(out) <- c(length(dists), 6, nperm)
                  r.vals$vario[, , -1] <- out
                }
                else {
                  dim(out) <- c(length(dists), 3, nperm)
                  r.vals$vario[, , -1] <- out
                }
            }
        }
        else {
            if (all) {
                if (median) {
                  dim(out) <- c(length(dists), 4, 3, nperm)
                  r.vals$vario[, , , -1] <- out
                }
                else {
                  dim(out) <- c(length(dists), 2, 3, nperm)
                  r.vals$vario[, , , -1] <- out
                }
            }
            else {
                if (median) {
                  dim(out) <- c(length(dists), 4, nperm)
                  r.vals$vario[, , -1] <- out
                }
                else {
                  dim(out) <- c(length(dists), 2, nperm)
                  r.vals$vario[, , -1] <- out
                }
            }
        }
    }
    r.vals$perm <- TRUE
    r.vals$vdists <- vobject$vario$Dist
    r.vals
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
