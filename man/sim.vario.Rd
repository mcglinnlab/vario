\name{sim.vario}
\alias{sim.vario}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.vario(PREFIX, DESTDIR = paste(getwd(), "/results_main", sep = ""), REPS = 10, REP.ID = NULL, M = 64, ENV = "fractal", K = 250, S = 10, NICH.WD = c(0.1, 0.5), DISP.WD = c(0.1, 0.5), DEATH = 0.5, BIRTH = 0.505, IMMIGRATE = 0.005, FEC = 1, FRACTAL = c(2.01, 2.25, 2.5, 2.75, 2.99), CYCLES = 1, TIME = 10000, MIG = 1, NPERM = 100, ALL = TRUE, NPAR = 5, RPARGS = c(TRUE, M/3, 1000, 1e+06, 0.01, 1), DESTDIRRAW = paste(getwd(), "/results_raw", sep = ""), SIMINPUT = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PREFIX}{
%%     ~~Describe \code{PREFIX} here~~
}
  \item{DESTDIR}{
%%     ~~Describe \code{DESTDIR} here~~
}
  \item{REPS}{
%%     ~~Describe \code{REPS} here~~
}
  \item{REP.ID}{
%%     ~~Describe \code{REP.ID} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{ENV}{
%%     ~~Describe \code{ENV} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{NICH.WD}{
%%     ~~Describe \code{NICH.WD} here~~
}
  \item{DISP.WD}{
%%     ~~Describe \code{DISP.WD} here~~
}
  \item{DEATH}{
%%     ~~Describe \code{DEATH} here~~
}
  \item{BIRTH}{
%%     ~~Describe \code{BIRTH} here~~
}
  \item{IMMIGRATE}{
%%     ~~Describe \code{IMMIGRATE} here~~
}
  \item{FEC}{
%%     ~~Describe \code{FEC} here~~
}
  \item{FRACTAL}{
%%     ~~Describe \code{FRACTAL} here~~
}
  \item{CYCLES}{
%%     ~~Describe \code{CYCLES} here~~
}
  \item{TIME}{
%%     ~~Describe \code{TIME} here~~
}
  \item{MIG}{
%%     ~~Describe \code{MIG} here~~
}
  \item{NPERM}{
%%     ~~Describe \code{NPERM} here~~
}
  \item{ALL}{
%%     ~~Describe \code{ALL} here~~
}
  \item{NPAR}{
%%     ~~Describe \code{NPAR} here~~
}
  \item{RPARGS}{
%%     ~~Describe \code{RPARGS} here~~
}
  \item{DESTDIRRAW}{
%%     ~~Describe \code{DESTDIRRAW} here~~
}
  \item{SIMINPUT}{
%%     ~~Describe \code{SIMINPUT} here~~
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
function (PREFIX, DESTDIR = paste(getwd(), "/results_main", sep = ""), 
    REPS = 10, REP.ID = NULL, M = 64, ENV = "fractal", K = 250, 
    S = 10, NICH.WD = c(0.1, 0.5), DISP.WD = c(0.1, 0.5), DEATH = 0.5, 
    BIRTH = 0.505, IMMIGRATE = 0.005, FEC = 1, FRACTAL = c(2.01, 
        2.25, 2.5, 2.75, 2.99), CYCLES = 1, TIME = 10000, MIG = 1, 
    NPERM = 100, ALL = TRUE, NPAR = 5, RPARGS = c(TRUE, M/3, 
        1000, 1e+06, 0.01, 1), DESTDIRRAW = paste(getwd(), "/results_raw", 
        sep = ""), SIMINPUT = NULL) 
{
    if (!is.null(SIMINPUT)) {
        M = SIMINPUT$M
        S = SIMINPUT$S
        K = SIMINPUT$K
        no.bbs = SIMINPUT$no.bbs
        ENV = SIMINPUT$envir
        grad.len = dist(range(ENV, na.rm = TRUE))[1]
        FRACTAL = round(D.vario.est(ENV), 2)
    }
    else grad.len = M
    coord <- cbind(rep(1:M, each = M), rep(1:M, times = M))
    Dist <- dist(coord)
    H <- round(Dist)
    hmax <- round(max(Dist)/2)
    H[H > hmax] <- NA
    H <- as.vector(H)
    Dist <- sapply(split(Dist, H), mean)
    OUT <- NULL
    OBS <- array(NA, dim = c(REPS, length(Dist), 4))
    MEAN <- array(NA, dim = c(REPS, length(Dist), 4))
    SD <- array(NA, dim = c(REPS, length(Dist), 4))
    Z <- array(NA, dim = c(REPS, length(Dist), 4))
    P.val <- array(NA, dim = c(REPS, length(Dist), 4))
    Reject.025 <- array(NA, dim = c(REPS, length(Dist), 4))
    Reject.975 <- array(NA, dim = c(REPS, length(Dist), 4))
    if (is.null(REP.ID)) 
        REP.ID <- 1:REPS
    for (D in 1:length(FRACTAL)) {
        for (s in 1:length(NICH.WD)) {
            for (u in 1:length(DISP.WD)) {
                for (m in 1:REPS) {
                  if (is.character(ENV)) {
                    HAB <- gen.env(M = M, env = ENV, D = FRACTAL[D])
                    OPTIM <- runif(S, 0.1, 0.9) * M
                  }
                  else {
                    HAB <- ifelse(is.na(ENV), 0, ENV)
                    OPTIM <- runif(S, range(ENV, na.rm = TRUE))
                  }
                  if (is.null(SIMINPUT)) 
                    D.est <- D.vario.est(HAB)
                  else D.est <- FRACTAL
                  SIM <- neut.sim.uni(M = M, K = K, S = S, s = NICH.WD[s], 
                    u = DISP.WD[u], b = BIRTH, d = DEATH, m = IMMIGRATE, 
                    habitat = HAB, fitness = OPTIM, fec = FEC, 
                    time = TIME, cycles = CYCLES, mig = MIG, 
                    grad.len = grad.len)
                  mat <- apply(census(SIM, snap = length(SIM$snaps)), 
                    c(1, 2), as.logical) * 1
                  mat <- mat[, apply(mat, 2, sum) > 0]
                  mat[as.vector(no.bbs), ] <- -999
                  V.rand <- vario(mat, SIM$coords, pos.neg = FALSE)
                  VP.rand <- null.perms(mat, V.rand, NPERM, coords = SIM$coords, 
                    meth = "random", npar = NPAR)
                  V.spat <- vario(mat, SIM$coords, pos.neg = TRUE)
                  VP.spat <- null.perms(mat, V.spat, NPERM, coords = SIM$coords, 
                    meth = "randpat", RPargs = RPARGS, npar = NPAR)
                  VP <- array(NA, dim = dim(VP.spat$vario) + 
                    c(0, 1, 0))
                  VP[, 1, ] <- VP.rand$vario[, 1, ]
                  VP[, -1, ] <- VP.spat$vario
                  VP[, 2, ] <- VP.spat$vario[, 2, ] + VP.spat$vario[, 
                    3, ]
                  if (!is.null(DESTDIRRAW)) {
                    OBJS <- paste(PREFIX, "M", M, "D", FRACTAL[D], 
                      "K", max(K), "S", S, "s", NICH.WD[s], "u", 
                      DISP.WD[u], "rep", REP.ID[m], ".Rdata", 
                      sep = "")
                    save(SIM, VP.rand, VP.spat, file = paste(DESTDIRRAW, 
                      "/", OBJS, sep = ""))
                  }
                  OBS[m, , ] <- VP[, , 1]
                  MEAN[m, , ] <- apply(VP, 1:2, mean)
                  SD[m, , ] <- apply(VP, 1:2, sd)
                  Z[m, , ] <- (OBS[m, , ] - MEAN[m, , ])/SD[m, 
                    , ]
                  for (k in 1:length(Dist)) {
                    for (i in 1:4) {
                      if (VP[k, i, 1] > MEAN[m, k, i]) 
                        P.val[m, k, i] <- sum(VP[k, i, 1] <= 
                          VP[k, i, ])/(NPERM + 1)
                      else P.val[m, k, i] <- sum(VP[k, i, 1] >= 
                        VP[k, i, ])/(NPERM + 1)
                    }
                  }
                  Reject.025[m, , ] <- apply(VP, 1:2, quantile, 
                    0.025)
                  Reject.975[m, , ] <- apply(VP, 1:2, quantile, 
                    0.975)
                  if (is.null(OUT)) {
                    len <- length(Dist) * 4
                    len2 <- length(Dist)
                    OUT <- data.frame(M = rep(M, len), K = rep(max(K), 
                      len), S = rep(S, len), S.rel = rep(ncol(mat), 
                      len), D = rep(FRACTAL[D], len), D.est = rep(D.est, 
                      len), s = rep(NICH.WD[s], len), s.rel = rep(round(NICH.WD[s] * 
                      grad.len, 2), len), u = rep(DISP.WD[u], 
                      len), u.rel = rep(round(DISP.WD[u] * M, 
                      2), len), rep = rep(REP.ID[m], len), mod = rep(c("rand", 
                      rep("spat", 3)), each = len2), type = rep(c("exp", 
                      "obs", "pos", "neg"), each = len2), dist = rep(Dist, 
                      4), var = as.vector(OBS[m, , ]), mean = as.vector(MEAN[m, 
                      , ]), z = as.vector(Z[m, , ]), p = as.vector(P.val[m, 
                      , ]), avg.rand = as.vector(MEAN[m, , ]), 
                      CI.low = as.vector(Reject.025[m, , ]), 
                      CI.up = as.vector(Reject.975[m, , ]))
                  }
                  else {
                    OUT <- rbind(OUT, data.frame(M = rep(M, len), 
                      K = rep(max(K), len), S = rep(S, len), 
                      S.rel = rep(ncol(mat), len), D = rep(FRACTAL[D], 
                        len), D.est = rep(D.est, len), s = rep(NICH.WD[s], 
                        len), s.rel = rep(NICH.WD[s] * grad.len, 
                        len), u = rep(DISP.WD[u], len), u.rel = rep(DISP.WD[u] * 
                        M, len), rep = rep(REP.ID[m], len), mod = rep(c("rand", 
                        rep("spat", 3)), each = len2), type = rep(c("exp", 
                        "obs", "pos", "neg"), each = len2), dist = rep(Dist, 
                        4), var = as.vector(OBS[m, , ]), mean = as.vector(MEAN[m, 
                        , ]), z = as.vector(Z[m, , ]), p = as.vector(P.val[m, 
                        , ]), avg.rand = as.vector(MEAN[m, , 
                        ]), CI.low = as.vector(Reject.025[m, 
                        , ]), CI.up = as.vector(Reject.975[m, 
                        , ])))
                  }
                }
                if (length(REP.ID) == 1) {
                  OBJS.sum1 <- paste(PREFIX, "M", M, "D", FRACTAL[D], 
                    "K", max(K), "S", S, "s", NICH.WD[s], "u", 
                    DISP.WD[u], "rep", REP.ID, ".Rdata", sep = "")
                  OBJS.sum2 <- paste(PREFIX, "M", M, "D", FRACTAL[D], 
                    "K", max(K), "S", S, "s", NICH.WD[s], "u", 
                    DISP.WD[u], "rep", REP.ID, ".csv", sep = "")
                }
                else {
                  OBJS.sum1 <- paste(PREFIX, "M", M, "D", FRACTAL[D], 
                    "K", max(K), "S", S, "s", NICH.WD[s], "u", 
                    DISP.WD[u], ".Rdata", sep = "")
                  OBJS.sum2 <- paste(PREFIX, "M", M, "D", FRACTAL[D], 
                    "K", max(K), "S", S, "s", NICH.WD[s], "u", 
                    DISP.WD[u], ".csv", sep = "")
                }
                save(OBS, MEAN, Z, SD, P.val, Reject.025, Reject.975, 
                  file = paste(DESTDIR, "/", OBJS.sum1, sep = ""))
                write.csv(OUT, file = paste(DESTDIR, "/", OBJS.sum2, 
                  sep = ""))
            }
        }
    }
    if (length(REP.ID) > 1) {
        OBJS.final <- paste(PREFIX, "M", M, "D", FRACTAL[D], 
            "K", max(K), "S", S, ".csv", sep = "")
        write.csv(OUT, file = paste(DESTDIR, "/", OBJS.final, 
            sep = ""))
    }
    OUT
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
