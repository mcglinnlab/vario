\name{SpatPerm2D}
\alias{SpatPerm2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SpatPerm2D(psp, shiftpos = NULL, rotate = NULL, meth = "shift", sp = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psp}{
%%     ~~Describe \code{psp} here~~
}
  \item{shiftpos}{
%%     ~~Describe \code{shiftpos} here~~
}
  \item{rotate}{
%%     ~~Describe \code{rotate} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{sp}{
%%     ~~Describe \code{sp} here~~
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
function (psp, shiftpos = NULL, rotate = NULL, meth = "shift", 
    sp = FALSE) 
{
    n <- dim(psp)[2]
    if (length(dim(psp)) == 3) {
        S <- dim(psp)[1]
        flag <- FALSE
    }
    else {
        S <- 1
        psp <- array(psp, dim = c(S, n, n))
        flag <- TRUE
    }
    Rpsp <- psp
    if (sp) {
        if (meth != "reflect") {
            if (meth != "random") {
                for (j in 1:S) {
                  if (is.null(shiftpos)) {
                    shift.x <- sample(n, size = 1)
                    shift.y <- sample(n, size = 1)
                  }
                  else {
                    shift.x <- shiftpos[1]
                    shift.y <- shiftpos[2]
                  }
                  if (shift.x == 1) 
                    ncoord.x <- 1:n
                  else ncoord.x <- c(shift.x:n, 1:(shift.x - 
                    1))
                  if (shift.y == 1) 
                    ncoord.y <- 1:n
                  else ncoord.y <- c(shift.y:n, 1:(shift.y - 
                    1))
                  if (meth == "shift") {
                    Rpsp[j, , ] <- psp[j, ncoord.x, ncoord.y]
                  }
                  if (meth == "both") {
                    if (is.null(rotate)) 
                      rotate <- sample(4, size = 1)
                    if (rotate == 2) {
                      for (x in 1:n) {
                        for (y in 1:n) {
                          Rpsp[j, (n - y) + 1, x] <- psp[j, x, 
                            y]
                        }
                      }
                    }
                    if (rotate == 3) {
                      for (x in 1:n) {
                        for (y in 1:n) {
                          Rpsp[j, (n - x) + 1, (n - y) + 1] <- psp[j, 
                            x, y]
                        }
                      }
                    }
                    if (rotate == 4) {
                      for (x in 1:n) {
                        for (y in 1:n) {
                          Rpsp[j, y, (n - x) + 1] <- psp[j, x, 
                            y]
                        }
                      }
                    }
                    flips <- sample(2, replace = TRUE)
                    if (flips[1] == 1) {
                      if (flips[2] == 1) 
                        Rpsp[j, n:1, n:1] <- Rpsp[j, ncoord.x, 
                          ncoord.y]
                      else Rpsp[j, n:1, ] <- Rpsp[j, ncoord.x, 
                        ncoord.y]
                    }
                    else {
                      if (flips[2] == 1) 
                        Rpsp[j, , n:1] <- Rpsp[j, ncoord.x, ncoord.y]
                      else Rpsp[j, , ] <- Rpsp[j, ncoord.x, ncoord.y]
                    }
                  }
                }
            }
        }
        if (meth == "reflect") {
            for (j in 1:S) {
                if (is.null(rotate)) 
                  rotate <- sample(4, size = 1)
                if (rotate == 2) {
                  for (x in 1:n) {
                    for (y in 1:n) {
                      Rpsp[j, (n - y) + 1, x] <- psp[j, x, y]
                    }
                  }
                }
                if (rotate == 3) {
                  for (x in 1:n) {
                    for (y in 1:n) {
                      Rpsp[j, (n - x) + 1, (n - y) + 1] <- psp[j, 
                        x, y]
                    }
                  }
                }
                if (rotate == 4) {
                  for (x in 1:n) {
                    for (y in 1:n) {
                      Rpsp[j, y, (n - x) + 1] <- psp[j, x, y]
                    }
                  }
                }
                flips <- sample(2, replace = TRUE)
                if (flips[1] == 1) {
                  if (flips[2] == 1) 
                    Rpsp[j, n:1, n:1] <- Rpsp[j, , ]
                  else Rpsp[j, n:1, ] <- Rpsp[j, , ]
                }
                else {
                  if (flips[2] == 1) 
                    Rpsp[j, , n:1] <- Rpsp[j, , ]
                }
            }
        }
        if (meth == "random") {
            for (j in 1:S) {
                take <- sample(n^2)
                Rpsp[j, , ] <- matrix(psp[j, , ][take], ncol = n, 
                  nrow = n)
            }
        }
    }
    else {
        if (meth == "reflect") {
            stop("Reflecting fixed species co-occurances w/o shifting is not meaningful")
        }
        if (meth != "random") {
            if (is.null(shiftpos)) {
                shift.x <- sample(n, size = 1)
                shift.y <- sample(n, size = 1)
            }
            else {
                shift.x <- shiftpos[1]
                shift.y <- shiftpos[2]
            }
            if (shift.x == 1) 
                ncoord.x <- 1:n
            else ncoord.x <- c(shift.x:n, 1:(shift.x - 1))
            if (shift.y == 1) 
                ncoord.y <- 1:n
            else ncoord.y <- c(shift.y:n, 1:(shift.y - 1))
            if (meth == "shift") {
                Rpsp <- psp[, ncoord.x, ncoord.y]
            }
            if (meth == "both") {
                if (is.null(rotate)) 
                  rotate <- sample(4, size = 1)
                if (rotate == 2) {
                  for (x in 1:n) {
                    for (y in 1:n) {
                      Rpsp[, (n - y) + 1, x] <- psp[, x, y]
                    }
                  }
                }
                if (rotate == 3) {
                  for (x in 1:n) {
                    for (y in 1:n) {
                      Rpsp[, (n - x) + 1, (n - y) + 1] <- psp[, 
                        x, y]
                    }
                  }
                }
                if (rotate == 4) {
                  for (x in 1:n) {
                    for (y in 1:n) {
                      Rpsp[, y, (n - x) + 1] <- psp[, x, y]
                    }
                  }
                }
                if (sample(2, size = 1) == 1) 
                  Rpsp[, n:1, n:1] <- psp[, ncoord.x, ncoord.y]
                else Rpsp <- psp[, ncoord.x, ncoord.y]
            }
        }
        else {
            take <- sample(n^2)
            for (j in 1:S) {
                Rpsp[j, , ] <- matrix(psp[j, , ][take], ncol = n, 
                  nrow = n)
            }
        }
    }
    if (flag) 
        Rpsp <- drop(Rpsp)
    Rpsp
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
