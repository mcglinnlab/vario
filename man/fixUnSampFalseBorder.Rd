\name{fixUnSampFalseBorder}
\alias{fixUnSampFalseBorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fixUnSampFalseBorder(oarray, rarray)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oarray}{
%%     ~~Describe \code{oarray} here~~
}
  \item{rarray}{
%%     ~~Describe \code{rarray} here~~
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
function (oarray, rarray) 
{
    rarray.tmp = rarray
    if (length(dim(oarray)) == 3) {
        n2 = dim(oarray)[2]
        if (-999 \%in\% oarray[1, -c(1, n2), -c(1, n2)]) {
            S = dim(oarray)[1]
            o.na = oarray == -999
            r.na = rarray == -999
            end.tmp = which(o.na[1, -c(1, n2), -c(1, n2)])
            for (i in 1:S) {
                start.tmp = which(r.na[i, -c(1, n2), -c(1, n2)])
                if (sum(which(!is.na(match(end.tmp, start.tmp)))) > 
                  0) {
                  end = end.tmp[-which(!is.na(match(end.tmp, 
                    start.tmp)))]
                  start = start.tmp[-which(!is.na(match(start.tmp, 
                    end.tmp)))]
                }
                else {
                  end = end.tmp
                  start = start.tmp
                }
                rarray[i, -c(1, n2), -c(1, n2)][end] = rarray.tmp[i, 
                  -c(1, n2), -c(1, n2)][start]
                rarray[i, -c(1, n2), -c(1, n2)][start] = rarray.tmp[i, 
                  -c(1, n2), -c(1, n2)][end]
            }
        }
    }
    else {
        n2 = dim(oarray)[1]
        if (-999 \%in\% oarray[-c(1, n2), -c(1, n2)]) {
            o.na = oarray == -999
            r.na = rarray == -999
            end.tmp = which(o.na[-c(1, n2), -c(1, n2)])
            start.tmp = which(r.na[-c(1, n2), -c(1, n2)])
            if (sum(which(!is.na(match(end.tmp, start.tmp)))) > 
                0) {
                end = end.tmp[-which(!is.na(match(end.tmp, start.tmp)))]
                start = start.tmp[-which(!is.na(match(start.tmp, 
                  end.tmp)))]
            }
            else {
                end = end.tmp
                start = start.tmp
            }
            rarray[-c(1, n2), -c(1, n2)][end] = rarray.tmp[-c(1, 
                n2), -c(1, n2)][start]
            rarray[-c(1, n2), -c(1, n2)][start] = rarray.tmp[-c(1, 
                n2), -c(1, n2)][end]
        }
    }
    return(rarray)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
