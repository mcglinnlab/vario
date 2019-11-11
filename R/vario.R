library.dynam("vario", package="vario", lib.loc=.libPaths()[1])

## Spatial analysis functions---------------------------------------------------
#' Calculate uni- and multi-variate variograms
#'
#' This code is largely modified from the 'vegan' function 'mso' by Helene
#' Wagner, also parts of this code were derived from the 'geoR' function 'vairog'
#' by Paulo J. Ribeiro Jr. and Peter J. Diggle. 
#' Citations:
#' Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre,
#'   Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M.
#'   Henry H. Stevens and Helene Wagner (2011). vegan: Community Ecology
#'   Package. R package version 2.0-2.
#'   http://CRAN.R-project.org/package=vegan
#' Paulo J. Ribeiro Jr & Peter J. Diggle geoR: a package for
#'   geostatistical analysis R-NEWS, 1(2):15-18. June, 2001
#' 
#' Note: that if some areas are unsampled that obs != pos + neg + exp
#' @param x an object of class 'sim' or may be a sitexsp matrix, a vector of
#'   values, if x is the later then missing samples should be coded as -999
#' @param coord the spatial coordinates
#' @param grain interval size for distance classes, only used if 'breaks' not
#'   supplied
#' @param breaks the spatial breaks that define the spatial lags to be compared
#' @param hmax the maximum spatial lag
#' @param round.up if TRUE the spatial lags are rounded up to nearest integer
#' @param pos.neg if TRUE the positive and negative parts of the covariance
#'   matrix are output
#' @param binary if TRUE and x is class sim then abundances are converted to 
#'   binary 0 1 values
#' @param snap indicates which generation from an object class 'sim' to draw
#'   from
#' @param median indicates if in addition to the mean the medians of the
#'   distance matrices are calculated
#' @param direction a numerical value for the directional (azimuth) angle. This 
#'   used to specify directional variograms. Default defines the omnidirectional
#'   variogram. The value must be in the interval [0, pi] radians ([0, 180]
#'   degrees).
#' @param tolerance numerical value for the tolerance angle, when computing 
#'   directional variograms. The value must be in the interval [0, pi/2] radians
#'   ([0, 90] degrees).  Defaults to pi/8.
#' @param unit.angle defines the unit for the specification of angles in the two
#'   previous arguments. Options are 'radians' and 'degrees', with default to
#'   'radians'.
#' @param distance.metric can be one of the speices turnover metrics listed by
#'   the vegan function vegdist(). This is only appropriate if pos.neg = FALSE. 
#'   Common options include, 'jaccard' and 'bray'. If computed on pres/abse data
#'   then soreson index is computed by 'bray'.
#' @param min minimum value
#' @param max maximum value
#' @export
#' @examples
#' x = matrix(runif(100), ncol=10, nrow=10)
#' coord = cbind(1:10, 1)
#' vario(x, coord)
#' \dontrun{
#' library(vegan)
#' data(mite)
#' data(mite.xy)
#' vario(as.matrix(mite), mite.xy)
#' vario(as.matrix(mite), mite.xy, pos.neg=T)
#' }
vario = function(x, coord, grain=1, breaks=NA, log=FALSE, hmin=NA,
                 hmax=NA, round.int=FALSE, pos.neg=FALSE, binary=TRUE,
                 snap=NA, median=FALSE, quants=NA, direction = 'omnidirectional',
                 tolerance = pi/8, unit.angle = c('radians', 'degrees'),
                 distance.metric = 'euclidean', univariate=FALSE){
    if (class(x) == "sim"){
        coord = x$coords
        if (is.na(snap))
            snap = length(sim$snaps)
    }
    else if (is.vector(x))
        x = matrix(x)
    else if (!is.matrix(x))
        stop('x must be either of class sim, a vector, or a matrix')
    #x = ifelse(x == -999, NA, x)  ## best to fix these before entering into the function
    if (univariate) {
        vobject = vario_uni(x, bisect=FALSE, coord, grain, breaks, log, hmin, hmax,
                            round.int, pos.neg, binary, snap, median, quants, direction,
                            tolerance, unit.angle, distance.metric)
    }
    else {
        if (distance.metric != 'euclidean') {
            if (pos.neg)
                stop("cannot commpute pos-neg covariance using a turnover metric")
            else
                require(vegan)
        }
        unit.angle = match.arg(unit.angle)
        check_vario_direction_args(direction, tolerance, unit.angle)
        Dist = dist(coord)
        maxDist = max(Dist)
        if (is.na(breaks[1])) {
            if (is.na(hmin))
                hmin = grain
            if (is.na(hmax))
                hmax = round((maxDist / 2) / grain) * grain
            H = round(Dist / grain) * grain
        }
        else {
            if (is.na(hmin))
                hmin = min(Dist)  ## potentially this should be set when breaks is NA as well
            if (is.na(hmax))
                hmax = maxDist / 2
            H = Dist
            breaks = get_breaks(breaks, hmin, hmax, maxDist, log)
            if (round.int)
                breaks = round(breaks)
            for (i in 1:(length(breaks) - 1)) {
                H[H >= breaks[i] & H < breaks[i + 1]] = breaks[i]
            }  
        }
        H[H < hmin] = NA
        H[H > hmax] = NA
        H = as.vector(H)
        if (is.vector(x)) {
            S = 1
            N = length(x)
        }
        else {
            S = ncol(x)
            N = nrow(x)
        } 
        vobject = list()
        class(vobject) = 'vario'
        vobject$parms = data.frame(grain, hmin, hmax, S=S, N=N, pos.neg, median, direction,
                                   tolerance, unit.angle, distance.metric, 
                                   quants = ifelse(is.na(quants[1]), NA, 
                                                   paste(quants* 100, collapse=", ")))
        if(class(x) == "sim"){
            if(binary)
                x = apply(census(x, snap=snap), c(1, 2), as.logical) * 1
            else
                x = census(x, snap=snap)
            vobject$parms = cbind(vobject$parms, niche.wid.rel=x$p$s.rel,
                                  disp.wid.rel=x$p$u.rel) 
        }
        ## geoR code from function variog with slight modifications starts here'
        if (direction != "omnidirectional") {
            ## note that the changes: 'u' has been changed for as.vector(Dist) and
            ## coords changed to coord
            u.ang = .C("tgangle", as.double(as.vector(coord[ , 1])),
                       as.double(as.vector(coord[ , 2])), as.integer(dim(coord)[1]),
                       res = as.double(rep(0, length(as.vector(Dist)))), PACKAGE = 'vario')$res
            if (any(is.na(u.ang)))
                stop("NA returned in angle calculations maybe due to co-located data")
            u.ang = atan(u.ang)
            u.ang[u.ang < 0] = u.ang[u.ang < 0] + pi
            
            if (unit.angle == "degrees") {
                ang.deg = direction
                ang.rad = (ang.deg * pi) / 180
                tol.deg = tolerance
                tol.rad = (tol.deg * pi) / 180
            }
            else {
                ang.rad = direction
                ang.deg = (ang.rad * 180) / pi
                tol.rad = tolerance
                tol.deg = (tol.rad * 180) / pi
            }
            ang.lower = ang.rad - tol.rad
            ang.upper = ang.rad + tol.rad
            if (ang.lower >= 0 & ang.upper < pi)
                ang.ind = (!is.na(u.ang) & ((u.ang >= ang.lower) & (u.ang <= ang.upper)))
            if (ang.lower < 0)
                ang.ind = (!is.na(u.ang) & ((u.ang < ang.upper) | (u.ang > (pi + ang.lower))))
            if (ang.upper >= pi)
                ang.ind = (!is.na(u.ang) & ((u.ang > ang.lower) | (u.ang < (ang.upper - pi))))
            Dist[!ang.ind] = NA
            H[!ang.ind] = NA
        }
        ## geoR code from function variog ends here
        Dist = sapply(split(Dist, H), mean, na.rm=TRUE)
        vobject$vario = data.frame(H = as.numeric(names(table(H))), Dist = Dist,
                                   n = as.numeric(table(H)))
        ## below 'exp.gamma' is the expected variogram if 'x' is a sitexsp
        ## pres/abse matrix. The expectation is based upon the assumption of zero
        ## sp x sp covariances
        if(distance.metric == 'euclidean')
            exp.split = split(dist(x)^2 * .5, H)
        else
            exp.split = split(vegdist(x, method=distance.metric), H)
        exp.gamma = sapply(exp.split, mean, na.rm=TRUE)
        if (!is.na(quants[1])) {
            exp.qt = sapply(exp.split, function(x) quantile(x, quants, na.rm=TRUE))
            exp.qt = t(exp.qt)
            colnames(exp.qt) = paste(quants * 100)
        }  
        vobject$vario = cbind(vobject$vario, exp.var=exp.gamma)
        if (median)
            exp.med = sapply(exp.split, median, na.rm=TRUE)
        if (!is.vector(x)) { ## i.e. x is a site x sp matrix and not simply a vector
            ## if 'x' is a sitexsp pres/abse matrix the following computes site species richness
            rich = apply(x, 1, sum)
            ## see equation 7 in Wagner, H. 2003. Spatial covariance in plant
            ## communities... Ecology 84:1045-1057 to see that observed multivariate
            ## variogram can be computed from the species richness vector
            obs.gamma = sapply(split(dist(rich)^2 * .5, H), mean, na.rm=TRUE)
            vobject$vario = cbind(vobject$vario, obs.var = obs.gamma,
                                  ratio = obs.gamma/exp.gamma)
            if (pos.neg) {
                cov.mat = getCovFractions(x)
                pos.split = split(cov.mat$pos, H)
                neg.split = split(cov.mat$neg, H)
                pos = sapply(pos.split, mean)
                neg = sapply(neg.split, mean)
                vobject$vario = cbind(vobject$vario, pos = pos, neg = neg)
                if (median) {
                    pos.med = sapply(pos.split, median)
                    neg.med = sapply(neg.split, median)
                    vobject$vario = cbind(vobject$vario, exp.med = exp.med, pos.med = pos.med,
                                          neg.med = neg.med)
                }
                ##Note: obs.var = exp.var + pos.var + neg.var
            }
            else{
                if (median) {
                    obs.gamma.med = sapply(split(dist(rich)^2 * .5, H),median, na.rm=TRUE)
                    vobject$vario = cbind(vobject$vario, obs.med = obs.gamma.med,
                                          exp.med = exp.med)
                }  
            }
        }
        if (!is.na(quants[1]))
            vobject$vario = cbind(vobject$vario, exp.qt = exp.qt)
        if (is.vector(x))
            vobject$p = sum(x) / length(x)
        else
            vobject$p = apply(x, 2, sum, na.rm=TRUE) / nrow(x)
        vobject$perm = FALSE
        ## the above line indicates to the function 'vGraph' if the variogram is the
        ## result of a randomization or not
        row.names(vobject$vario) = NULL
    }
    return(vobject)
}

#' Computes univariate variograms for evey column of input matrix in addition to
#' the multivariate variogram
#' 
#' @param x site x sp matrix
#' @param bisect if the bisection style vairogram should be computed
#' @param ... arguments supplied to the function vario() Developer note speed 
#'   gains would be significant if partitioning of computation between species 
#'   was carried out within the vario function after computing the distance 
#'   matrix because that is a time intensive step
#' @export
#' @examples
#' x = matrix(runif(100), ncol=10, nrow=10)
#' coord = cbind(1:10, 1)
#' vario_uni(x, coord=coord) 
#' # alternatively call from vario() 
#' vario(x, coord, univariate=T)
vario_uni = function(x, bisect=FALSE, ...){
    require(snowfall)
    S = ncol(x)
    if (bisect)
        v = vario_bisect(x, ...)
    else
        v = vario(x, ...)
    n_cpus = length(suppressMessages(sfGetCluster()))
    if (n_cpus > 0) {
        sfLibrary(vario)
        sfLibrary(vegan)
        if (bisect) 
            exp_var = sfSapply(1:S, function(sp) vario_bisect(x[ , sp], ...)$vario$exp.var)
        else
            exp_var = sfSapply(1:S, function(sp) vario(x[ , sp], ...)$vario$exp.var)
    }
    else {
        if (bisect)
            exp_var = sapply(1:S, function(sp) vario_bisect(x[ , sp], ...)$vario$exp.var) 
        else
            exp_var = sapply(1:S, function(sp) vario(x[ , sp], ...)$vario$exp.var) 
    }
    colnames(exp_var) = paste('sp', 1:S, sep='')
    v$exp.var = exp_var  
    return(v)
}

#' Compute spatial breaks of distance bins used by the function vario()
#' 
#' @param breaks either a vector of breaks or an integer number of breaks to 
#'   compute. If a vector of breaks is supplied then the function returns that 
#'   exact same vector back.
#' @param hmin minimum distance of interest
#' @param hmax maximum distance of interest
#' @param maxDist maximum possible distance (greater or equal to hmax)
#' @param log boolean, if true the breaks are equidistance on a log scale
#' @export
#' @examples
#' get_breaks(4, 1, 16, 16)
#' get_breaks(4, 1, 16, 16, log=T)
get_breaks = function(breaks, hmin, hmax, maxDist, log=FALSE) {
    if (length(breaks) == 1) {
        if (log) {
            if (round(hmax, 2) == round(maxDist / 2, 2)) {
                incre = (hmax - hmin) / breaks
                hmax = hmax + incre
            }
            breaks = exp(seq(log(hmin), log(hmax), length.out=breaks))
        }
        else 
            breaks = seq(hmin, hmax, length.out=breaks)    
    }
    return(breaks)
}

#' Internal function that carries out checks on the directional arguments that
#' are supplied to the function vario(), if these checks are failed then vario()
#' will stop with an error message
#' Note: this code was copied from the package geoR in the function variog
check_vario_direction_args = function(direction = 'omnidirectional',
                                      tolerance = pi/8,
                                      unit.angle = c('radians', 'degrees')) {
    if (mode(direction) == "numeric") {
        if (length(direction) > 1)
            stop("only one direction is allowed")
        if (length(tolerance) > 1)
            stop("only one tolerance value is allowed")
        if (unit.angle == "degrees") {
            ang.deg = direction
            ang.rad = (ang.deg * pi) / 180
            tol.deg = tolerance
            tol.rad = (tol.deg * pi) / 180
        }
        else {
            ang.rad = direction
            ang.deg = (ang.rad * 180) / pi
            tol.rad = tolerance
            tol.deg = (tol.rad * 180) / pi
        }
        if (ang.rad > pi | ang.rad < 0)
            stop("direction must be an angle in the interval [0,pi[ radians")
        if (tol.rad > pi/2 | tol.rad < 0)
            stop("tolerance must be an angle in the interval [0,pi/2] radians")
        if (tol.deg >= 90) {
            direction = "omnidirectional"
        }
    }  
}

#' Calculate the lower diagonal of a sp covariance matrix
#' 
#' @return positive and negative fractions of covariance as two lower triangular
#'   matrices, each in vector format
#' @param x is a sitexsp matrix (sp as columns) of real numbers rows are the 
#'   sites, columns are the species
#' @export
#' @examples
#' x = matrix(runif(100), 10, 10)
#' getCovFractions(x)
getCovFractions = function(x) {
    N = as.integer(nrow(x))
    S = as.integer(ncol(x)) 
    x = as.double(ifelse(is.na(x) | x == -999,-99999,x))
    pos = as.double(rep(0,(N*(N-1))/2))
    neg = as.double(rep(0,(N*(N-1))/2))
    result = .C('loopcovreal',x,N,S,pos,neg, PACKAGE = 'vario')
    out = list()
    out$pos = result[[4]]
    out$neg = result[[5]]
    return(out)
} 

#' Estimates the fractal dimension D from a 2-dimensional grid
#'
#' @param x is a matrix of real numbers
#' @param coord the spatial coordinates
#' @export
#' @examples
#' x = matrix(runif(100), 10, 10)
#' estimateD(x)
estimateD = function(x, coord=NULL) {
    n = dim(x)[1]
    if(is.null(coord)) 
        coord= expand.grid(1:n, 1:n)
    v = vario(as.vector(x), coord)$vario
    mod = lm(log(v$exp) ~ log(v$Dist))
    m = coef(mod)[2]
    D = (6 - m) / 2
    return(D)
}

## Permutation functions -------------------------------------------------------

#' Internal function that maintains the spatial locations of the unsampled
#' pixels in a random realization of a two dimentioal spatial array. -999 is the
#' identifier for unsampled cells, in this case oarray and rarray DO have a
#' false border of -999
#' @param oarray observed starting array
#' @param rarray randomized array
fixUnSampFalseBorder = function(oarray,rarray){
    rarray.tmp = rarray
    if(length(dim(oarray)) == 3){ ## if multiple species then
        n2 = dim(oarray)[2]
        if(-999 %in% oarray[1,-c(1,n2),-c(1,n2)]){ 
            S = dim(oarray)[1]
            o.na = oarray == -999
            r.na = rarray == -999
            end.tmp = which(o.na[1,-c(1,n2),-c(1,n2)])  
            for(i in 1:S){
                start.tmp = which(r.na[i,-c(1,n2),-c(1,n2)]) 
                if(sum(which(!is.na(match(end.tmp,start.tmp)))) > 0){
                    ## drop ones in which end.tmp and start.tmp match
                    end = end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
                    start = start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
                }
                else{
                    end = end.tmp
                    start = start.tmp
                }
                rarray[i,-c(1,n2),-c(1,n2)][end] = rarray.tmp[i,-c(1,n2),-c(1,n2)][start]
                rarray[i,-c(1,n2),-c(1,n2)][start] = rarray.tmp[i,-c(1,n2),-c(1,n2)][end]
            } } }
    else{  ## only a single species
        n2 = dim(oarray)[1]
        if(-999 %in% oarray[-c(1,n2),-c(1,n2)]){ 
            o.na = oarray == -999
            r.na = rarray == -999
            end.tmp = which(o.na[-c(1,n2),-c(1,n2)])  
            start.tmp = which(r.na[-c(1,n2),-c(1,n2)])  
            if(sum(which(!is.na(match(end.tmp,start.tmp)))) > 0){
                ## drop ones in which end.tmp and start.tmp match
                end = end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
                start = start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
            }
            else{
                end = end.tmp
                start = start.tmp
            }
            rarray[-c(1,n2),-c(1,n2)][end] = rarray.tmp[-c(1,n2),-c(1,n2)][start]
            rarray[-c(1,n2),-c(1,n2)][start] = rarray.tmp[-c(1,n2),-c(1,n2)][end]
        } }
    return(rarray)
}

#' Internal function that maintains the spatial locations of the unsampled 
#' pixels in a random realization of a two dimentioal spatial array. -999 is the
#' identifier for unsampled cells, in this case oarray and rarray DO NOT have a 
#' false border of -999
#' @param oarray observed starting array
#' @param rarray randomized array
fixUnSampTrueBorder = function(oarray,rarray){
    if(-999 %in% oarray){  
        rarray.tmp = rarray
        if(length(dim(oarray)) == 3){  ## if multiple species then
            S = dim(oarray)[1]
            n = dim(oarray)[2]
            o.na = oarray == -999
            r.na = rarray == -999
            end.tmp = which(o.na[1,,]) 
            for(i in 1:S){
                start.tmp = which(r.na[i,,]) 
                if(sum(which(!is.na(match(end.tmp,start.tmp)))) > 0){
                    ## drop ones in which end.tmp and start.tmp match
                    end = end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
                    start = start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
                }
                else{
                    end = end.tmp
                    start = start.tmp
                }
                rarray[i,,][end] = rarray.tmp[i,,][start]
                rarray[i,,][start] = rarray.tmp[i,,][end]
            } }
        else{  ## only a single species
            n2 = dim(oarray)[1]
            o.na = oarray == -999
            r.na = rarray == -999
            end.tmp = which(o.na)  
            start.tmp = which(r.na)  
            if(sum(which(!is.na(match(end.tmp,start.tmp)))) > 0){
                ## drop ones in which end.tmp and start.tmp match
                end = end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
                start = start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
            }
            else{
                end = end.tmp
                start = start.tmp
            }
            rarray[end] = rarray.tmp[start]
            rarray[start] = rarray.tmp[end]
        } }
    return(rarray)
}

#' Mediate the generation of statistical null values for the variograms
#' 
#' This function is to be used in a parrallel processing loop which will 
#' generate a population of null values
#' 
#' @param pop the species x row x col array where row and column refer to 
#'   spatial location
#' @param vobject the output of the vario function with serves as the basis for 
#'   empirical comparisons
#' @param coords the spatial coordinates of a M^2 x S matrix
#' @param meth the type of permutation to use, options include: reflect: random 
#'   reflection/rotations of species (only makes sence when sp are not fixed 
#'   shift: random torodial shifting with or with sp fixed both: both reflection
#'   and shifting random: random shuffle randpat: random patterns algo of 
#'   Roxburgh and Chesson 1998, must parameterize RPargs (See below)
#' @param sp boolean, if FALSE then obs composition of quadrats is fixed to the 
#'   observed pattern, if TRUE then species are each shuffled independently
#' @param all boolean, if TRUE then all relevant nulls are calculated
#' @param RPargs a list of arguments that must be supplied if the random 
#'   patterns algo is desired, the arguments of RPargs are input into the 
#'   function 'randPatPar', they include: 1) allRP if TRUE & all = TRUE, then 
#'   Random Patterns algo used as the spatial null, 2)nstrata, 3)mtrials1, 
#'   4)mtrials2, 5)alpha, 6)npar
#' @param median: if TRUE then means and medians are calculated Note: 'meth' and
#'   'sp' are arguments to randomization function 'spatPerm2D'
#' @export
nullGen = function(pop,vobject,coords,meth,sp,all=FALSE,RPargs=FALSE,median=FALSE) {
    S = dim(pop)[1]
    n = dim(pop)[2]
    n2 = n+2
    grain = vobject$parms$grain
    hmax = vobject$parms$hmax
    pos.neg = vobject$parms$pos.neg
    median = vobject$parms$median
    if(class(vobject$parms$direction) == "factor") 
        direction = as.character(vobject$parms$direction)
    else
        direction = as.numeric(vobject$parms$direction)
    tolerance = vobject$parms$tolerance
    unit.angle = as.character(vobject$parms$unit.angle)
    if(all){
        rmat = apply(pop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
        rv = vobject$vario
        if(pos.neg){
            if(median) 
                r.vals = array(NA,dim=c(length(rv$Dist),6,3)) ##dists,6 results, 3 methods
            else
                r.vals = array(NA,dim=c(length(rv$Dist),3,3)) ##dists,3 results, 3 methods
        }
        else
            r.vals = array(NA,dim=c(length(rv$Dist),2,3)) ##dists, 2 results, 3 methods
        ###need to add in logicals for handling median values here
        for(j in 1:3){ #cycle through permutation methods
            if(j == 1){
                rpop = spatPerm2D(pop,meth='random',sp=TRUE)##random model
                rpop = fixUnSampTrueBorder(pop,rpop)
            }
            if(j == 2){
                if(RPargs[[1]]){ ##use the random pattern algo for the spatial null
                    out = randPatPar(psp=pop,nstrata=RPargs[[2]],mtrials1=RPargs[[3]],mtrials2=RPargs[[4]],alpha=RPargs[[5]],npar=RPargs[[6]])
                    rpop = array(0,dim=c(S,n2,n2))
                    for(i in 1:S){
                        rpop[i,,] = array(out[-(1:5),i],dim=c(n2,n2))
                    }
                    rpop = rpop[,-c(1,n2),-c(1,n2)]
                }
                else
                    rpop = spatPerm2D(pop,meth='both',sp=TRUE)##Random Shifts & Reflections spatial null
            }
            if(j == 3){
                rpop = spatPerm2D(pop,meth='random',sp=FALSE) ##species model
                rpop = fixUnSampTrueBorder(pop,rpop)
            }
            rmat = apply(rpop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
            rv<-vario(x=rmat,coord=coords,grain=grain,breaks=breaks,hmin=hmin,hmax=hmax,
                      pos.neg=pos.neg,median=median,direction=direction,tolerance=tolerance,
                      unit.angle=unit.angle,distance.metric=distance.metric)$vario
            if(pos.neg){
                if(median)
                    r.vals[,,j] = as.matrix(rv[,c(5,7:11)])
                else
                    r.vals[,,j] = as.matrix(rv[,c(5,7:8)])
            }
            else
                r.vals[,,j] = as.matrix(rv[,4:5])
        }}
    else{ ##not using every null model, only one null
        if(RPargs[[1]]){ ##use the random pattern algo for the spatial null
            out = randPatPar(psp=pop,nstrata=RPargs[[2]],mtrials1=RPargs[[3]],
                             mtrials2=RPargs[[4]],alpha=RPargs[[5]],npar=RPargs[[6]])
            rpop = array(0,dim=c(S,n2,n2))
            for(i in 1:S){
                rpop[i,,] = array(out[-(1:5),i],dim=c(n+2,n+2))
            }
            rpop = rpop[,-c(1,n+2),-c(1,n+2)]
        } 
        else{
            rpop = spatPerm2D(pop,meth=meth,sp=sp)
            rpop = fixUnSampTrueBorder(pop,rpop)
        }
        rmat = apply(rpop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
        rv<-vario(x=rmat,coord=coords,grain=grain,breaks=breaks,hmin=hmin,hmax=hmax,
                  pos.neg=pos.neg,median=median,direction=direction,tolerance=tolerance,
                  unit.angle=unit.angle,distance.metric=distance.metric)$vario
        if(pos.neg){
            if(median)
                r.vals = as.matrix(rv[,c(5,7:11)])
            else
                r.vals = as.matrix(rv[,c(5,7:8)])
        }
        else
            r.vals = as.matrix(rv[,4:5])
    }
    return(r.vals)
}

#' Generate statistical null expectations for the variograms
#' 
#' @param x is either an output of class 'sim' that is the output of 
#'   'sim.neut.uni' OR an site x species matrix
#' @param vobject is the output of 'vario', this informs the function of what 
#'   parameterization of vario to use specifically it indiates if the pos.neg 
#'   components and median should be calculated
#' @param nperm is the number of permutations
#' @param coords the spatial coordinates of the sites, not needed if x is of 
#'   class 'sim'
#' @param meth the type of permutation to use, options include: reflect: random 
#'   reflection/rotations of species (only makes sence when sp are not fixed 
#'   shift: random torodial shifting with or with sp fixed both: both reflection
#'   and shifting random: random shuffle randpat: random patterns algo of 
#'   Roxburgh and Chesson 1998, must parameterize RPargs (See below)
#' @param sp boolean, if FALSE then obs composition of quadrats is fixed to the 
#'   observed pattern, if TRUE then species are each shuffled independently
#' @param all boolean, if TRUE then all relevant nulls calculated
#' @param snap is the time period of the simulation to analyze,defaults to NULL 
#'   val, if not set gets internally set to last time period
#' @param npar number of processors to run the function on
#' @param median if TRUE then median is also calculated in addition to mean
#' @param linux if TRUE then function assumes you are on a linux cluster and 
#'   therefore exports a different compiled code
#' @param RPargs is a vector of arguments that are needed to perform the Random 
#'   Patterns spatial null model
#' @param breaks gives either the number or the position of the breaks for the 
#'   function vario Notes: 1) see the notes associated with the function 
#'   'nullGen' that indicate how 'RPargs" should be parameterized 2) 'meth' and 
#'   'sp' are arguments to randomization function 'spatPerm2D'
#' @export
nullPerm = function(x, vobject, nperm, coords=NULL, meth='both',
                    sp=TRUE, all=FALSE, snap=NULL, npar=1, linux=FALSE,
                    RPargs=FALSE, breaks=NA) {
    dists = vobject$vario$Dist
    grain = vobject$parms$grain
    hmin = vobject$parms$hmin
    hmax = vobject$parms$hmax
    pos.neg = vobject$parms$pos.neg
    median = vobject$parms$median
    if (class(vobject$parms$direction) == "factor") 
        direction = as.character(vobject$parms$direction)
    else
        direction = as.numeric(vobject$parms$direction)
    tolerance = vobject$parms$tolerance
    unit.angle = as.character(vobject$parms$unit.angle)
    if (is.na(unit.angle))
        unit.angle = 'degrees'
    distance.metric = as.character(vobject$parms$distance.metric)
    if (class(x) == 'sim') {
        coords = x$coords
        if (is.null(snap))
            snap = length(x$snaps)
    }
    else if (is.null(coords))
        stop('need to supply spatial coordinates if not a simulation product')
    r.vals = list()
    r.vals$parms = vobject$parms
    r.vals$p = vobject$p
    if (median & !pos.neg)
        stop("if computing medians must also compute pos.neg fractions, set pos.neg=TRUE")
    if (pos.neg) { ##pos and neg fractions
        if (all) { ##all relevant null models
            if (median) { ##will compute mean and median
                r.vals$vario = array(0, dim=c(length(dists), 6, 3, nperm + 1))##dists, results, methods, perms
                r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos',
                                                                      'neg', 'exp.med',
                                                                      'pos.med', 'neg.med')]) 
            }
            else { ##only compute means
                r.vals$vario = array(0, dim=c(length(dists), 3, 3, nperm+1))##dists, results, methods, perms
                r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos', 'neg')])
            }
        }  
        else { ##only a single null used
            if (median) {
                r.vals$vario = array(0, dim=c(length(dists), 6, nperm + 1))##dists, results, perms
                r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos',
                                                                    'neg', 'exp.med',
                                                                    'pos.med', 'neg.med')]) 
            }
            else {
                r.vals$vario = array(0, dim=c(length(dists), 3, nperm + 1))##dists, results, perms
                r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos', 'neg')])
            }
        }
    }
    else {##only exp and obs fractions
        if (all) {
            if (median) { ##will compute mean and median
                r.vals$vario = array(0, dim=c(length(dists), 4, 3, nperm + 1))
                r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var',
                                                                      'obs.med', 'exp.med')])
            }
            else {
                r.vals$vario = array(0, dim=c(length(dists), 2, 3, nperm + 1))
                r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var')])
            }
        }
        else { ##only a single null used
            if (median) {
                r.vals$vario = array(0, dim=c(length(dists), 4, nperm + 1))
                r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var',
                                                                    'obs.med', 'exp.med')])
            }
            else { 
                r.vals$vario = array(0, dim=c(length(dists), 2, nperm + 1))
                r.vals$vario[ , , 1] = as.matrix(vobject$vario[, c('exp.var', 'obs.var')])
            }
        }
    }
    if (class(x) == 'sim') {
        pop = as.logical(x$snaps[[snap]]) ##converts it to a pres/absence vector
        dim(pop) = c(x$p$S, x$p$M, x$p$M)
    }
    else {
        pop = array(x, dim=c(sqrt(nrow(x)), sqrt(nrow(x)), ncol(x)))
        pop = aperm(pop, c(3, 1, 2))
    }
    if (RPargs[[1]] & npar == 1) {
        r.vals$p.conv1 = 0 ##average proportion of species that converged with strata swaps
        r.vals$p.conv2 = 0 ##average proportion of species that converged with pixel swaps
    }
    if (npar == 1) { ##all permutations option not yet implemented for 1 processor
        pb = txtProgressBar(min = 0, max = nperm, style = 3)
        for (i in 1:nperm) {
            if (RPargs[[1]]) { ##use the random pattern algo for the spatial null
                out = RandPatPar(psp=pop, nstrata=RPargs[[2]], mtrials1=RPargs[[3]],
                                 mtrials2=RPargs[[4]], alpha=RPargs[[5]], npar=RPargs[[6]])
                S = dim(pop)[1]
                n2 = dim(pop)[2] + 2
                r.vals$p.conv1 = r.vals$p.conv1 + sum(out[2, ] <= RPargs[[5]]) / S / nperm
                r.vals$p.conv2 = r.vals$p.conv2 + sum(out[4, ] <= RPargs[[5]], na.rm=TRUE) / S / nperm
                rpop = array(0, dim=c(S, n2, n2))
                for (k in 1:S) {
                    rpop[k, , ] = array(out[-(1:5), k], dim=c(n2, n2))
                }
                rpop = rpop[ , -c(1, n2), -c(1, n2)]
            }
            else {
                rpop = SpatPerm2D(pop, meth=meth, sp=sp)
                rpop = FixUnSamp2(pop, rpop)
            }
            rmat = apply(rpop, 1, as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
            rv = vario(x=rmat, coord=coords, grain=grain, breaks=breaks, hmin=hmin,
                       hmax=hmax, pos.neg=pos.neg, median=median, direction=direction,
                       tolerance=tolerance, unit.angle=unit.angle,
                       distance.metric=distance.metric)$vario
            if (pos.neg) {
                if (all) {
                    if (median) {
                        r.vals$vario[ , , , i + 1] = as.matrix(rv[ , ,
                                                                  c('obs.var','pos','neg',
                                                                    'exp.med','pos.med',
                                                                    'neg.med')])
                    }
                    else {
                        r.vals$vario[ , , , i + 1] = as.matrix(rv[ , , c('obs.var', 'pos', 'neg')])
                    }
                } 
                else {
                    if (median) {
                        r.vals$vario[ , , i + 1] = as.matrix(rv[ , ,
                                                                c('obs.var','pos','neg',
                                                                  'exp.med','pos.med',
                                                                  'neg.med')])
                    }  
                    else {
                        r.vals$vario[ , , i + 1] = as.matrix(rv[ , c('obs.var', 'pos', 'neg')])
                    }  
                }
            }
            else {
                if (all) {
                    if (median) {
                        r.vals$vario[ , , , i + 1] = as.matrix(rv[ , , c('exp.var', 'obs.var',
                                                                         'obs.med', 'exp.med')])
                    }
                    else {
                        r.vals$vario[ , , , i + 1] = as.matrix(rv[ , , c('exp.var', 'obs.var')])
                    }
                }
                else {
                    if (median) {
                        r.vals$vario[ , , i + 1] = as.matrix(rv[ , c('exp.var', 'obs.var',
                                                                     'obs.med', 'exp.med')])
                    }
                    else {
                        r.vals$vario[ , , i + 1] = as.matrix(rv[ , c('exp.var', 'obs.var')])
                    }
                }
            }
            #print(i)
            Sys.sleep(0.1)
            # update progress bar
            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    else { ##computing in parallel
        require(snowfall)
        sfInit(parallel=TRUE, cpus=npar, type="SOCK")
        sfClusterSetupRNG()
        sfExport("pop", "vobject", "coords", "meth", "all", "sp", "RPargs",
                 "breaks", "RandPatPar", "RandPat", "FixUnSamp", "FixUnSamp2",
                 "SpatPerm2D", "SpatPerm2D.str", "vario", "getCovFractions",
                 "null.gen", "check_vario_direction_args")
        if (linux)
            sfClusterEval(dyn.load("danspkg.so"))
        else {
            sfLibrary(danspkg)
        }
        out = unlist(sfLapply(1:nperm, function(...) null.gen(pop=pop, vobject=vobject,
                                                              coords=coords, meth=meth,
                                                              sp=sp, all=all, RPargs=RPargs,
                                                              breaks=breaks)))
        sfStop()
        if (pos.neg) {
            if (all) {
                if (median) {
                    dim(out) = c(length(dists), 6, 3, nperm)
                    r.vals$vario[ , , , -1] = out 
                }
                else {
                    dim(out) = c(length(dists), 3, 3, nperm)
                    r.vals$vario[ , , , -1] = out 
                }
            }
            else {
                if (median) {
                    dim(out) = c(length(dists), 6, nperm)  
                    r.vals$vario[ , , -1] = out     
                }
                else {
                    dim(out)=c(length(dists), 3, nperm)  
                    r.vals$vario[ , , -1] = out     
                }
            }
        }
        else {
            if (all) {
                if (median) {
                    dim(out) = c(length(dists), 4, 3, nperm)  
                    r.vals$vario[ , , , -1] = out 
                }
                else {
                    dim(out) = c(length(dists), 2, 3, nperm)  
                    r.vals$vario[ , , , -1] = out 
                }
            }
            else {
                if (median) {
                    dim(out) = c(length(dists), 4, nperm)  
                    r.vals$vario[ , , -1] = out 
                }
                else {
                    dim(out) = c(length(dists), 2, nperm)  
                    r.vals$vario[ , , -1] = out 
                }
            }
        }
    }
    r.vals$perm = TRUE
    r.vals$vdists = vobject$vario$Dist
    return(r.vals)
}

#' Random Patterns spatial permutation algorithm
#' 
#' This function evaulates the .C function 'randpatpar' which is the random
#' patterns algo of Roxburgh and Chesson 1998. Returns species index, phi stat,
#' number of actual swaps, and the randomized presences as a single vector of
#' numbers. To be called in serial or parallel by function "randPatPar'
#' @param i the ith species index
#' @param psp multidimenstional S x (n+2) x (n+2) array
#' @param rpsp a randomized version of psp
#' @param n the size of the orginal 2-D array along one spatial axis (i.e.,
#'   without extra rows and columns)
#' @param pl the places in rpsp that can be swaped
#' @param mtrials1 the number of times to attempt a swap at the strata level
#' @param mtrials2 the number of times to attempt a swap at the pixel level
#' @param alpha the cutoff value for the phi statistic of Roxburgh and Chesson
#'   1998
#' @export
randPat = function(i,psp,rpsp,n,nstrata,pl,mtrials1=1e3,mtrials2=1e6,
                   alpha=0.01){
    psp = psp[i,,]
    rpsp = rpsp[i,,]
    n2 = n+2
    ## PART I
    ## begin permuting the blocks defined by nstrata
    rpsp.tmp = rpsp[-c(1,n2),-c(1,n2)]
    coords = cbind(rep(1:nstrata,nstrata),rep(1:nstrata,each=nstrata))
    rcoords = coords[sample(nstrata^2),]
    for(j in 1:nstrata^2){
        rows = ((coords[j,1]-1)*n/nstrata+1) : ((coords[j,1]-1)*n/nstrata+n/nstrata)
        cols = ((coords[j,2]-1)*n/nstrata+1) : ((coords[j,2]-1)*n/nstrata+n/nstrata)
        rrows = ((rcoords[j,1]-1)*n/nstrata+1) : ((rcoords[j,1]-1)*n/nstrata+n/nstrata)
        rcols = ((rcoords[j,2]-1)*n/nstrata+1) : ((rcoords[j,2]-1)*n/nstrata+n/nstrata)
        rpsp[-c(1,n2),-c(1,n2)][rows,cols] = rpsp.tmp[rrows,rcols]
    }
    rpsp = fixUnSampFalseBorder(psp,rpsp)
    ostat = .C('spatstat',as.double(as.vector(psp)),as.integer(n),
               as.double(rep(0,4)),PACKAGE = "vario")[[3]]
    nstat = .C('spatstat',as.double(as.vector(rpsp)),as.integer(n),
               as.double(rep(0,4)),PACKAGE = "vario")[[3]]
    phi = .C('calcphi',as.double(nstat),as.double(ostat),as.double(0),
             PACKAGE = "vario")[[3]]
    ## now begin random swapping of blocks defined by strata 
    ntrials = 0 ; gtrials = 0
    rpsp.tmp1 = rpsp
    while(phi > alpha & ntrials < mtrials1){
        rpsp.tmp2 = rpsp.tmp1[-c(1,n2),-c(1,n2)]
        rcoords = coords[sample(nstrata^2,2),]
        startrows = ((rcoords[1,1]-1)*n/nstrata+1) : 
            ((rcoords[1,1]-1)*n/nstrata+n/nstrata)
        startcols = ((rcoords[1,2]-1)*n/nstrata+1) : 
            ((rcoords[1,2]-1)*n/nstrata+n/nstrata)
        endrows = ((rcoords[2,1]-1)*n/nstrata+1) : 
            ((rcoords[2,1]-1)*n/nstrata+n/nstrata)
        endcols = ((rcoords[2,2]-1)*n/nstrata+1) : 
            ((rcoords[2,2]-1)*n/nstrata+n/nstrata)
        rpsp.tmp1[-c(1,n2),-c(1,n2)][startrows,startcols] = rpsp.tmp2[endrows,endcols]
        rpsp.tmp1[-c(1,n2),-c(1,n2)][endrows,endcols] = rpsp.tmp2[startrows,startcols]
        rpsp.tmp1 = fixUnSampFalseBorder(psp,rpsp.tmp1)
        nstat = .C('spatstat',as.double(as.vector(rpsp.tmp1)),as.integer(n),
                   as.double(rep(0,4)),PACKAGE = "vario")[[3]]
        phiTemp = .C('calcphi',as.double(nstat),as.double(ostat),as.double(0),
                     PACKAGE = "vario")[[3]]
        if(phiTemp < phi){
            phi = phiTemp
            gtrials = gtrials +1
            ## make change permanent
            rpsp = rpsp.tmp1
        }
        else{
            ## start back with orginal random mat
            rpsp.tmp1 = rpsp
        }
        ntrials = ntrials+1
    } 
    if(phi > alpha & mtrials2>0){
        ## Part II
        ## carry out individual pixel swapping
        psp = as.vector(psp)
        rpsp = as.vector(rpsp)
        tmp = .C('randpatpar',psp = as.double(psp),
                 rpsp = as.double(rpsp), n = as.integer(n),
                 ostat = as.double(rep(0,4)), nstat = as.double(rep(0,4)), 
                 phi = as.double(0), phiTemp = as.double(0),
                 alpha = as.double(alpha), pl = as.integer(pl-1),
                 nplaces = as.integer(length(pl)-1),ntrials = as.double(0),
                 gtrials = as.double(0), mtrials = as.double(mtrials2),
                 PACKAGE = "vario")
        out = c(i,phi,gtrials,tmp$phi,tmp$gtrials,tmp$rpsp)
    }
    else
        out = c(i,phi,gtrials,NA,NA,as.vector(rpsp))
    return(out)
}

#' Random Patterns spatial permutation algorithm
#' 
#' A convience function for working with randPat which calls the .C function
#' 'randpatpar'. This function allows you to specify the number of processors to
#' run on adding processsors only helps if working with many species as each
#' species is evaulated on a different processor. Returns a (5+(n+2)^2) x S
#' matrix, the first five rows of which are species index, phi strata stat,
#' number of strata swaps, phi pixel swap, and number of pixel swaps, and then
#' the remaining rows are the presences/abundances in the randomized occurances
#' @param psp multidimenstional S x (n+2) x (n+2) array
#' @param n the size of the orginal 2-D array along one spatial axis (i.e.,
#'   without extra rows and columns)
#' @param pl the places in rpsp that can be swaped
#' @param mtrials the numbef of times to attempt a swap
#' @param alpha the cutoff value for the phi statistic of Roxburgh and Chesson
#'   1998
#' @param npar the number of processors to run the code on
#' @export
randPatPar = function(psp,nstrata,mtrials1=1e3,mtrials2=1e6,alpha=0.01,npar=1){
    n = dim(psp)[2]
    if(length(dim(psp)) == 3)
        S = dim(psp)[1]
    else{
        S = 1
        psp = array(psp,dim=c(S,n,n))
    }
    ## first prepare psp for the randomization process
    ## add border of -999
    n2 = n+2
    psp.temp = array(0,dim=c(S,n2,n2))  ## create an array with boundary cells
    psp.temp[,-c(1,n2),-c(1,n2)] = psp  ## populate the array with the input data
    psp.temp[,1,-c(1,n2)] = -999  ## x of 0
    psp.temp[,n2,-c(1,n2)] = -999  ## x of n+1
    psp.temp[,-c(1,n2),1] = -999  ## y of 0
    psp.temp[,-c(1,n2),n2] = -999  ## y of n+1
    psp.temp[,1,1] = -999
    psp.temp[,1,n2] = -999
    psp.temp[,n2,1] = -999
    psp.temp[,n2,n2] = -999
    psp = psp.temp
    pl = 1:n2^2
    c1 = NA
    c2 = NA
    skip = pl[-999 == as.vector(psp[1,,])]
    pl = pl[is.na(match(pl,skip))]
    ##now read to begin intital randomization
    rpsp = psp
    ##inital reflection/rotation for each species independently
    rpsp[,-c(1,n2),-c(1,n2)] = spatPermStrata(psp[,-c(1,n2),-c(1,n2)],meth='reflect',
                                              sp=TRUE,nstrata=nstrata)
    nplaces = length(pl)
    nswaps = (nplaces*(nplaces-1))/2
    if(npar>1){
        require(snowfall)
        sfInit(parallel=TRUE, cpus=npar, type='SOCK')
        sfClusterSetupRNG()
        sfExport('randPat','fixUnSampFalseBorder','psp','rpsp','n','nstrata','pl',
                 'mtrials1','mtrials2','alpha')
        sfClusterEval(dyn.load('randompatternspar.dll'))
        out = unlist(sfSapply(1:S,function(i){
            randPat(i=i,psp=psp,rpsp=rpsp,n=n,nstrata=nstrata,pl=pl,mtrials1=mtrials1,
                    mtrials2=mtrials2,alpha=alpha)}))
        sfStop()
    }
    else{
        out = NULL
        for(i in 1:S)
            out = cbind(out,randPat(i=i,psp=psp,rpsp=rpsp,n=n,nstrata=nstrata,pl=pl,
                                    mtrials1=mtrials1,mtrials2=mtrials2,alpha=alpha))
    }
    return(out)
}

#' Permute an array of occurances under a given set of constraints in 2-D space
#' 
#' @param psp the sp x row x col array, where rows and columns specify where on 
#'   the spatial grid the sample was located
#' @param shiftpos two numbers that are the x and y places to shift the grid, 
#'   this is generated randomly if needed
#' @param rotate a single number 1-4 that indicates how many counterclockwise 
#'   rotations to perform, generated randomly
#' @param meth the type of permutation to use, options include: 1) reflect:
#'   random, 2) reflection/rotations of species (only makes sence when sp are
#'   not fixed 3) shift: random torodial shifting with or with sp fixed 4) both:
#'   both reflection and shifting 5) random: random shuffle
#' @param sp if FALSE then obs composition of quadrats is fixed to the observed 
#'   pattern if TRUE then species are each shuffled independently
#' @export
spatPerm2D = function(psp,shiftpos=NULL,rotate=NULL,meth='shift',sp=FALSE) {
    n = dim(psp)[2]
    if(length(dim(psp)) == 3){
        S = dim(psp)[1]
        flag = FALSE
    }
    else{
        S = 1
        psp = array(psp,dim=c(S,n,n))
        flag = TRUE
    }
    Rpsp = psp
    if(sp){  ## then between sp associations nullified
        if(meth != 'reflect'){
            if(meth != 'random'){
                ## generate vectors of random shifts, one for the x- and one for y-coord
                for(j in 1:S){
                    if(is.null(shiftpos)){
                        shift.x = sample(n,size=1) ; shift.y = sample(n,size=1);
                    }
                    else{
                        shift.x = shiftpos[1] ; shift.y = shiftpos[2]
                    }
                    #gen new coords
                    if(shift.x == 1)
                        ncoord.x = 1:n
                    else
                        ncoord.x = c(shift.x:n,1:(shift.x-1))
                    if(shift.y == 1)
                        ncoord.y = 1:n
                    else
                        ncoord.y = c(shift.y:n,1:(shift.y-1))
                    if(meth == 'shift'){
                        ## begin rearranging the rows of matrix for jth sp
                        Rpsp[j,,] = psp[j,ncoord.x,ncoord.y]
                    }
                    if(meth == 'both'){  ## reflecting/rotating and shifting
                        if(is.null(rotate))
                            rotate = sample(4,size=1)  ## how many counterclockwise rotatations to make
                        if(rotate == 2){
                            for(x in 1:n){
                                for(y in 1:n){
                                    Rpsp[j,(n-y)+1,x] = psp[j,x,y]
                                } } }
                        if(rotate == 3){
                            for(x in 1:n){
                                for(y in 1:n){
                                    Rpsp[j,(n-x)+1,(n-y)+1] = psp[j,x,y]
                                } } }
                        if(rotate == 4){
                            for(x in 1:n){
                                for(y in 1:n){
                                    Rpsp[j,y,(n-x)+1] = psp[j,x,y]
                                } } }
                        flips = sample(2,replace=TRUE)  ## generates two coin flips
                        if(flips[1] == 1){  ## reflect along x-axis
                            if(flips[2] == 1)  ## reflect along y-axis
                                Rpsp[j,n:1,n:1] = Rpsp[j,ncoord.x,ncoord.y]
                            else  ## not reflected along y-axis
                                Rpsp[j,n:1,] = Rpsp[j,ncoord.x,ncoord.y]
                        }
                        else{  ## not reflected along x-axis
                            if(flips[2] == 1)  ## reflected along y-axis
                                Rpsp[j,,n:1] = Rpsp[j,ncoord.x,ncoord.y]
                            else  ## not reflected along either axis
                                Rpsp[j,,] = Rpsp[j,ncoord.x,ncoord.y]
                        } } } } } 
        if(meth == 'reflect'){  ## if only want reflecting/rotating
            for(j in 1:S){
                if(is.null(rotate))
                    rotate = sample(4,size=1)  ## how many counterclockwise rotatations to make
                if(rotate == 2){
                    for(x in 1:n){
                        for(y in 1:n){
                            Rpsp[j,(n-y)+1,x] = psp[j,x,y]
                        } } }
                if(rotate == 3){
                    for(x in 1:n){
                        for(y in 1:n){
                            Rpsp[j,(n-x)+1,(n-y)+1] = psp[j,x,y]
                        } } }
                if(rotate == 4){
                    for(x in 1:n){
                        for(y in 1:n){
                            Rpsp[j,y,(n-x)+1] = psp[j,x,y]
                        } } }
                flips = sample(2,replace=TRUE) ##generates two coin flips
                if(flips[1] == 1){ #reflect along x-axis
                    if(flips[2] == 1) #reflect along y-axis
                        Rpsp[j,n:1,n:1] = Rpsp[j,,]
                    else #not reflected along y-axis
                        Rpsp[j,n:1,] = Rpsp[j,,]
                }
                else{ #not reflected along x-axis
                    if(flips[2] == 1) #reflected along y-axis
                        Rpsp[j,,n:1] = Rpsp[j,,]
                } } }
        if(meth == 'random'){
            for(j in 1:S){
                take = sample(n^2) #sample w/o replacement
                Rpsp[j,,] = matrix(psp[j,,][take],ncol=n,nrow=n)
            } } }
    else{  ## species co-occurances are fixed
        ## this only makes sense for "shift" or "both" meth
        if(meth == 'reflect'){
            stop('Reflecting fixed species co-occurances w/o shifting is not meaningful')
        }
        ## generate vector of random shifts
        if(meth != 'random'){
            if(is.null(shiftpos)){
                shift.x = sample(n,size=1) ; shift.y = sample(n,size=1);
            }
            else{
                shift.x = shiftpos[1] ; shift.y = shiftpos[2]
            }
            ## gen new coords
            if(shift.x == 1) 
                ncoord.x = 1:n
            else 
                ncoord.x = c(shift.x:n,1:(shift.x-1))
            if(shift.y == 1)
                ncoord.y = 1:n
            else
                ncoord.y = c(shift.y:n,1:(shift.y-1))
            if(meth == 'shift'){
                ## begin rearranging the rows of matrix for jth sp
                Rpsp = psp[,ncoord.x,ncoord.y]
            }
            if(meth == 'both'){  ## reflecting/rotating and shifting
                if(is.null(rotate))
                    rotate = sample(4,size=1)  ## how many counterclockwise rotatations to make
                if(rotate == 2){
                    for(x in 1:n){
                        for(y in 1:n){
                            Rpsp[,(n-y)+1,x] = psp[,x,y]
                        } } }
                if(rotate == 3){
                    for(x in 1:n){
                        for(y in 1:n){
                            Rpsp[,(n-x)+1,(n-y)+1] = psp[,x,y]
                        } } }
                if(rotate == 4){
                    for(x in 1:n){
                        for(y in 1:n){
                            Rpsp[,y,(n-x)+1] = psp[,x,y]
                        } } }
                if(sample(2,size=1) == 1)  ## equivalent to a coin flip, if 1 then reflect and shift
                    Rpsp[,n:1,n:1] = psp[,ncoord.x,ncoord.y]
                else  ## just shift
                    Rpsp = psp[,ncoord.x,ncoord.y]
            } }
        else{  ## meth is random and sp columns are fixed
            take = sample(n^2)  ## sample w/o replacement
            for(j in 1:S){
                Rpsp[j,,] = matrix(psp[j,,][take],ncol=n,nrow=n)
            } } }
    if(flag)
        Rpsp = drop(Rpsp)
    return(Rpsp)
}

#' Permute an array of occurances under a given set of constraints in
#' 2-dimensions of space with defined spatial strata, see 'nstrata'
#' 
#' @param psp the sp x row x col array, where rows and columns specify where on 
#'   the spatial grid the sample was located
#' @param shiftpos: two numbers that are the x and y places to shift the grid,
#'   this is generated randomly if needed
#' @param rotate: a single number 1-4 that indicates how many counterclockwise 
#'   rotations to perform, generated randomly
#' @param meth: the type of permutation to use, options include: 1) reflect:
#'   random reflection/rotations of species (only makes sence when sp are not
#'   fixed 2) shift: random torodial shifting with or with sp fixed 3) both:
#'   both reflection and shifting 4) random: random shuffle
#' @param sp: if FALSE then obs composition of quadrats is fixed to the observed
#'   pattern if TRUE then species are each shuffled independently
#' @param nstrata: is the number of strata along a single spatial axis within
#'   which to randomize
#' @export
spatPermStrata = function(psp,shiftpos=NULL,rotate=NULL,meth='shift',
                          sp=FALSE,nstrata=1){
    n = dim(psp)[2]
    if(length(dim(psp)) == 3){
        S = dim(psp)[1]
        flag = FALSE
    }
    else{
        S = 1
        psp = array(psp,dim=c(S,n,n))
        flag = TRUE
    }
    strata.size = n/nstrata
    if(round(strata.size) != strata.size){
        stop('Number of strata must be evenly divisable by the linear dimension
             of the grid')
    }
    ## now simply apply the function spatPerm2D on subsets of the orginal matrix
    ## and append all the pieces together at end
    Rpsp = psp
    brks = seq(1,n,strata.size)
    for(i in 1:nstrata){
        sub.i = brks[i]:(brks[i]+strata.size-1)
        for(j in 1:nstrata){
            sub.j = brks[j]:(brks[j]+strata.size-1)
            Rpsp[,sub.i,sub.j] = spatPerm2D(psp[,sub.i,sub.j],meth=meth,sp=sp,
                                            shiftpos=shiftpos,rotate=rotate)
        } }
    if(flag)
        Rpsp = drop(Rpsp)
    return(Rpsp)
}

#' Shuffle a community site x species matrix
#' @param comm site x species matrix with abundance or pres/absen data
#' @param swap: two options: 'indiv' or 'quad' for individual or quadrat-based shuffling
#' @export
#' @examples
#' x = matrix(1:100, 10, 10)
#' shuffle_comm(x, 'quad')
#' shuffle_comm(x, 'indiv')
shuffle_comm = function(comm, swap) {
    if (swap == 'indiv') {
        nquad = nrow(comm)
        comm_shuffled = comm
        for (j in 1:ncol(comm)) {
            rand_samp = sample(nquad, size=sum(comm[ , j]), replace=T)
            comm_shuffled[ , j] = as.numeric(table(c(rand_samp, 1:nquad)) - 1)
        }
    }
    else if (swap == 'quad')
        comm_shuffled = comm[sample(nrow(comm)), ]
    else
        stop('swap must be "indiv" or "quad"')
    return(comm_shuffled)
}

#' Generate individual or sample-based random shuffle statistical null 
#' expectations for variograms.
#' 
#' @param x is either an output of class 'sim' that is the output of
#'   'sim.neut.uni' OR an site x species matrix
#' @param vobject is the output of 'vario', this informs the function of what
#'   parameterization of vario to use specifically it indiates if the pos.neg
#'   components and median should be calculated
#' @param swap two options: 'indiv' or 'quad' for individual or quadrat-based
#'   shuffling repsectively.
#' @param nperm is the number of permutations
#' @param npar is the number of processors to run the analysis over
#' @param coords the spatial coordinates of the sites, not needed if x is of
#'   class 'sim'
#' @param breaks what spatial breaks to use
#' @note The sample-based shuffling the mean of the null distribution is the
#'   same as the average variance across all distance classes. So shuffling is
#'   only useful in deriving the variance around that expectation.
#' @note if x is a presence-absence matrix then the individual and sample-based 
#'   approaches will generate identical within-species variograms
#' @note this function requires snowfall and rlecuyer packages are loaded
#' @export
random_shuffle = function(x, vobject, swap, nperm, npar, coords=NULL, breaks=NA) {
    dists = vobject$vario$Dist
    grain = vobject$parms$grain
    hmin = vobject$parms$hmin
    hmax = vobject$parms$hmax
    pos.neg = vobject$parms$pos.neg
    median = vobject$parms$median
    if (class(vobject$parms$direction) == "factor") 
        direction = as.character(vobject$parms$direction)
    else
        direction = as.numeric(vobject$parms$direction)
    tolerance = vobject$parms$tolerance
    unit.angle = as.character(vobject$parms$unit.angle)
    if (is.na(unit.angle))
        unit.angle = 'degrees'
    distance.metric = as.character(vobject$parms$distance.metric)
    if (is.null(coords))
        stop('need to supply spatial coordinates if not a simulation product')
    r.vals = list()
    r.vals$parms = vobject$parms
    r.vals$p = vobject$p
    if (npar > 1) {
        sfClusterSetupRNG()
        sfSource('./scripts/spat_functions.R')
        sfExport('x', 'swap', 'coords', 'distance.metric')
        if (direction == 'bisection') {
            rv = sfSapply(1:nperm, function(...)
                vario_bisect(shuffle_comm(x, swap), coords,
                             distance.metric=distance.metric)$vario[ , 'exp.var'])
        }  
        else {
            sfExport('grain','breaks','hmin','hmax','pos.neg','median','direction',
                     'tolerance','unit.angle')
            rv = sfLapply(1:nperm, function(...)
                vario(shuffle_comm(x, swap), coords, grain=grain, breaks=breaks, hmin=hmin,
                      hmax=hmax, pos.neg=pos.neg, median=median, direction=direction,
                      tolerance=tolerance, unit.angle=unit.angle,
                      distance.metric=distance.metric)$vario[ , c('exp.var', 'obs.var')])
        }  
    }
    else {
        if (direction == 'bisection')
            rv = sapply(1:nperm, function(...)
                vario_bisect(shuffle_comm(x, swap), coords,
                             distance.metric=distance.metric)$vario[ , 'exp.var'])
        else
            rv = lapply(1:nperm, function(...)
                vario(shuffle_comm(x, swap), coords, grain=grain, breaks=breaks, hmin=hmin,
                      hmax=hmax, pos.neg=pos.neg, median=median, direction=direction,
                      tolerance=tolerance, unit.angle=unit.angle,
                      distance.metric=distance.metric)$vario[ , c('exp.var', 'obs.var')])
    }
    if (is.matrix(rv)) 
        rv_dim = dim(rv)
    else if (is.list(rv))
        rv_dim = dim(rv[[1]])
    r.vals$vario = array(NA, dim= rv_dim + c(0, 1))
    if (direction == 'bisection') {
        r.vals$vario[ , 1] = as.matrix(vobject$vario[ , 'exp.var'])
        r.vals$vario[ , -1] = rv
    }
    else {
        r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var')])
        r.vals$vario[ , , -1] = array(unlist(rv), dim=c(rv_dim, nperm))  
    }
    r.vals$perm = TRUE
    r.vals$vdists = vobject$vario$Dist
    return(r.vals)
}




## Plotting Functions-----------------------------------------------------------
#' Plot community variograms, the results of function 'vario'
#' 
#' @param vobject the output of the function 'vario' or the function 'nullPerm'
#' @param optim is the location of the species niches
#' @param exp.only if TRUE then only the two expected (spatial and nonspatial) 
#'   components of variance displayed
#' @param flip.neg if TRUE then negative fraction is expressed as a positive
#'   value
#' @param ylim the limits on the y-axis
#' @param xlab the label for the x-axis
#' @param ylab the label of the y-axis
#' @param cls the colors for the curves
#' @export
#' @examples
#' x = matrix(runif(100), 10, 10)
#' coord = cbind(1:10, 1)
#' v = vario(x, coord)
#' vGraph(v)
vGraph = function(vobject,optim=NA,exp.only=FALSE,flip.neg=FALSE,ylim=NULL,
                  xlab=NULL,ylab=NULL,cls=NULL) {
    n = sqrt(vobject$parms$N)
    N = n^2
    S = vobject$parms$S
    coords = cbind(rep(1:n,each=n),rep(1:n,times=n))
    v = vobject$vario
    p = vobject$p
    if(is.null(cls)){
        ##green, purple, orange with transparency
        cls = c(rgb(127, 201, 127,alpha=255*.5,maxColorValue = 255),
                rgb(190, 174, 212,alpha=255*.5,maxColorValue = 255), 
                rgb(253, 192, 134,alpha=255*.5,maxColorValue = 255))
    }
    if(is.null(xlab))
        xlab = 'Spatial Lag'
    if(is.null(ylab))
        ylab = 'Variance'
    if(!is.na(optim[1])){
        y = 1:S
        plot(sort(optim),y,axes=F,frame.plot=F,xlab='grad',ylab='',xlim=c(0,max(optim)))
        axis(side=1)
        axis(side=2,at=c(1,S))
        arrows(sort(optim)-vobject$parms$niche.wid.rel/2,y,
               sort(optim)+vobject$parms$niche.wid.rel/2,y,angle=90,len=.01,code=3)
    }
    if(vobject$perm){ ##if vobject is result of permutations
        d = vobject$vdists ##assign distances
        ##for each distance class I want to calculate the upper 95 and lower 95
        ##for both the neg and pos curves
        quants = apply(v,c(1,2),quantile,c(.025,.975))
        if(is.null(ylim))
            ylim=range(list(quants,v[,,1]))
        if(vobject$parms$pos.neg){
            plot(d,v[,1,1],col=1,pch=19,type='o',ylab=ylab,xlab=xlab,ylim=ylim)
            points(d,v[,1,1]-v[,2,1]-v[,3,1],col='green3',type='o',pch=19)
            lines(d,v[,2,1],col='purple',lwd=2)
            lines(d,v[,3,1],col='orange',lwd=2)
        }
        else{
            plot(d,v[,1,1],pch=19,type='o',col='green3',ylab=ylab,xlab=xlab,ylim=ylim) #expected values
            lines(d,v[,2,1],col=1,lwd=2,type='o',pch=19) ##observed variance
        }
        lines(d,rep(sum(p*(1-p)),length(d)),col='blue',lwd=2)
        for(i in 1:dim(quants)[3])
            polygon(c(d,d[length(d):1]),c(quants[1,,i],quants[2,length(d):1,i]),
                    col=cls[i],border = NA)
    }
    else{
        if(is.null(ylim)){
            if(vobject$parms$pos.neg){
                if(flip.neg) v[,8] = v[,8]*-1
                ylim=range(list(v[,c(4:5,7:8)],sum(p*(1-p))))
            }
            else{
                if(exp.only)
                    ylim=range(list(v[,4],sum(p*(1-p))))
                else
                    ylim=range(list(v[,4:5],sum(p*(1-p))))
            }}
        plot(v$Dist,v$obs.var,col=1,pch=19,type='n',ylab=ylab,xlab=xlab,ylim=ylim)
        points(v$Dist,v$exp.var,col='green3',type='o',pch=19)
        if(!exp.only)
            points(v$Dist,v$obs.var,col=1,pch=19,type='o')
        if(vobject$parms$pos.neg){
            lines(v$Dist,v$pos,col='purple',lwd=2)
            lines(v$Dist,v$neg,col='orange',lwd=2)
        }
        lines(v$Dist,rep(sum(p*(1-p)),length(v$Dist)),col='blue',lwd=2)
    }
}

#' Plot community variograms, when permutations have been carried out
#'
#' @param vrand output of nullPerm
#' @param vspat output of vario (?)
#' @param obs.var boolean, if TRUE the obs var is plotted
#' @param ylims
#' @param xlims
#' @param ylab
#' @param xlab 
#' @param cls
#' @param lwd
#' @param plot.new
#' @export
vGraphPerm = function(vrand=NULL,vspat=NULL,obs.var=FALSE,ylims=NA,xlims=NA,
                      ylab='variance',xlab='lag',cls=NA,lwd=1,plot.new=TRUE) {
    rflag = !is.null(vrand)
    sflag = !is.null(vspat)
    if( rflag ){
        dr = vrand$vdists ##assign distances
        vr = vrand$vario[,1,] ##just the expected variance component
    }
    if( sflag){
        ds = vspat$vdists
        vs = vspat$vario
    }
    if(is.na(xlims[1]))
        xlims = range(dr,ds)
    if(is.na(cls)){
        ##purple, blue, green
        cls = c(rgb(190, 174, 212,alpha=255*.5,maxColorValue = 255),
                '#99CCFF',
                rgb(127, 201, 127,alpha=255*.5,maxColorValue = 255))
    }
    if(rflag & sflag & plot.new){
        if(obs.var) 
            par(mfrow=c(1,3))
        else
            par(mfrow=c(1,2))
    }
    if(rflag){
        q.rand = apply(vr,1,quantile,c(.025,.975))
        ylims=range(q.rand,vr[,1])
        plot(dr,vr[,1],type='n',ylab=ylab,xlab=xlab,ylim=ylims,xlim=xlims,
             main='Within-species Agg.')
        polygon(c(dr,dr[length(dr):1]),c(q.rand[1,],q.rand[2,length(dr):1]),
                border=NA,col=cls[1])
        lines(dr,vr[,1],col=1,lwd=lwd)
    }
    if(sflag){
        if(obs.var){
            ovar = vs[,2,] + vs[,3,] ##pos + neg fractions
            q.obs = apply(ovar,1,quantile,c(.025,.975))
            ylims=range(q.obs,ovar)
            plot(ds,ovar[,1],type='n',ylab=ylab,xlab=xlab,ylim=ylims,xlim=xlims,
                 main='Total Between-species Agg.')
            polygon(c(ds,ds[length(ds):1]),c(q.obs[1,],q.obs[2,length(ds):1]),
                    border=NA,col=cls[1])
            lines(ds,ovar[,1],col=1,lwd=lwd)
        }
        q.spat = apply(vs,1:2,quantile,c(.025,.975))
        ylims=range(q.spat,vs[,,1])
        plot(ds,vs[,2,1],type='n',ylab=ylab,xlab=xlab,ylim=ylims,xlim=xlims,
             main='Pos/Neg Between-species Agg.')
        polygon(c(ds,ds[length(ds):1]),c(q.spat[1,,2],q.spat[2,length(ds):1,2]),
                border=NA,col=cls[2])
        polygon(c(ds,ds[length(ds):1]),c(q.spat[1,,3],q.spat[2,length(ds):1,3]),
                border=NA,col=cls[3])
        lines(ds,vs[,2,1],col='dodgerblue',lwd=lwd)
        lines(ds,vs[,3,1],col='green3',lwd=lwd)
    }}

