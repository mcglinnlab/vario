'vario' = function(x,coord,grain=1,breaks=NA,hmax=NA,round.up=FALSE,pos.neg=FALSE,
                   binary=TRUE,snap=NA,median=FALSE,direction = 'omnidirectional',
                   tolerance = pi/8, unit.angle = c('radians', 'degrees'),
                   distance.metric = 'euclidean')
{
  ## Purpose: calculates uni- and multi-variate variograms
  ##
  ## This code is largely modified from the 'vegan' function 'mso' by Helene
  ## Wagner, also parts of this code were derived from the 'geoR' function 'vairog'
  ## by Paulo J. Ribeiro Jr. and Peter J. Diggle. 
  ## Citations:
  ## Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre,
  ##   Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M.
  ##   Henry H. Stevens and Helene Wagner (2011). vegan: Community Ecology
  ##   Package. R package version 2.0-2.
  ##   http://CRAN.R-project.org/package=vegan
  ## Paulo J. Ribeiro Jr & Peter J. Diggle geoR: a package for
  ##   geostatistical analysis R-NEWS, 1(2):15-18. June, 2001
  ##
  ## Note: that if some areas are unsampled that obs != pos + neg + exp
  ## Arguments:
  ## x: an object of class 'sim' or may be a sitexsp matrix, a vector of values,
  ##    if x is the later then missing samples should be coded as -999
  ## coord: the spatial coordinates
  ## grain: interval size for distance classes, only used if 'breaks' not supplied
  ## breaks: the spatial breaks that define the spatial lags to be compared
  ## hmax: the maximum spatial lag
  ## round.up: if TRUE the spatial lags are rounded up to nearest integer
  ## pos.neg: if TRUE the positive and negative parts of the covariance matrix
  ##          are output
  ## binary: if TRUE and x is class sim then abundances are converted to
  ##         binary 0 1 values
  ## snap: indicates which generation from an object class 'sim' to draw from
  ## median: indicates if in addition to the mean the medians of the distance
  ##         matrices are calculated
  ## direction: a numerical value for the directional (azimuth) angle. This
  ##            used to specify directional variograms. Default defines the
  ##            omnidirectional variogram. The value must be in the interval
  ##            [0, pi] radians ([0, 180] degrees).
  ## tolerance: numerical value for the tolerance angle, when computing
  ##            directional variograms. The value must be in the interval [0,
  ##            pi/2] radians ([0, 90] degrees).  Defaults to pi/8.
  ## unit.angle: defines the unit for the specification of angles in the two
  ##             previous arguments. Options are 'radians' and
  ##             'degrees', with default to 'radians'.
  ## distance.metric': can be one of the speices turnover metrics listed by the 
  ##                   vegan function vegdist(). This is only appropriate if
  ##                   pos.neg = FALSE. Common options include, 'jaccard' and
  ##                   'bray'. If computed on pres/abse data then soreson index
  ##                   is computed by 'bray'. 
  if(distance.metric != 'euclidean'){
    if(pos.neg)
      stop("cannot commpute pos-neg covariance using a turnover metric")
  }
  else
    require(vegan)
  ## geoR code from function variog starts here'
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
  } }
  ## geoR code from function variog ends here'
  if(class(x) == "sim"){
    coord = x$coords
    if(is.na(snap))  
      snap = length(sim$snaps)
  }
  else
    x = ifelse(x == -999,NA,x)
  Dist = dist(coord)
  if(is.na(breaks[1])){
    if (round.up) 
      H = ceiling(Dist/grain) * grain
    else 
      H = round(Dist/grain) * grain
  }
  else{
    H = Dist
    for(i in 1:(length(breaks)-1))
      H[H>=breaks[i] & H<breaks[i+1]] = breaks[i]
  }
  if (is.na(hmax))
    hmax = round((max(Dist)/2)/grain) * grain  
  vobject = list()
  if(class(x) == "sim"){
    if(binary)
      resp = apply(census(x,snap=snap),c(1,2),as.logical)*1
    else
     resp = census(x,snap=snap)
    vobject$parms = data.frame(grain,hmax,S=ncol(resp),N=nrow(resp),pos.neg,
                               median,niche.wid.rel=x$p$s.rel,
                               disp.wid.rel=x$p$u.rel, direction,tolerance,
                               unit.angle,distance.metric)
  } 
  else{
    resp = as.matrix(x)
    vobject$parms = data.frame(grain,hmax,S=ncol(resp),N=nrow(resp),pos.neg,
                               median,direction,tolerance,unit.angle,
                               distance.metric)
  }
  H[H > hmax] = NA
  H = as.vector(H)
  ## geoR code from function variog with slight modifications starts here'
  if (direction != "omnidirectional") {
    ## note that the changes: 'u' has been changed for as.vector(Dist) and 
    ## coords changed to coord
    u.ang = .C("tgangle", as.double(as.vector(coord[, 1])), 
              as.double(as.vector(coord[, 2])), as.integer(dim(coord)[1]), 
              res = as.double(rep(0, length(as.vector(Dist)))),PACKAGE = vario)$res
    if (any(is.na(u.ang))) 
      stop("NA returned in angle calculations maybe due to co-located data")
    u.ang = atan(u.ang)
    u.ang[u.ang < 0] = u.ang[u.ang < 0] + pi
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
  ## geoR code from function variog ends here'
  Dist = sapply(split(Dist, H), mean,na.rm=TRUE)
  vobject$vario = data.frame(H = names(table(H)), Dist = Dist,
                             n = as.numeric(table(H)))
  ## below 'exp.gamma' is the expected variogram if 'resp' is a sitexsp
  ## pres/abse matrix. The expectation is based upon the assumption of zero
  ## sp x sp covariances
  if(distance.metric == 'euclidean')
    exp.split = split(dist(resp)^2*.5,H)
  else
    exp.split = split(vegdist(resp,method=distance.metric),H) 
  exp.gamma = sapply(exp.split,mean,na.rm=TRUE) 
  vobject$vario = cbind(vobject$vario, exp.var = exp.gamma)
  if(median)
    exp.med = sapply(exp.split,median,na.rm=TRUE) 
  if(ncol(resp) > 1){  ## i.e. resp is a site x sp matrix and not simply a vector
    ## if 'resp' is a sitexsp pres/abse matrix the following computes site species richness
    rich = apply(resp,1,sum) 
    ## see equation 7 in Wagner, H. 2003. Spatial covariance in plant
    ## communities... Ecology 84:1045-1057 to see that observed multivariate
    ## variogram can be computed from the species richness vector
    obs.gamma = sapply(split(dist(rich)^2*.5,H),mean,na.rm=TRUE) 
    if(pos.neg){
      cov.mat = getCovFractions(resp)
      pos.split = split(cov.mat$pos,H)
      neg.split = split(cov.mat$neg,H)
      pos = sapply(pos.split,mean)
      neg = sapply(neg.split,mean)
      if(median){
        pos.med = sapply(pos.split,median)
        neg.med = sapply(neg.split,median) 
        vobject$vario = cbind(vobject$vario, obs.var = obs.gamma, 
                              ratio = obs.gamma/exp.gamma, pos = pos, neg = neg,
                              exp.med = exp.med, pos.med = pos.med,
                              neg.med = neg.med)
      }
      else
        vobject$vario = cbind(vobject$vario, obs.var = obs.gamma,
                              ratio = obs.gamma/exp.gamma, pos = pos, neg = neg)
    ##Note: obs.var = exp.var + pos.var + neg.var
    }
    else{
      if(median){
        obs.gamma.med = sapply(split(dist(rich)^2*.5,H),median,na.rm=TRUE) 
        vobject$vario = cbind(vobject$vario, obs.var = obs.gamma, 
                              ratio = obs.gamma/exp.gamma, exp.med = exp.med,
                              obs.med = obs.gamma.med)
      }
      else
        vobject$vario = cbind(vobject$vario, obs.var = obs.gamma,
                              ratio = obs.gamma/exp.gamma)
  } } 
  vobject$p = apply(resp,2,sum,na.rm=TRUE)/nrow(resp)
  vobject$perm = FALSE 
  ## the above line indicates to the function 'vGraph' if the variogram is the
  ## result of a randomization or not  
  return(vobject)
}