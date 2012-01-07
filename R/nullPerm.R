'nullPerm' = function(x,vobject,nperm,coords=NULL,meth='both',sp=TRUE,all=FALSE,
                       snap=NULL,npar=1,RPargs=FALSE)
{
 ## Purpose: to generate statistical null expectations for the variograms
 ## Arguments:
 ## x: is either an output of class 'sim' that is the output of 'sim.neut.uni'
 ##    OR an site x species matrix
 ## vobject: is the output of 'vario', this informs the function of what
 ##          parameterization of vario to use specifically it indiates if the
 ##          pos.neg components and median should be calculated
 ## nperm: is the number of permutations
 ## coords: the spatial coordinates of the sites, not needed if x is of class 'sim'
 ## meth: the type of permutation to use, options include:
 ##       reflect: random reflection/rotations of species (only makes sence
 ##                when sp are not fixed
 ##       shift: random torodial shifting with or with sp fixed
 ##       both: both reflection and shifting
 ##       random: random shuffle
 ##       randpat: random patterns algo of Roxburgh and Chesson 1998, must
 ##                parameterize RPargs (See below)
 ## sp: boolean, if FALSE then obs composition of quadrats is fixed to the
 ##     observed pattern, if TRUE then species are each shuffled independently
 ## all: boolean, if TRUE then all relevant nulls calculated
 ## snap: is the time period of the simulation to analyze,defaults to NULL val,
 ##       if not set gets internally set to last time period
 ## npar: number of processors to run the function on
 ## median: if TRUE then median is also calculated in addition to mean
 ## RPargs: is a vector of arguments that are needed to perform the Random
 ##         Patterns spatial null model
 ## Notes: 
 ## 1) see the notes associated with the function 'nullGen' that indicate how 
 ##    'RPargs" should be parameterized
 ## 2) 'meth' and 'sp' are arguments to randomization function 'spatPerm2D'
 dists = vobject$vario$Dist
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
 hmax = vobject$parms$hmax
 if(class(x) == 'sim'){
  coords = x$coords
  if(is.null(snap)) snap = length(x$snaps)
 }
 else{
  if(is.null(coords))
   stop('need to supply spatial coordinates if not a simulation product')
 }
 r.vals = list()
 r.vals$parms = vobject$parms
 r.vals$p = vobject$p
 if(median&!pos.neg)
  stop("if computing medians must also compute pos.neg fractions, set pos.neg=TRUE")
 if(pos.neg){ ##pos and neg fractions
  if(all){ ##all relevant null models
   if(median){ ##will compute mean and median
    r.vals$vario = array(0,dim=c(length(dists),6,3,nperm+1))##dists,results,methods,perms
    r.vals$vario[,,,1] = as.matrix(vobject$vario[,c(5,7:11)])   
   }
   else{ ##only compute means
    r.vals$vario = array(0,dim=c(length(dists),3,3,nperm+1))##dists,results,methods,perms
    r.vals$vario[,,,1] = as.matrix(vobject$vario[,c(5,7:8)])
  }}
  else{ ##only a single null used
   if(median){
    r.vals$vario = array(0,dim=c(length(dists),6,nperm+1))##dists,results,perms
    r.vals$vario[,,1] = as.matrix(vobject$vario[,c(5,7:11)])
   }
   else{
    r.vals$vario = array(0,dim=c(length(dists),3,nperm+1))##dists,results,perms
    r.vals$vario[,,1] = as.matrix(vobject$vario[,c(5,7:8)])
 }}}
 else{##only exp and obs fractions
  if(all){
   if(median){ ##will compute mean and median
    r.vals$vario = array(0,dim=c(length(dists),4,3,nperm+1))
    r.vals$vario[,,,1] = as.matrix(vobject$v[,c(4:5,7:8)])
   }
   else{
    r.vals$vario = array(0,dim=c(length(dists),2,3,nperm+1))
    r.vals$vario[,,,1] = as.matrix(vobject$v[,4:5])
  }}
  else{ ##only a single null used
   if(median){
    r.vals$vario = array(0,dim=c(length(dists),4,nperm+1))
    r.vals$vario[,,1] = as.matrix(vobject$v[,c(4:5,7:8)])
   }
   else{ 
    r.vals$vario = array(0,dim=c(length(dists),2,nperm+1))
    r.vals$vario[,,1] = as.matrix(vobject$v[,4:5])
 }}}
 if(class(x) == 'sim'){
  pop = as.logical(x$snaps[[snap]]) ##converts it to a pres/absence vector
  dim(pop) = c(x$p$S, x$p$M, x$p$M)
 }
 else{
  pop = array(x,dim=c(sqrt(nrow(x)),sqrt(nrow(x)),ncol(x)))
  pop = aperm(pop,c(3,1,2))
 }
 if(RPargs[[1]]&npar == 1){
  r.vals$p.conv1 = 0 ##average proportion of species that converged with strata swaps
  r.vals$p.conv2 = 0 ##average proportion of species that converged with pixel swaps
 }
 if(npar == 1){ ##all permutations option not yet implemented for 1 processor
  for(i in 1:nperm){
   if(RPargs[[1]]){ ##use the random pattern algo for the spatial null
    out = randPatPar(psp=pop,nstrata=RPargs[[2]],mtrials1=RPargs[[3]],
                     mtrials2=RPargs[[4]],alpha=RPargs[[5]],npar=RPargs[[6]])
    S = dim(pop)[1]
    n2 = dim(pop)[2]+2
    r.vals$p.conv1 = r.vals$p.conv1 + sum(out[2,]<=RPargs[[5]])/S/nperm
    r.vals$p.conv2 = r.vals$p.conv2 + sum(out[4,]<=RPargs[[5]],na.rm=TRUE)/S/nperm
    rpop = array(0,dim=c(S,n2,n2))
    for(k in 1:S){
     rpop[k,,] = array(out[-(1:5),k],dim=c(n2,n2))
    }
    rpop = rpop[,-c(1,n2),-c(1,n2)]
   }
   else{
    rpop = spatPerm2D(pop,meth=meth,sp=sp)
    rpop = fixUnSampTrueBorder(pop,rpop)
   }
   rmat = apply(rpop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
   rv = vario(x=rmat,coord=coords,grain=grain,hmax=hmax,pos.neg=pos.neg,
              median=median,direction=direction,tolerance=tolerance,
              unit.angle=unit.angle)$vario
   if(pos.neg){
    if(all){
     if(median)
      r.vals$vario[,,,i+1] = as.matrix(rv[,,c(5,7:11)])
     else
      r.vals$vario[,,,i+1] = as.matrix(rv[,,c(5,7:8)])
    } 
    else{
     if(median)
      r.vals$vario[,,i+1] = as.matrix(rv[,c(5,7:11)])
     else
      r.vals$vario[,,i+1] = as.matrix(rv[,c(5,7:8)])
   }}
   else{
    if(all){
     if(median)
      r.vals$vario[,,,i+1] = as.matrix(rv[,,c(4:5,7:8)])
     else
       r.vals$vario[,,,i+1] = as.matrix(rv[,,4:5])
    }
    else{
     if(median)
      r.vals$vario[,,i+1] = as.matrix(rv[,c(4:5,7:8)])
     else
      r.vals$vario[,,i+1] = as.matrix(rv[,4:5])
   }}
   print(i)
  }
 }
 else{ ##computing in parallel
  require(snowfall)
  sfInit(parallel=TRUE, cpus=npar, type="SOCK")
  sfClusterSetupRNG()
  sfExport('pop', 'vobject', 'coords', 'meth', 'all', 'sp', 'RPargs',
           'randPatPar', 'randPat', 'fixUnSampFalseBorder','fixUnSampTrueBorder',
           'spatPerm2D', 'spatPermStrata', 'vario', 'getCovFractions', 'nullGen')
  sfLibrary(vario)
  out = unlist(sfLapply(1:nperm,function(...) nullGen(pop=pop,vobject=vobject,
               coords=coords,meth=meth,sp=sp,all=all,RPargs=RPargs)))
  sfStop()
  if(pos.neg){
   if(all){
    if(median){
     dim(out) = c(length(dists),6,3,nperm)
     r.vals$vario[,,,-1] = out 
    }
    else{
     dim(out) = c(length(dists),3,3,nperm)
     r.vals$vario[,,,-1] = out 
   }}
   else{
    if(median){
     dim(out) = c(length(dists),6,nperm)  
     r.vals$vario[,,-1] = out     
    }
    else{
     dim(out) = c(length(dists),3,nperm)  
     r.vals$vario[,,-1] = out     
  }}}
  else{
   if(all){
    if(median){
     dim(out) = c(length(dists),4,3,nperm)  
     r.vals$vario[,,,-1] = out 
    }
    else{
     dim(out) = c(length(dists),2,3,nperm)  
     r.vals$vario[,,,-1] = out 
   }}
   else{
    if(median){
     dim(out) = c(length(dists),4,nperm)  
     r.vals$vario[,,-1] = out 
    }
    else{
     dim(out) = c(length(dists),2,nperm)  
     r.vals$vario[,,-1] = out 
 }}}}
 r.vals$perm = TRUE
 r.vals$vdists = vobject$vario$Dist
 return(r.vals)
} 