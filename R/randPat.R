'randPat' = function(i,psp,rpsp,n,nstrata,pl,mtrials1=1e3,mtrials2=1e6,
                     alpha=0.01)
{
  ## Purpose: to be called in serial or parallel by function "randPatPar'
  ## this function evaulates the .C function 'randpatpar' which is the random patterns algo of
  ## Roxburgh and Chesson 1998. Returns species index, phi stat, number of actual swaps, and the 
  ## randomized presences as a single vector of numbers
  ## Arguments:
  ## i: the ith species index
  ## psp: multidimenstional S x (n+2) x (n+2) array
  ## rpsp: a randomized version of psp
  ## n: the size of the orginal 2-D array along one spatial axis (i.e., without extra rows and columns)
  ## pl: the places in rpsp that can be swaped
  ## mtrials1: the number of times to attempt a swap at the strata level
  ## mtrials2: the number of times to attempt a swap at the pixel level
  ## alpha: the cutoff value for the phi statistic of Roxburgh and Chesson 1998
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
             as.double(rep(0,4)),PACKAGE = vario)[[3]]
  nstat = .C('spatstat',as.double(as.vector(rpsp)),as.integer(n),
             as.double(rep(0,4)),PACKAGE = vario)[[3]]
  phi = .C('calcphi',as.double(nstat),as.double(ostat),as.double(0),
           PACKAGE = vario)[[3]]
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
               as.double(rep(0,4)),PACKAGE = vario)[[3]]
    phiTemp = .C('calcphi',as.double(nstat),as.double(ostat),as.double(0),
                 PACKAGE = vario)[[3]]
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
             PACKAGE = vario)
    out = c(i,phi,gtrials,tmp$phi,tmp$gtrials,tmp$rpsp)
  }
  else
    out = c(i,phi,gtrials,NA,NA,as.vector(rpsp))
  return(out)
}