'randPatPar' = function(psp,nstrata,mtrials1=1e3,mtrials2=1e6,alpha=0.01,npar=1)
{
 ## Purpose: convience function for working with randPat which calls the .C
 ## function 'randpatpar'. This function allows you to specify the number of
 ## processors to run on adding processsors only helps if working with many
 ## species as each species is evaulated on a different processor. Returns a 
 ## (5+(n+2)^2) x S matrix, the first five rows of which are species index, phi
 ## strata stat, number of strata swaps, phi pixel swap, and number of pixel
 ## swaps, and then the remaining rows are the presences/abundances in the
 ## randomized occurances
 ## Arguments:
 ## psp: multidimenstional S x (n+2) x (n+2) array
 ## n: the size of the orginal 2-D array along one spatial axis (i.e., without
 ##    extra rows and columns)
 ## pl: the places in rpsp that can be swaped
 ## mtrials: the numbef of times to attempt a swap
 ## alpha: the cutoff value for the phi statistic of Roxburgh and Chesson 1998
 ## npar: the number of processors to run the code on
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