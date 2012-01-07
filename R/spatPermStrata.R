'spatPermStrata' = function(psp,shiftpos=NULL,rotate=NULL,meth='shift',
                            sp=FALSE,nstrata=1)
{
  ## Purpose: to permute an array of occurances under a given set of constraints
  ## in 2-dimensions of spacewith defined spatial strata, see 'nstrata'
  ## argument below
  ## Arguments:
  ## psp: the sp x row x col array, where rows and columns specify where on
  ##      the spatial grid the sample was located
  ## shiftpos: two numbers that are the x and y places to shift the grid, this
  ##           is generated randomly if needed
  ## rotate: a single number 1-4 that indicates how many counterclockwise
  ##         rotations to perform, generated randomly
  ## meth: the type of permutation to use, options include:
  ##   reflect: random reflection/rotations of species (only makes sence when
  ##              sp are not fixed
  ##   shift: random torodial shifting with or with sp fixed
  ##   both: both reflection and shifting
  ##   random: random shuffle
  ## sp: if FALSE then obs composition of quadrats is fixed to the observed pattern
  ##     if TRUE then species are each shuffled independently
  ## nstrata: is the number of strata along a single spatial axis within which
  ##          to randomize
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