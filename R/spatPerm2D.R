
'spatPerm2D' = function(psp,shiftpos=NULL,rotate=NULL,meth='shift',sp=FALSE)
{
  ## Purpose: to permute an array of occurances under a given set of constraints
  ## in 2-dimensions of space
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
