'fixUnSampFalseBorder' = function(oarray,rarray)
{
  ## Purpose: to maintain the spatial locations of the unsampled pixels in
  ## rarray which is a random realization of oarray, -999 is the identifier for 
  ## unsampled cells, in this case oarray and rarray DO have a false border of -999
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