'fixUnSampTrueBorder' = function(oarray,rarray)
{
  ## purpose: to maintain the spatial locations
  ## of the unsampled pixels in rarray which is a random
  ## realization of oarray, -999 is the identifier for 
  ## unsampled cells, in this case oarray and rarray DO NOT have a false border of -999
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