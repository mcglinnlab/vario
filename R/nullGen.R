'nullGen' = function(pop,vobject,coords,meth,sp,all=FALSE,RPargs=FALSE,median=FALSE)
{
 ## Purpose: to mediate the generation of statistical null values for the variograms 
 ## to be used in a parrallel processing loop which will generate a population
 ## of null values
 ## Arguments:
 ## pop: the species x row x col array where row and column refer to spatial location
 ## vobject: the output of the vario function with serves as the basis for 
 ##          empirical comparisons
 ## coords: the spatial coordinates of a M^2 x S matrix
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
 ## all: boolean, if TRUE then all relevant nulls are calculated
 ## RPargs: a list of arguments that must be supplied if the random patterns
 ##         algo is desired, the arguments of RPargs are input into the
 ##         function 'randPatPar', they include:
 ##         1)allRP if TRUE & all = TRUE, then Random Patterns algo used as the
 ##         spatial null, 2)nstrata, 3)mtrials1, 4)mtrials2, 5)alpha, 6)npar
 ## median: if TRUE then means and medians are calculated
 ## Note: 'meth' and 'sp' are arguments to randomization function 'spatPerm2D'
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
   rv = vario(x=rmat,coord=coords,grain=grain,hmax=hmax,pos.neg=pos.neg,
              median=median,direction=direction,tolerance=tolerance,
              unit.angle=unit.angle)$vario
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
  rv = vario(x=rmat,coord=coords,grain=grain,hmax=hmax,pos.neg=pos.neg,median=median,
            direction=direction,tolerance=tolerance,unit.angle=unit.angle)$vario
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