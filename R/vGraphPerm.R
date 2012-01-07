'vGraphPerm' = function(vrand=NULL,vspat=NULL,obs.var=FALSE,ylims=NA,xlims=NA,
                        ylab='variance',xlab='lag',cls=NA,lwd=1,plot.new=TRUE)
{
 ## Purpose: to graph community variograms, the results of function 'vario'
 ## when all permutations have been run
 ## Arguments
 ## 'vobject' the output of the function 'vario' or the function 'nullPerm'
 ## 'draw' may be 'exp','total', 'both' or 'pos-neg'
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
