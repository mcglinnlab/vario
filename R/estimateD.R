'estimateD' = function(z)
{
 ##esimates the fractal dimension D from a 2-dimensional grid
 ##z is a matrix of real numbers
 n = dim(z)[1]
 gcords = expand.grid(1:n,1:n)
 v = vario(as.vector(z),gcords)$vario
 mod = lm(log(v$exp)~log(v$Dist))
 m = coef(mod)[2]
 D = (6-m)/2
 return(D)
}