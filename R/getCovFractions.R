'getCovFractions' = function(x)
{
  ## Purpose: calculates the lower diagonal of a sp covariance matrix
  ## to provide the positive and negative fractions of covariance
  ## output is two lower triangular matrices, each in vector format
  ## Called within the function 'vario'
  ## Arguments: 
  ## x is a sitexsp matrix (sp as columns) of real numbers
  ## rows are the sites, columns are the species
  N = as.integer(nrow(x))
  S = as.integer(ncol(x)) 
  x = as.double(ifelse(is.na(x) | x == -999,-99999,x))
  pos = as.double(rep(0,(N*(N-1))/2))
  neg = as.double(rep(0,(N*(N-1))/2))
  result = .C('loopcovreal',x,N,S,pos,neg,PACKAGE = vario)
  out = list()
  out$pos = result[[4]]
  out$neg = result[[5]]
  return(out)
} 