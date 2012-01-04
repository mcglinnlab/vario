/*
  Code Citation:  
    Paulo J. Ribeiro Jr & Peter J. Diggle geoR: a package for
      geostatistical analysis R-NEWS, 1(2):15-18. June, 2001
*/

#include <stdio.h>
#include <R.h>

#define Integer int
#define Real double


void tgangle(Real *xloc, Real *yloc, Integer *nl, Real *res) 
     /* 
	This function computes the tangent of the (azimuth) 
        angle between pairs of locations
	
        xloc, yloc     : xy coordinates of the locations
	nl,            : number of locations
	res            : stores the results to be returned, 
	                 a vector with tangent of the angles
     */   
     
{ 
  Integer i,j, ind;
  Real dx,dy;
  ind = 0;
  for (j=0; j<*nl; j++) {  
    for (i=j+1; i<*nl; i++) {
      dx = (xloc[i] - xloc[j]) ;
      dy = (yloc[i] - yloc[j]) ;
      /* Care here: considering the AZIMUTH angle, therefore dx/dy */
      res[ind] = dx/dy ;
      ind++ ;
    }
  }
}
