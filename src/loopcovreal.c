/* x is a sitexsp matrix (sp as columns) of real numbers
   n is number of sites (i.e. rows)
   s is number of species (i.e. columns)
   pos and neg are lists of zeros the length of the unique pairwise combinations of species and site 
*/

#include <R.h>

void loopcovreal(double *x, int *n, int *s, double *pos, double *neg){
 int i, j, k, l, icount ; 
 double d ;
 icount = 0 ; 
 for (i = 0; i < (*n-1); i++) { /*loop through sites */
  while(x[i] < -99998) /* assuming that if the first species is -99999 (i.e. unsampled) that all the rest are as well */
   i++ ;
  for (j = (i+1); j < *n; j++) {  
   while(x[j] < -99998)
    j++ ;
   for (k = 0; k < (*s-1); k++) { /*loop through species */
    for (l = (k+1); l < *s; l++) {
     d = (x[i+(k**n)]-x[j+(k**n)]) * (x[i+(l**n)]-x[j+(l**n)]) ;
     if (d > 0) {
      pos[icount] += d ;
     } 
     else { 
      neg[icount] += d ;
   }}}
   icount++ ;
}}}
