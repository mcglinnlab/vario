/*
  AUTHOR: Dan McGlinn
  CONTACT: danmcglinn@gmail.com
  PURPOSE: To compute the positive and negative fractions of covariance between
  the columns (i.e., species) of a matrix   
  ARUGMENTS: 
    x:   one dimensional representation of a sitexsp matrix (sp as columns) 
         of real numbers, unsampled rows (i.e. sites) have value -99999
    n:   number of sites (i.e. rows)
    s:   number of species (i.e. columns)
    pos: list of zeros the length of the unique pairwise combinations of
         species and site 
    neg: list of zeros the length of the unique pairwise combinations of
         species and site 
  OUTPUT: 
    pos: vector of positive covariances
    neg: vector of negative covariances
*/

#include <R.h>

void loopcovreal(double *x, int *n, int *s, double *pos, double *neg){
  int i, j, k, l, icount ; 
  double d ;
  icount = 0 ; 
  for (i = 0 ; i < (*n-1); i++) { /*loop through sites */
    while(x[i] < -99998) 
      i++ ;
    /* the above while loop assumes that if the first species is -99999
      (i.e. unsampled) that all the rest of the species are unsampled as well 
    */
    for (j = (i+1) ; j < *n ; j++) {  
      while(x[j] < -99998)
        j++ ;
      for (k = 0 ; k < (*s-1) ; k++) { /*loop through species */
        for (l = (k+1) ; l < *s ; l++) {
          d = (x[i+(k**n)]-x[j+(k**n)]) * (x[i+(l**n)]-x[j+(l**n)]) ;
          if (d > 0) {
            pos[icount] += d ;
          } 
          else { 
            neg[icount] += d ;
          }
        }  
      }   
    icount++ ;
    }
  }
} 