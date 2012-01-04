/* 
  AUTHOR: Dan McGlinn
  CONTACT: danmcglinn@gmail.com
  PURPOSE: The purpose of this alogrithim is to implement the random patterns
  algo described in Roxburgh and Chesson (1998).  Towards this end this algo
  accepts an initial observed arrangement of occurances and an initially
  nullfied arrangedment of occurances.  The output is a null arrangement of
  occurances that approximates the observed arrangement in terms of their first
  order spatial statistics. These statistics include: E = # number of neighbors
  in cardinal directions, C = # number of neighbors in intermediate directions
  (i.e., diagonal directions), S = # number of times all sampled neighbors are
  occupied, and  O = # of times all the sampled neighbors are unoccupied. There
  are important differences in the implementation of this algo and that
  described in Roxburgh & Chesson.  The primary difference is that here we do
  not have special counting rules of neighbors if you are along the edges of 
  the landscape. Instead we have artifically extended the landscape one pixel
  in each direction prior to inputing it to this function. This allows the algo
  to search outside of the observed spatial domain by one pixel. My primary
  reason for doing this is that this drops the reqirement of treating edges of
  the observed grid differently from internal cells.  Additionally becuase this
  algo will be given sets of cells in which there are no observations it would
  be quite complex if we started counting the spatial stats specially for each
  place we have no sampled data (i.e., edge of a grid or empty spot). Unsampled
  grid cells have an abundance less than -998, the algorithim takes this into
  consideration in the sub function spatstat.
  REFERENCES:
    Roxburgh, S.H. Chesson, P.A. New Method for Detecting Species Associations
      with Spatially Autocorrelated Data. Ecology. 79: 2180-2192.
  ARGUMENTS:
    psp:     a vector of real numbers that specifies presences/abundances in a 2-D
             spatial sampling grid, that has been modified by adding sites around
             the grid for the purposes of this algo, -998 indicates an unsampled
             pixel
    rpsp:    an intially randomized version of psp 
    n:       number of sites along a single spatial dimension of the n x n grid
    ostat:   the observed spatial statistics, length 4
    nstat:   the null spatial statistics, length 4
    phi:     the ratio of sum(abs(null/obs -1)) for four spatial stats
    phiTemp: a temporary version of the above argument, simply for convience
    alpha:   the cutoff value for the phi statistic of Roxburgh and Chesson 1998
    pl:      the places in rpsp that can be swaped 
    nplaces: single value that is the length of the pl vector  
    ntrials: a counter for the number of times an attempted swap was made 
    gtrials: a counter for the number of times an attempted swap was a 'good'
             swap (i.e., it lowered phi)
    mtrials: the maximum number of times to attempt swapping
  OUTPUT:
    rpsp:    the randomized sitexsp version of psp
    phi:     the goodness of fit statistic for the swapping algo
    gtrials: how many good swaps occurred
  PROGRAM FLOW AND OUTLINE:
    calc the initial spatial stats for observed and null presences
    calc the intial value of phi
    begin swapping algo
      select two random places to swap 
      swap the cells that correspond to those chosen in previous step
      calc null spatial stat
      calc temporary phi value with obs and new null stats
      if temp phi < phi and we haven't tried to swap more than mtrials times
         then keep the swap and update phi
      else 
        change the swap back and do not update phi
*/
 
#include <R.h>
#include <Rmath.h> /*contains the runif function*/
#include <stdlib.h> /* necessary to use the function frand() */
#include <math.h> /* necessary for fabs which gives abs of a double */

void randpatpar 
  (double *psp, double *rpsp, int *n, double *ostat, double *nstat, 
   double *phi, double *phiTemp, double *alpha, int *pl, int *nplaces, 
   double *ntrials, double *gtrials, double *mtrials)
{ 
  void spatstat(double *psp, int *n, double *outstat) ;
  void calcphi(double *nstat, double *ostat, double *phi) ;
  GetRNGstate() ;  /* Initialize the R random number generator */
  int size, rint1, rint2, swap1, swap2, k ;
  size = (*n+2)*(*n+2) ; 
  double rpspTemp[size], rdouble ;
  spatstat (psp, n, ostat) ;
  spatstat (rpsp, n, nstat) ;
  calcphi (nstat, ostat, phi) ;
  /* now begin random patts algo */      
  while ((*phi > *alpha) && (*ntrials < *mtrials)){ 
    /*choose random integer using R function runif */
    rdouble = runif(0,*nplaces) ;  /*note that we want integers from 0 to nplaces-1 */
    rint1 = floor(rdouble) ;
    rint2 = rint1 ; 
    while(rint2 == rint1){
      rdouble = runif(0,*nplaces) ; 
      rint2 = floor(rdouble) ;
    }
    if(rint1 < rint2){
      swap1 = pl[rint1] ;
      swap2 = pl[rint2] ;
    }
    else{
      swap1 = pl[rint2] ;
      swap2 = pl[rint1] ;
    }
    /* create a temporary holding vector for the random array */
    for(k = 0 ; k < size ; k++){ 
      rpspTemp[k] = rpsp[k] ;
    }
    /*carry out the swaps */
    rpsp[swap1] = rpspTemp[swap2] ;
    rpsp[swap2] = rpspTemp[swap1] ; 
    spatstat (rpsp, n, nstat) ;
    calcphi (nstat, ostat, phiTemp) ;
    if(*phiTemp < *phi){  /*then keep the swap */
      *phi = *phiTemp ;      
      *gtrials = *gtrials + 1 ;
    }
    else{  /*swap back to orginal orientation*/
      rpsp[swap1] = rpspTemp[swap1] ;
      rpsp[swap2] = rpspTemp[swap2] ;  
    }
    *ntrials = *ntrials + 1 ; 
  }
  /*calc spatstat one last time to return the best null stats*/
  spatstat(rpsp, n, nstat) ;  
  PutRNGstate() ;  /* cleanup the R random number generator */
}

void spatstat (double *psp, int *n, double *outstat) 
{
  int j, col, row, Epot, Cpot ;  
  double Etemp, Ctemp ; 
  for(j = 0 ; j < 4 ; j++){
    outstat[j] = 0 ;
  }
  for(row = 1 ; row < (*n+1) ; row++){
    for(col = 1 ; col < (*n+1) ; col++){
      if (psp[col+(row*(*n+2))] > 0) { /* if cell occupied */
        Etemp = 0 ;
        Ctemp = 0 ;
        Epot = 4 ;
        Cpot = 4 ;
        /* calculate number of edge contacts */
        if(psp[(col+1)+(row*(*n+2))] < -998)
          Epot-- ;
        else
          Etemp += psp[(col+1)+(row*(*n+2))] ; 
        if(psp[(col-1)+(row*(*n+2))] < -998)
          Epot-- ;
        else
          Etemp += psp[(col-1)+(row*(*n+2))] ;
        if(psp[col+((row+1)*(*n+2))] < -998)
          Epot-- ;
        else
          Etemp += psp[col+((row+1)*(*n+2))] ;
        if(psp[col+((row-1)*(*n+2))] < -998)
          Epot-- ; 
        else
          Etemp += psp[col+((row-1)*(*n+2))] ;
        outstat[0] += Etemp ;
        /* calculate number of corner contacts */
        if(psp[(col-1)+((row-1)*(*n+2))] < -998)
          Cpot-- ;
        else
          Ctemp += psp[(col-1)+((row-1)*(*n+2))] ;
        if(psp[(col-1)+((row+1)*(*n+2))] < -998)
          Cpot-- ;
        else
          Ctemp += psp[(col-1)+((row+1)*(*n+2))] ;
        if(psp[(col+1)+((row-1)*(*n+2))] < -998)
          Cpot-- ;
        else
          Ctemp += psp[(col+1)+((row-1)*(*n+2))] ;
        if(psp[(col+1)+((row+1)*(*n+2))] < -998)
          Cpot-- ;
        else
          Ctemp += psp[(col+1)+((row+1)*(*n+2))] ;
        outstat[1] += Ctemp  ;
        if((Epot + Cpot) > 0){
          outstat[2] += (Etemp + Ctemp)/(Epot + Cpot) ;
          if((Etemp + Ctemp) == 0) /* if all neighbors empty add one to open */
            outstat[3] += psp[col+(row*(*n+2))] ;
        }
      }  
    }
  }
}

void calcphi (double *nstat, double *ostat, double *phi)
{
  int j ;
  *phi = 0 ;
  for (j = 0 ; j < 4 ; j++){
    if(ostat[j]>0){
      *phi += fabs(nstat[j] / ostat[j] - 1) ;
      /*Rprintf("phi is %f\n",*phi);*/
    }
  }
}
