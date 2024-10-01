#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"

#define EULER 0.57721566
#define MAXIT 100
#define FPMIN 1.0e-30
#define EPS 6.0e-8     

void ju_ei(double *x, double *eix){

  /* Computes the exponential intergral Ei(x) for x>0 (hauptwert)
  // AND x<0 (Riemann-Integral)
  // for small |x| it uses the powerseries:
  // Ei(x) = C + ln|x| + x + x^2/2*2! + x^3/3*3!
  // for large x (x > |lnEPS|, EPS required relative error) uses 
  // asymptotic expansion:
  // Ei(x) ~ exp(x)/x*(1 + 1!/x + 2!/x^2 + ...)
  // for large neg x (|x| > 0.5*|lnEPS|)

  // ju-changed MAXIT from 100 => 1000, EPS from 6e-8 => 6e-7

  // Modified NumRecipes function "ei.c"
  // Julia Adolphs, 28.05.04             */

  void nrerror(char error_text[]);

  int k;
  double fact,prev,sum,term;
  

  /* --- x < 0 -------------------------------------------------------------- 
  // Hier existiert das Riemann-Integral, aber keine Stammfunktion, daher
  // Reihenentwicklung */
  
  if (*x < 0.0){ 

    if (*x > -FPMIN){
      *eix = log(fabs(*x))+EULER;
    }  

    if( *x >= 0.5*log(EPS) ){
      sum = 0.0;
      fact = 1.0;
      for (k=1; k<=MAXIT; k++) {
	fact *= (*x)/k;
	term = fact/k;
	sum += term; 
	if ( fabs(term) < fabs(EPS*sum) ) 
	  break;

      }
      if (k > MAXIT) 
	nrerror("Series failed in ei");
      
      *eix = sum+log(fabs(*x))+EULER;
    } 
    else {
      sum = 0.0;
      term = 1.0;
      for (k=1; k<=MAXIT; k++) {
	prev = term;
	term *= k/(*x);
	if ( fabs(term) < EPS ) 
	  break;
	if ( fabs(term) < fabs(prev) ) 
	  sum += term;
	else {
	  sum -= prev;
	  break;
	}
      }    
      *eix = exp(*x)*(1.0+sum)/(*x);
    }
    
  } /* endif x<0


  -----  x > 0  (Cauchy Hauptwertintegral)  ----------------------------- */

  if (*x > 0.0){ 

    if (*x < FPMIN){         /* special case: invoid failure on convergence */
      *eix = log(*x)+EULER;
    }  
    

    if( *x <= -log(EPS) ){
      sum = 0.0;
      fact = 1.0;
      for (k=1; k<=MAXIT; k++) {
	fact *= (*x)/k;
	term = fact/k;
	sum += term;
	if (term < EPS*sum) 
	  break;
      }
      if (k > MAXIT) 
	nrerror("Series failed in ei");
      
      *eix = sum+log(*x)+EULER;      
    } 
    else {
      sum = 0.0;
      term = 1.0;
      for (k=1; k<=MAXIT; k++) {
	prev = term;
	term *= k/(*x);
	if (term < EPS) 
	  break;
	if (term < prev) 
	  sum += term;
	else {
	  sum -= prev;
	  break;
	}
      }    
      *eix = exp(*x)*(1.0+sum)/(*x);
    }
  } /* endif x>0 */
}

#undef EPS
#undef EULER
#undef MAXIT
#undef FPMIN
