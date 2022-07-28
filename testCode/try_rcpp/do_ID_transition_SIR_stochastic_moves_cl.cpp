#include <Rcpp.h>
using namespace Rcpp;

//' Make stochastic transition moves
//'
//' @param curstateR current state
//' @param transmatrixR transition matrix
//' @return something
// [[Rcpp::export]]
SEXP do_ID_transition_SIR_stochastic_moves_cl (SEXP curstateR,
                                               SEXP transmatrixR) {
  int i,j,n;
  int len = length(curstateR);
  int *curstate = INTEGER(curstateR);
  int *drw = malloc(len * sizeof(int));
  SEXP newstate;


  double accum = 0.0;
  double pp;

  PROTECT(newstate = allocVector(INTSXP, len));

  /*Initialize new state to 0....might not be necessary ultimately*/
  for (int i=0; i<len; i++) {
    INTEGER(newstate)[i] = 0;
  }

  /* /\*Loop over all states*\/ */
  GetRNGstate();

  for (i=0; i<len; i++) {
    /*We do the multinomial as a series of binomial draws...this
     allows us to not draw mortality...which is just the remainder*/
    n = curstate[i];
    if (n == 0) continue; //avioid useless draws
    accum=0.0;


    for (j=0;j<len;j++) {
      pp=((accum>=1.)?(0.0):(REAL(transmatrixR)[i* (len) + j]/(1-accum)));
      drw[j] = ((pp<1.0)?((pp>0)?(int)rbinom((double)n,pp):0):n);
      n -=drw[j];
      accum = accum + REAL(transmatrixR)[i* (len) + j];
    }


    for (j=0;j<len;j++) {
      INTEGER(newstate)[j] += (int)drw[j];
    }
  }


  UNPROTECT(1);

  free(drw);

  return newstate;


}
