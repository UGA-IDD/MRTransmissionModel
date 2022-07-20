#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/*Avoid the need to deal with mortality*/
void do_ID_transition_SIR_stochastic_moves(int *curstate,
                                           int *len,
                                           double *transmatrix,
                                           int *newstate) {
  int i,j,n;
  /*int *drw_probs = malloc(*tm_width * sizeof(double));*/
  int *drw = malloc(*len * sizeof(int));
  double accum = 0.0;
  double pp;

  /*Initialize new state to 0....might not be necessary ultimately*/
  for (i=0; i< *len; i++) {
    newstate[i] = 0;
  }


  /*Loop over all states*/
  GetRNGstate();
  for (i = 0; i < *len; i++) {
    /*We do the multinomial as a series of binomial draws...this
     allows us to not draw mortality...which is just the remainder*/
    n = curstate[i];

    accum=0.0;
    for (j=0;j<*len;j++) {
      pp=((accum>=1.)?(0.0):(transmatrix[i* (*len) + j]/(1-accum)));
      drw[j] = ((pp<1.0)?(int)rbinom((double)n,pp):n);
      n -=drw[j];
      accum = accum + transmatrix[i* (*len) + j];
    }


    /*update the new state...assumes state forwardness*/
    for (j=0;j<*len;j++) {
      newstate[j] = newstate[j]+(int)drw[j];
    }
  }
  PutRNGstate();

  free(drw);
}
