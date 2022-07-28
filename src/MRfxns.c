#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/*******************************************************/
/* function for getting the next state of the model from*/
/* this one stochastically                              */
/*                                                      */
/* IMPORTANT: ASSUMES STATE FORWARDNESS, THAT IS...NO   */
/* ONE CAN TRANSITION TO A STATE THAT COMES BEFORE THEIR*/
/* POSITION IN THE MATRIX                              */
/*                                                      */
/*@param curstate the curent state of the system        */
/*@param len  the length of the base states             */
/*@param transmatrix the transistion and motality       */
/*     probabilities.                                   */
/*@param tm_ht the ht of the transition matrix           */
/*@param newstate the newstate of the system after the call */
/*@param n_states the number of states forward to look  */
/*    for a transition. Based on the idea that only     */
/*    transitions of one age class forward are possible. */
/*    If negative 1, then this is ignored               */
/********************************************************/


/*.Call version*/
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

/*Updates the ageing and survival matrix
 for the MSIRV simulations.

 @param age_surv_matrix the matrix to update
 @param age_spec_vacc_prob the age specific probability of vaccination
 @param v_inds index for vaccination
 @param m_inds index for maternal antibodies
 @param s_inds index for suscetibles
 @param i_inds index for infectious*/
SEXP update_age_surv_MSIRV(SEXP age_surv_matrix,
                           SEXP sz,
                           SEXP age_spec_vacc_prob,
                           SEXP v_indsR,
                           SEXP m_indsR,
                           SEXP s_indsR,
                           SEXP i_indsR) {


  int i,j;
  int ind;
  //int I = INTEGER(GET_DIM(age_surv_matrix))[0];
  int I = INTEGER(sz)[0];
  int len = length(v_indsR);
  int *v_inds = INTEGER(v_indsR);
  int *m_inds = INTEGER(m_indsR);
  int *s_inds = INTEGER(s_indsR);
  int *i_inds = INTEGER(i_indsR);
  SEXP rc;


  //Rprintf("%d\n",I);

  age_spec_vacc_prob = coerceVector(age_spec_vacc_prob, REALSXP);
  age_surv_matrix = coerceVector(age_surv_matrix, REALSXP);

  PROTECT(rc=allocMatrix(REALSXP,I,I));
  for (i=0; i<I*I; i++) {
    REAL(rc)[i] = REAL(age_surv_matrix)[i];
  }

  for (j=0;j<len;j++) {
    for (i=0; i<len;i++) {
      //M->V transition
      ind = (v_inds[i]-1)+I * (m_inds[j]-1);
      //REAL(age_surv_matrix)[ind] =
      REAL(rc)[ind] =
        REAL(age_surv_matrix)[ind] *  REAL(age_spec_vacc_prob)[i];



      //M->M transition
      ind = (m_inds[i]-1)+I*(m_inds[j]-1);
      REAL(rc)[ind] =
        REAL(age_surv_matrix)[ind] *  (1-REAL(age_spec_vacc_prob)[i]);

      //M->S transition
      ind = (s_inds[i]-1)+I*(m_inds[j]-1);
      REAL(rc)[ind] =
        REAL(age_surv_matrix)[ind] *  (1-REAL(age_spec_vacc_prob)[i]);

      //S->V
      ind = (v_inds[i]-1)+I*(s_inds[j]-1);
      REAL(rc)[ind] =
        REAL(age_surv_matrix)[ind] *  (REAL(age_spec_vacc_prob)[i]);

      //S->S
      ind = (s_inds[i]-1)+I*(s_inds[j]-1);
      REAL(rc)[ind] =
        REAL(age_surv_matrix)[ind] *  (1-REAL(age_spec_vacc_prob)[i]);

      //S->I
      ind = (i_inds[i]-1)+I*(s_inds[j]-1);
      REAL(rc)[ind] =
        REAL(age_surv_matrix)[ind] *  (1-REAL(age_spec_vacc_prob)[i]);

    }
  }
  //Rprintf("\n");

  /* for (i=0;i<10;i++) { */
  /*   for (j=0; j<10;j++) { */
  /*     ind = I*j + i; */
  /*     //Rprintf("(%d,%d,%d)%1.6f",i,j,ind,REAL(age_surv_matrix)[ind]); */
  /*     Rprintf("%1.6f",REAL(rc)[ind]); */
  /*   } */
  /*   Rprintf("\n"); */
  /* } */

  UNPROTECT(1);
  return(rc);
}


/*Function calculates the phi matrix and then
 updates all of the transitions*/
SEXP calc_phi_and_update_tran (SEXP waifw,
                               SEXP state,
                               SEXP s_indsR,
                               SEXP i_indsR,
                               SEXP exponentR,
                               SEXP denomR,
                               SEXP tran_matrix) {

  int i,j;
  int ind;
  int n_age_class = length(i_indsR);
  int *i_inds = INTEGER(i_indsR);
  int *s_inds = INTEGER(s_indsR);
  double exponent = REAL(exponentR)[0];
  double denom = REAL(denomR)[0];
  int N = INTEGER(GET_DIM(tran_matrix))[0];

  //SEXP rc;

  waifw = coerceVector(waifw, REALSXP);

  //do not need phi outside of function
  double *phi = malloc(n_age_class * sizeof(double));

  for (i=0; i<n_age_class; i++) {
    phi[i] = 0.0;

    for (j=0; j<n_age_class; j++) {
      phi[i] += REAL(waifw)[i + j* n_age_class] *
        pow(REAL(state)[i_inds[j]-1], exponent)/denom;
    }
  }

  for (i=0; i<n_age_class; i++) {
    phi[i] =  1-exp(-phi[i]);
  }

  //Now that phi is calculated, update the
  //tran matrix
  PROTECT(tran_matrix = coerceVector(tran_matrix, REALSXP));
  //PROTECT(rc = allocMatrix(REALSXP,N,N));

  //for (i=0; i<N*N; i++) {
  //   REAL(rc)[i] = REAL(tran_matrix)[i];
  //}

  for (i=0; i<n_age_class; i++) {
    for (j=0; j<n_age_class; j++) {

      //people who stay susceptible
      ind = (s_inds[i] - 1) + N * (s_inds[j]-1);
      //REAL(rc)[ind] =
      REAL(tran_matrix)[ind] *=(1 - phi[j]);

      //people who become infected
      ind = (i_inds[i] - 1) + N * (s_inds[j]-1);
      //REAL(rc)[ind] =
      REAL(tran_matrix)[ind] *= phi[j];

      //no one stays infected
      ind = (i_inds[i] - 1) + N * (i_inds[j]-1);
      //REAL(rc)[ind] = 0;
      REAL(tran_matrix)[ind] = 0;
    }
  }

  free(phi);

  UNPROTECT(1);
  return tran_matrix;
}
