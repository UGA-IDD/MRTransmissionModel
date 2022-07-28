#include <Rcpp.h>
using namespace Rcpp;

//' Update age and survival matrix for MSIRV sims
//'
//' @param age_surv_matrix the matrix to update
//' @param sz something
//' @param age_spec_vacc_prob the age specific probability of vaccination
//' @param v_indsR index for vaccination
//' @param m_indsR index for maternal antibodies
//' @param s_indsR index for susceptibles
//' @param i_indsR index for infectious
//' @return something
// [[Rcpp::export]]
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
  int len = Rf_length(v_indsR);
  int *v_inds = INTEGER(v_indsR);
  int *m_inds = INTEGER(m_indsR);
  int *s_inds = INTEGER(s_indsR);
  int *i_inds = INTEGER(i_indsR);
  SEXP rc;

  //Rprintf("%d\n",I);

  age_spec_vacc_prob = Rf_coerceVector(age_spec_vacc_prob, REALSXP);
  age_surv_matrix = Rf_coerceVector(age_surv_matrix, REALSXP);

  PROTECT(rc=Rf_allocMatrix(REALSXP,I,I));
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
