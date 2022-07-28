#include <Rcpp.h>
using namespace Rcpp;

//' Update transitions
//'
//' @param waifw who acquires immunity from whom matrix
//' @param state state
//' @param s_indsR something
//' @param i_indsR something
//' @param exponentR something
//' @param denomR something
//' @param tran_matrix transition matrix
//' @param N number of rows in transition matrix
//' @return transition matrix
// [[Rcpp::export]]
SEXP calc_phi_and_update_tran (SEXP waifw,
                               SEXP state,
                               SEXP s_indsR,
                               SEXP i_indsR,
                               SEXP exponentR,
                               SEXP denomR,
                               SEXP tran_matrix) {

  int i,j;
  int ind;
  int n_age_class = Rf_length(i_indsR);
  int *i_inds = INTEGER(i_indsR);
  int *s_inds = INTEGER(s_indsR);
  double exponent = REAL(exponentR)[0];
  double denom = REAL(denomR)[0];
  int N = INTEGER(GET_DIM(tran_matrix))[0];

  //SEXP rc;

  waifw = Rf_coerceVector(waifw, REALSXP);

  //do not need phi outside of function
  //double *phi = malloc(n_age_class * sizeof(double));
  NumericVector phi(n_age_class);

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
  PROTECT(tran_matrix = Rf_coerceVector(tran_matrix, REALSXP));
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

