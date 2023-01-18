#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP calc_phi_and_update_tran(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP do_ID_transition_SIR_stochastic_moves_cl(SEXP, SEXP);
extern SEXP update_age_surv_MSIRV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"calc_phi_and_update_tran",                 (DL_FUNC) &calc_phi_and_update_tran,                 7},
  {"do_ID_transition_SIR_stochastic_moves_cl", (DL_FUNC) &do_ID_transition_SIR_stochastic_moves_cl, 2},
  {"update_age_surv_MSIRV",                    (DL_FUNC) &update_age_surv_MSIRV,                    7},
  {NULL, NULL, 0}
};

void R_init_MRTransmissionModel(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
