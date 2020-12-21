#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP RFGLS_BFcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RFGLS_invZcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RFGLSpredicttree_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RFGLStree_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"RFGLS_BFcpp",          (DL_FUNC) &RFGLS_BFcpp,          10},
  {"RFGLS_invZcpp",        (DL_FUNC) &RFGLS_invZcpp,         7},
  {"RFGLSpredicttree_cpp", (DL_FUNC) &RFGLSpredicttree_cpp,  9},
  {"RFGLStree_cpp",        (DL_FUNC) &RFGLStree_cpp,        20},
  {NULL, NULL, 0}
};

void R_init_RandomForestsGLS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
