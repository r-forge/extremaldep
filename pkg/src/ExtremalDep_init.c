#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void chistup(void *, void *, void *, void *);
extern void chistlo(void *, void *, void *, void *);
extern void dmextst(void *, void *, void *, void *, void *);
extern void desn(void *, void *, void *, void *, void *, void *);
extern void dest(void *, void *, void *, void *, void *, void *, void *);
extern void dmesn(void *, void *, void *, void *, void *, void *);
extern void dmesn3(void *, void *, void *, void *, void *, void *);
extern void dmest(void *, void *, void *, void *, void *, void *, void *);
extern void dmest3(void *, void *, void *, void *, void *, void *, void *);
extern void pmextst(void *, void *, void *, void *, void *);
extern void llHRmax(void *, void *, void *, void *);
extern void llETmax(void *, void *, void *, void *);
extern void llextst(void *, void *, void *, void *, void *, void *);
extern void pesn(void *, void *, void *, void *, void *);
extern void pest(void *, void *, void *, void *, void *);
extern void bivpkst(void *, void *, void *, void *, void *);
extern void trivpkst(void *, void *, void *, void *, void *);
extern void pmesn(void *, void *, void *, void *, void *);
extern void pmesn3(void *, void *, void *, void *, void *);
extern void pmest(void *, void *, void *, void *, void *);
extern void pmest3(void *, void *, void *, void *, void *);
extern void doCubature(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"chistup",    (DL_FUNC) &chistup,    4},
    {"chistlo",    (DL_FUNC) &chistlo,    4},
    {"dmextst",    (DL_FUNC) &dmextst,    5},
    {"desn",       (DL_FUNC) &desn,       6},
    {"dest",       (DL_FUNC) &dest,       7},
    {"dmesn",      (DL_FUNC) &dmesn,      6},
    {"dmesn3",     (DL_FUNC) &dmesn3,     6},
    {"dmest",      (DL_FUNC) &dmest,      7},
    {"dmest3",     (DL_FUNC) &dmest3,     7},
    {"pmextst",    (DL_FUNC) &pmextst,    5},
    {"llHRmax",    (DL_FUNC) &llHRmax,    4},
    {"llETmax",    (DL_FUNC) &llETmax,    4},
    {"llextst",    (DL_FUNC) &llextst,    6},
    {"pesn",       (DL_FUNC) &pesn,       5},
    {"pest",       (DL_FUNC) &pest,       5},
    {"bivpkst",    (DL_FUNC) &bivpkst,    5},
    {"trivpkst",   (DL_FUNC) &trivpkst,   5},
    {"pmesn",      (DL_FUNC) &pmesn ,     5},
    {"pmesn3",     (DL_FUNC) &pmesn3,     5},
    {"pmest",      (DL_FUNC) &pmest,      5},
    {"pmest3",     (DL_FUNC) &pmest3,     5},
    {"doCubature", (DL_FUNC) &doCubature, 8},
    {NULL, NULL, 0}
};

void R_init_cubing(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
