#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
/* #include "spatial.h" */

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
extern void pHuslerReiss(void *, void *, void *);
extern void dHuslerReiss(void *, void *, void *);
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

/*----------------------------------------------------------------
File name: CompRandFld.c
 ---------------------------------------------------------------*/

extern void GevLogLik(double *data, int *ndata, double *par, double *res);
extern void dGEV(double *x, double *loc, double *scale, double *shape, double *res);
extern void pGEV(double *x, double *loc, double *scale, double *shape, double *res);
extern void qGEV(double *x, double *loc, double *scale, double *shape, double *res);

/*----------------------------------------------------------------
File name: simschlather.c
Description: Simulation from Schlather max-stable model
Start
 ---------------------------------------------------------------*/

extern void rschlathertbm(double *coord, int *nObs, int *nSites, int *dim,
                     int *covmod, int *grid, double *nugget, double *range,
                     double *smooth, double *uBound, int *nlines,
                     double *ans);
extern void rschlatherdirect(double *coord, int *nObs, int *nSites, int *dim,
                        int *covmod, int *grid, double *nugget, double *range,
                        double *smooth, double *uBound, double *ans, int *ans2);
extern void rschlathercirc(int *nObs, int *ngrid, double *steps, int *dim,
                      int *covmod, double *nugget, double *range,
                      double *smooth, double *uBound, double *ans);

/*----------------------------------------------------------------
File name: simgeometric.c
Description: Simulation from Geometric Gaussian max-stable model
Start
 ---------------------------------------------------------------*/

extern void rgeomtbm(double *coord, int *nObs, int *nSite, int *dim,
                int *covmod, int *grid, double *sigma2, double *nugget,
                double *range, double *smooth, double *uBound, int *nlines,
                double *ans);
extern void rgeomdirect(double *coord, int *nObs, int *nSite, int *dim,
                   int *covmod, int *grid, double *sigma2, double *nugget,
                   double *range, double *smooth, double *uBound,
                   double *ans, int *ans2);
extern void rgeomcirc(int *nObs, int *ngrid, double *steps, int *dim,
                 int *covmod, double *sigma2, double *nugget, double *range,
                 double *smooth, double *uBound, double *ans);

/*----------------------------------------------------------------
File name: simBrownResnick.c
Description: Simulation from Brown-Resnick max-stable model
Start
 ---------------------------------------------------------------*/

extern void rbrowndirect(double *coord, double *bounds, int *nObs, int *nSite,
                    int *dim, int *grid, double *range, double *smooth,
                    double *uBound, int *method, int *maxSim, int *nPP,
                    int *idxsubOrig, int *nsubOrig, double *ans, int *ans2);

/*----------------------------------------------------------------
File name: simextremalt.c
Description: Simulation from Extremal (skew)-t max-stable model
Start
 ---------------------------------------------------------------*/

extern void rextremalttbm(double *coord, int *nObs, int *nSite, int *dim,
                     int *covmod, int *grid, double *nugget, double *range,
                     double *smooth, double *DoF, double *uBound, int *nlines,
                     double *ans);
extern void rextremaltdirect(double *coord, int *nObs, int *nSite, int *dim,
                        int *covmod, int *grid, double *nugget, double *range,
                        double *smooth, double *DoF, double *uBound, double *ans, int *ans2);
extern void rextremalskewtdirect(double *coord, int *nObs, int *nSite, int *dim,
                            int *covmod, int *grid, double *nugget, double *range,
                            double *smooth, double *DoF, double *alpha, double *uBound, double *ans, int *ans2);
extern void rextremaltcirc(int *nObs, int *ngrid, double *steps, int *dim,
                      int *covmod, double *nugget, double *range,
                      double *smooth, double *DoF, double *uBound, double *ans);

/*----------------------------------------------------------------
File name: maxStableExactSim.c
Description: Exact Simulation from max-stable models
Start
 ---------------------------------------------------------------*/

extern void rbrownexact(double *coord, int *nObs, int *nSite, int *dim,
                   int *grid, double *range, double *smooth,
                   double *ans, int *ans2);
extern void rextremaltexact(double *coord, int *nObs, int *nSite, int *dim,
                       int *covmod, int *grid, double *nugget, double *range,
                       double *smooth, double *DoF, int *cholsky, double *ans, int *ans2);
extern void rextremalskewtexact(double *coord, int *nObs, int *nSite, int *dim,
                           int *covmod, int *grid, double *nugget, double *range,
                           double *smooth, double *DoF, double *alpha, int *cholsky, double *ans, int *ans2);
extern void rschlatherexact(void *, void *, void *, void *,
                       void *, void *, void *, void *,
                       void *, void *, void *);
extern void rgeomexact(double *coord, int *nObs, int *nSite, int *dim,
                  int *covmod, int *grid, double *sigma2, double *nugget,
                  double *range, double *smooth, double *ans, int *ans2);

/*----------------------------------------------------------------
File name: maxStableExactSim.c
Description: Exact Simulation from max-stable models
Start
---------------------------------------------------------------*/

extern void mypmvnorm(double *upper, int *d, double *chol,
                      int *Nmax, int *Nmin, double *eps, int *logeps, double *out);

extern void mypmvt(double *upper, int *d, double *chol, double *df,
                   int *Nmax, int *Nmin, double *eps, int *logeps, double *out);


static const R_CMethodDef CEntries[] = {
    {"chistup",         (DL_FUNC) &chistup,          4},
    {"chistlo",         (DL_FUNC) &chistlo,          4},
    {"dmextst",         (DL_FUNC) &dmextst,          5},
    {"desn",            (DL_FUNC) &desn,             6},
    {"dest",            (DL_FUNC) &dest,             7},
    {"dmesn",           (DL_FUNC) &dmesn,            6},
    {"dmesn3",          (DL_FUNC) &dmesn3,           6},
    {"dmest",           (DL_FUNC) &dmest,            7},
    {"dmest3",          (DL_FUNC) &dmest3,           7},
    {"pmextst",         (DL_FUNC) &pmextst,          5},
    {"pHuslerReiss",    (DL_FUNC) &pHuslerReiss,     3},
    {"dHuslerReiss",    (DL_FUNC) &dHuslerReiss,     3},
    {"llHRmax",         (DL_FUNC) &llHRmax,          4},
    {"llETmax",         (DL_FUNC) &llETmax,          4},
    {"llextst",         (DL_FUNC) &llextst,          6},
    {"pesn",            (DL_FUNC) &pesn,             5},
    {"pest",            (DL_FUNC) &pest,             5},
    {"bivpkst",         (DL_FUNC) &bivpkst,          5},
    {"trivpkst",        (DL_FUNC) &trivpkst,         5},
    {"pmesn",           (DL_FUNC) &pmesn ,           5},
    {"pmesn3",          (DL_FUNC) &pmesn3,           5},
    {"pmest",           (DL_FUNC) &pmest,            5},
    {"pmest3",          (DL_FUNC) &pmest3,           5},
    {"doCubature",      (DL_FUNC) &doCubature,       8},
    /* ------------- CompRandFld.c ------*/
    {"GevLogLik",       (DL_FUNC) &GevLogLik,        4},
    {"dGEV",            (DL_FUNC) &dGEV,             5},
    {"pGEV",            (DL_FUNC) &pGEV,             5},
    {"qGEV",            (DL_FUNC) &qGEV,             5},
    /* ------------- simschlather.c ------*/
    {"rschlathertbm",            (DL_FUNC) &rschlathertbm,            12},
    {"rschlatherdirect",         (DL_FUNC) &rschlatherdirect,         12},
    {"rschlatherexact",          (DL_FUNC) &rschlatherexact,          11},
    /* ------------- simgeometric.c ------*/
    {"rgeomtbm",                 (DL_FUNC) &rgeomtbm,                 13},
    {"rgeomdirect",              (DL_FUNC) &rgeomdirect,              13},
    {"rgeomcirc",                (DL_FUNC) &rgeomcirc,                11},
    /* ------------- simBrownResnick.c ------*/
    {"rbrowndirect",             (DL_FUNC) &rbrowndirect,             16},
    /* ------------- simextremalt.c ------*/
    {"rextremalttbm",            (DL_FUNC) &rextremalttbm,            13},
    {"rextremaltdirect",         (DL_FUNC) &rextremaltdirect,         13},
    {"rextremalskewtdirect",     (DL_FUNC) &rextremalskewtdirect,     14},
    {"rextremaltcirc",           (DL_FUNC) &rextremaltcirc,           11},
    /* ------------- maxStableExactSim.c ------*/
    {"rbrownexact",              (DL_FUNC) &rbrownexact,               9},
    {"rextremaltexact",          (DL_FUNC) &rextremaltexact,          13},
    {"rextremalskewtexact",      (DL_FUNC) &rextremalskewtexact,      14},
    {"rschlatherexact",          (DL_FUNC) &rschlatherexact,          11},
    {"rgeomexact",               (DL_FUNC) &rgeomexact,               12},
    /* ------------- mypmvnorm.c ------*/
    {"mypmvnorm",                (DL_FUNC) &mypmvnorm,                 8},
    {"mypmt",                    (DL_FUNC) &mypmvt,                    9},
    {NULL, NULL, 0}
};

void R_init_cubing(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
