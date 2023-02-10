/* Adaptive multidimensional integration of a vector of integrands.
 *
 * Copyright (c) 2005-2009 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <R.h>
#include <Rmath.h>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <time.h>

#define LOW -1.0e15
#define MINF -1.0e15
#define EPS DBL_EPSILON

/*----------------------------------------------------------------
File name: simschlather.c
Description: Simulation from Schlather max-stable model
Start
 ---------------------------------------------------------------*/

  void rschlathertbm(double *coord, int *nObs, int *nSites, int *dim,
                     int *covmod, int *grid, double *nugget, double *range,
                     double *smooth, double *uBound, int *nlines,
                     double *ans);
  void rschlatherdirect(double *coord, int *nObs, int *nSites, int *dim,
                        int *covmod, int *grid, double *nugget, double *range,
                        double *smooth, double *uBound, double *ans, int *ans2);
  void rschlathercirc(int *nObs, int *ngrid, double *steps, int *dim,
                      int *covmod, double *nugget, double *range,
                      double *smooth, double *uBound, double *ans);
  void tbmcore(int *nsite, int *neffSite, int *dim, int *covmod,
               int *grid, double *coord, double *nugget, double *sill,
               double *range, double *smooth, int *nlines, double *lines,
               double *ans);
  void circcore(double *rho, double *a, double *ia, int m, int halfM, int mdag,
                int mdagbar, int ngrid, int nbar, double isqrtMbar, double nugget,
                double *ans);

/*----------------------------------------------------------------
File name: simgeometric.c
Description: Simulation from Geometric Gaussian max-stable model
Start
 ---------------------------------------------------------------*/

  void rgeomtbm(double *coord, int *nObs, int *nSite, int *dim,
                int *covmod, int *grid, double *sigma2, double *nugget,
                double *range, double *smooth, double *uBound, int *nlines,
                double *ans);
  void rgeomdirect(double *coord, int *nObs, int *nSite, int *dim,
                   int *covmod, int *grid, double *sigma2, double *nugget,
                   double *range, double *smooth, double *uBound,
                   double *ans, int *ans2);
  void rgeomcirc(int *nObs, int *ngrid, double *steps, int *dim,
                 int *covmod, double *sigma2, double *nugget, double *range,
                 double *smooth, double *uBound, double *ans);

/*----------------------------------------------------------------
File name: simBrownResnick.c
Description: Simulation from Brown-Resnick max-stable model
Start
 ---------------------------------------------------------------*/

  void rbrowndirect(double *coord, double *bounds, int *nObs, int *nSite,
                    int *dim, int *grid, double *range, double *smooth,
                    double *uBound, int *method, int *maxSim, int *nPP,
                    int *idxsubOrig, int *nsubOrig, double *ans, int *ans2);

/*----------------------------------------------------------------
File name: simextremalt.c
Description: Simulation from Extremal (skew)-t max-stable model
Start
 ---------------------------------------------------------------*/

  void rextremalttbm(double *coord, int *nObs, int *nSite, int *dim,
                     int *covmod, int *grid, double *nugget, double *range,
                     double *smooth, double *DoF, double *uBound, int *nlines,
                     double *ans);
  void rextremaltdirect(double *coord, int *nObs, int *nSite, int *dim,
                        int *covmod, int *grid, double *nugget, double *range,
                        double *smooth, double *DoF, double *uBound, double *ans, int *ans2);
  void rextremalskewtdirect(double *coord, int *nObs, int *nSite, int *dim,
                            int *covmod, int *grid, double *nugget, double *range,
                            double *smooth, double *DoF, double *alpha, double *uBound, double *ans, int *ans2);
  void rextremaltcirc(int *nObs, int *ngrid, double *steps, int *dim,
                      int *covmod, double *nugget, double *range,
                      double *smooth, double *DoF, double *uBound, double *ans);

/*----------------------------------------------------------------
File name: maxStableExactSim.c
Description: Exact Simulation from max-stable models
Start
 ---------------------------------------------------------------*/

  void rbrownexact(double *coord, int *nObs, int *nSite, int *dim,
                   int *grid, double *range, double *smooth,
                   double *ans, int *ans2);
  void rextremaltexact(double *coord, int *nObs, int *nSite, int *dim,
                       int *covmod, int *grid, double *nugget, double *range,
                       double *smooth, double *DoF, int *cholsky, double *ans, int *ans2);
  void rextremalskewtexact(double *coord, int *nObs, int *nSite, int *dim,
                           int *covmod, int *grid, double *nugget, double *range,
                           double *smooth, double *DoF, double *alpha, int *cholsky, double *ans, int *ans2);
  void rschlatherexact(double *coord, int *nObs, int *nSite, int *dim,
                       int *covmod, int *grid, double *nugget, double *range,
                       double *smooth, double *ans, int *ans2);
  void rgeomexact(double *coord, int *nObs, int *nSite, int *dim,
                  int *covmod, int *grid, double *sigma2, double *nugget,
                  double *range, double *smooth, double *ans, int *ans2);

/*----------------------------------------------------------------
File name: spatial.c
Description: Exact Simulation from max-stable models
Start
---------------------------------------------------------------*/

  void buildcovmat(int *nSite, int *grid, int *covmod, double *coord, int *dim,
                   double *nugget, double *sill, double *range, double *smooth,
                   double *covMat);
  void distance(double *coord, int *nDim, int *nSite, int *vec, double *dist);
  void distance2orig(double *coord, int n, int dim, double *dist, int grid);
  void getSiteIndex(int currentPair, int nSite, int *site1, int *site2);
  /* LinAlg - Linear algebra functions: column-major order */
  void R_smult(double *A, double *B, int *m, int *n, int *p, double *C);
  void R_qform(double *A, double *B, int *m, int *p, double *C);
  void R_xmult(double *A, double *B, int *m, int *n, int *p, double *C);
  /* covariance */
  double whittleMatern(double *dist, int n, double nugget, double sill, double range,
  		                 double smooth, double *rho);
  double cauchy(double *dist, int n, double nugget, double sill, double range,
  	            double smooth, double *rho);
  double powerExp(double *dist, int n, double nugget, double sill, double range,
          		    double smooth, double *rho);
  double bessel(double *dist, int n, int dim, double nugget, double sill,
          	    double range, double smooth, double *rho);
  double fbm(double *coord, double *dist, int dim, int nSite, double sill, double range,
        	   double smooth, double *rho);
  /* randomlines */
  void vandercorput(int *n, double *coord);
  void rotation(double *coord, int *n, double *u, double *v, double *w,
  	            double *angle);

/*----------------------------------------------------------------
File name: fft.c
Description: Exact Simulation from max-stable models
Start
---------------------------------------------------------------*/

  void fft_factor(int n, int *pmaxf, int *pmaxp);
  Rboolean fft_work(double *a, double *b, int nseg, int n, int nspn, int isn,
  		              double *work, int *iwork);
