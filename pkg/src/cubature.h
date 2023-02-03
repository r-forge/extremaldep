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

#define LOW -1.0e15

#ifndef CUBATURE_H
#define CUBATURE_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* USAGE: Call adapt_integrate with your function as described below.

	  To compile a test program, compile cubature.c with
	  -DTEST_INTEGRATOR as described at the end. */

/* a vector integrand - evaluates the function at the given point x
  (an array of length ndim) and returns the result in fval (an array
  of length fdim).   The void* parameter is there in case you have
  to pass any additional data through to your function (it corresponds
  to the fdata parameter you pass to adapt_integrate). */
  typedef void (*integrand) (unsigned ndim, const double *x, void *,
			   unsigned fdim, double *fval);


/* Integrate the function f from xmin[dim] to xmax[dim], with at most
  maxEval function evaluations (0 for no limit), until the given
  absolute or relative error is achieved.  val returns the integral,
  and err returns the estimate for the absolute error in val; both
  of these are arrays of length fdim, the dimension of the vector
  integrand f(x). The return value of the function is 0 on success
  and non-zero if there  was an error. */
  int adapt_integrate(unsigned fdim, integrand f, void *fdata,
                        unsigned dim, const double *xmin, const double *xmax,
                        unsigned maxEval, double reqAbsError, double reqRelError,
                        double *val, double *err);

  // univariate Skew-normal functions:
  double desn_int(double x, double mu, double omega, double alpha, double tau);
  void desn_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
  void pesn(double *xmin, double *xmax, double *par, double *val, double *err);
  void desn(double *x, double *mu, double *omega, double *alpha, double *tau, double *res);
  // multivariate Skew-normal functions:
  double dmn_int(double x[2], double rho, double std);
  double dmn_int3(double x[3], double rho[3], double std[3]);
  double dmesn_int(double x[2], double mu[2], double omega[4], double alpha[2], double tau);
  double dmesn_int3(double x[3], double mu[3], double omega[9], double alpha[3], double tau);
  void dmesn(double *x, double *mu, double *omega, double *alpha, double *tau, double *res);
  void dmesn3(double *x, double *mu, double *omega, double *alpha, double *tau, double *res);
  void dmesn_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
  void dmesn_t3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
  void pmesn(double *xmin, double *xmax, double *par, double *val, double *err);
  void pmesn3(double *xmin, double *xmax, double *par, double *val, double *err);
  double pmesn_int(double *par);
  double pmesn_int3(double *par);
  // univariate Skew-t functions:
  double dest_int(double x, double mu, double omega, double nu, double alpha, double tau);
  void dest_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
  void pest(double *xmin, double *xmax, double *par, double *val, double *err);
  double pest_int(double *par);
  void dest(double *x, double *mu, double *omega, double *nu,double *alpha, double *tau, double *res);
  // multivariate Skew-t functions:
  double dmt_int(double nu, double omr2, double Q, double std);
  double dmt_int3(double nu, double Q, double det);
  double dmest_int(double x[2], double mu[2], double omega[4], double nu, double alpha[2], double tau);
  double dmest_int3(double x[3], double mu[3], double omega[9], double nu, double alpha[3], double tau);
  void dmest(double *x, double *mu, double *omega, double *nu, double *alpha, double *tau, double *res);
  void dmest3(double *x, double *mu, double *omega, double *nu, double *alpha, double *tau, double *res);
  void dmest_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
  void dmest_t3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
  void pmest(double *xmin, double *xmax, double *par, double *val, double *err);
  void pmest3(double *xmin, double *xmax, double *par, double *val, double *err);
  double pmest_int(double *par);
  double pmest_int3(double *par);
  // bivariate Pickand dependence function for the extremal skew-t model:
  void chistup(double *omega, double *nu, double *alpha, double *res);
  void chistlo(double *omega, double *nu, double *alpha, double *res);
  void bivpkst(double *x, double *omega, double *nu, double *alpha, double *res);
  void trivpkst(double *x, double *omega, double *nu, double *alpha, double *res);
  // extremal skew-t functions
  void dmextst(double *x, double *omega, double *nu, double *alpha, double *res);
  void pmextst(double *x, double *omega, double *nu, double *alpha, double *res);
  void llextst(double *x, int *n, double *omega, double *nu, double *alpha, double *res);
  double dmextst_int(double *x, double *omega, double *nu, double *alpha);
  double d1x_dt(double x, double df);

  double HuslerReiss(double a, double *x);
  void pHuslerReiss(double *x, double *lambda, double *res);
  void dHuslerReiss(double *x, double *lambda, double *res);
  void llHRmax(double *x, double *lambda, int *n, double *res);
  double ExtremalT(double *data, double df, double rho);
  void llETmax(double *x, double *par, int *n, double *res);

/*----------------------------------------------------------------
File name: CompRandFld.c
Description: procedures for the computation of useful distributions
Start
 ---------------------------------------------------------------*/

  extern double dgev(double x, double loc, double scale, double shape);
  extern double pgev(double x, double loc, double scale, double shape);
  extern double qgev(double x, double loc, double scale, double shape);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* CUBATURE_H */
