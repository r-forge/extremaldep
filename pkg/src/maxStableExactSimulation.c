#include "spatial.h"

void rbrownexact(double *coord, int *nObs, int *nSite, int *dim,
		 int *grid, double *range, double *smooth,
		 double *ans, int *ans2){

  /*
    This function generates random fields for the Brown-Resnick model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double zero = 0, one = 1, irange = 1 / *range;
  int covmod = 6;//i.e, fractional Brownian motion

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
    *vario = malloc(neffSite * sizeof(double)),
    *shiftedCoord = malloc(*nSite * *dim * sizeof(double)),
    *orig = malloc(*dim * sizeof(double)),
    *poisson = malloc(*nObs * sizeof(double));

  int *iter = malloc(*nObs * sizeof(double));
  for (int i=0; i<*nObs; i++) iter[i] = 0;

  buildcovmat(nSite, grid, &covmod, coord, dim, &zero, &one, range,
	      smooth, covmat);

  /* Compute the Cholesky decomposition of the covariance matrix once for all */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info FCONE);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();
  for (int j=0;j<neffSite;j++){
    // Set the origin
    if (*grid){
      int idx1 = j / *nSite, idx2 = j % *nSite;//works only for 2d grid
      orig[0] = coord[idx1];
      orig[1] = coord[*nSite + idx2];
    }
    else {
      for (int d=0;d<*dim;d++)
	orig[d] = coord[j + d * *nSite];
    }

    // Compute the variogram gamma(s - origin)
    for (int l=0; l<*nSite;l++)
      for (int d=0; d<*dim; d++)
	shiftedCoord[d * *nSite + l] = coord[d * *nSite + l] - orig[d];

    distance2orig(shiftedCoord, *nSite, *dim, vario, *grid);

    for (int l=0; l<neffSite; l++)
      vario[l] = R_pow(vario[l] * irange, *smooth);

    for (int i=0; i<*nObs; i++){

      poisson[i] = exp_rand();
      double ipoisson = -log(poisson[i]);

      while (ans[j * lagj + i * lagi] < ipoisson){
		  iter[i]++;
	R_CheckUserInterrupt();

	// Generate a proposal extremal function
	for (int l=0;l<neffSite;l++)
	  gp[l] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt

      FCONE FCONE FCONE);

	double dummy = gp[j];
	for (int l=0;l<neffSite;l++)
	  gp[l] -= dummy + vario[l];

	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if ((ipoisson + gp[l]) > ans[l * lagj + i * lagi]){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++) {
		if((ipoisson + gp[l]) > ans[l * lagj + i * lagi]) ans2[l * lagj + i * lagi] = iter[i];
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], ipoisson + gp[l]);
	  }

	// Update the "norm" of the spectral function (log-scale)
	poisson[i] += exp_rand();
	ipoisson = -log(poisson[i]);
      }
    }
  }

  // Switch to unit Frechet margins
  for (int i=0;i<(*nObs * neffSite);i++)
    ans[i] = exp(ans[i]);

  PutRNGstate();
  free(covmat); free(gp); free(vario); free(shiftedCoord); free(orig); free(poisson);
  free(iter);

  return;
}


void rextremaltexact(double *coord, int *nObs, int *nSite, int *dim,
		     int *covmod, int *grid, double *nugget, double *range,
		     double *smooth, double *DoF, int *cholsky, double *ans, int *ans2){

  /*
    This function generates random fields for the extremal-t model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    covmod: the covariance model
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    DoF: the degree of freedom nu (NOT nu_d = nu + 1)
	cholsky: Cholsky (TRUE) or eigen (FALSE) decomposition
    ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double sill = 1 - *nugget,
    miDoF = - 1.0 / *DoF;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *scalemat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
	*gp2 = malloc(neffSite * sizeof(double)),
    *mu = malloc(neffSite * sizeof(double)),
    *poisson = malloc(*nObs * sizeof(double));

  int lwork = 10*neffSite + 2*neffSite*neffSite;
  double *evals = malloc(neffSite * sizeof(double)),
	*wkvec = malloc(lwork * sizeof(double));

  int *iter = malloc(*nObs * sizeof(double));
  for (int i=0; i<*nObs; i++) iter[i] = 0;

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);

  GetRNGstate();
  for (int j=0;j<neffSite;j++){
    // Get the mean vector of the extremal function
    for (int l=0;l<neffSite;l++)
      mu[l] = covmat[l + j * neffSite];

    // Compute the scale matrix of the extremal function
    for (int l=0;l<neffSite;l++)
      for (int m=l;m<neffSite;m++)
	scalemat[l + m * neffSite] = scalemat[m + l * neffSite] =
	  (covmat[l + m * neffSite] - covmat[j + l * neffSite] * covmat[j + m * neffSite]) /
	  (1.0 + *DoF);

	/*
	Rprintf("Original Matrix\n");
	for (int l=0;l<neffSite;l++) {
      for (int m=0;m<neffSite;m++) {
		  Rprintf("%f ", scalemat[l + m * neffSite]);
	  }
	  Rprintf("\n");
	}
	*/


	if(*cholsky == 1) {

      // By construction this matrix is singular since the j-th row is 0
      // so we regularize it
      scalemat[j + j * neffSite] = 1e-12;

      // Compute the Cholesky decomposition of this matrix
      int info = 0;
      F77_CALL(dpotrf)("U", &neffSite, scalemat, &neffSite, &info FCONE);

      if (info != 0)
        error("error code %d from Lapack routine '%s'", info, "dpotrf");

      // Reset it to the degenerate case
      scalemat[j + j * neffSite] = 0;
	} else {

	  // Compute the eigenvalue decomposition of this matrix
	  // Uses only the upper triangle of scalemat, but output uses
	  // whole of scalemat as the eigenvector matrix
	  int info = 0;
      F77_CALL(dsyev)("V", "U", &neffSite, scalemat, &neffSite, evals, wkvec, &lwork, &info
        FCONE FCONE);

      if (info != 0)
        error("error code %d from Lapack routine '%s'", info, "dsyev");

      for (int m=0;m<neffSite;m++)
	    evals[m] = fmax2(evals[m], 0.0);

      for (int l=0;l<neffSite;l++) {
        for (int m=0;m<neffSite;m++) {
		  scalemat[l + m * neffSite] = scalemat[l + m * neffSite] * sqrt(evals[m]);
	    }
      }
	}

	/*
	Rprintf("Decomposition Matrix\n");
	for (int l=0;l<neffSite;l++) {
      for (int m=0;m<neffSite;m++) {
		  Rprintf("%f ", scalemat[l + m * neffSite]);
	  }
	  Rprintf("\n");
	}
	*/


    for (int i=0; i<*nObs; i++){
      poisson[i] = exp_rand();
      double ipoissonDoF = R_pow(poisson[i], miDoF);

      while (ans[j * lagj + i * lagi] < ipoissonDoF){
		  iter[i]++;
	R_CheckUserInterrupt();

	// Generate a proposal extremal function (log-scale)
	for (int l=0;l<neffSite;l++)
	  gp[l] = norm_rand();

    /*
    Rprintf("Original Vector\n");
	for (int l=0;l<neffSite;l++) {
		  Rprintf("%f ", gp[l]);
	}
	Rprintf("\n");
	*/

    if(*cholsky == 1) {
	    F77_CALL(dtrmv)("U", "T", "N", &neffSite, scalemat, &neffSite, gp, &oneInt
        FCONE FCONE FCONE);
	} else {
		R_smult(scalemat, gp, &neffSite, &oneInt, &neffSite, gp2);
		for (int l=0;l<neffSite;l++) gp[l] = gp2[l];
	}

	/*
	Rprintf("Final Vector\n");
	for (int l=0;l<neffSite;l++) {
		  Rprintf("%f ", gp[l]);
	}
	Rprintf("\n");
	*/

	double scaleDoF = sqrt((1 + *DoF) / rchisq(1 + *DoF));
	for (int l=0;l<neffSite;l++)
	  gp[l] = mu[l] + gp[l] * scaleDoF;

	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if ((ipoissonDoF * gp[l]) > ans[l * lagj + i * lagi]){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++) {
		if((gp[l] * ipoissonDoF) > ans[l * lagj + i * lagi]) ans2[l * lagj + i * lagi] = iter[i];
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], gp[l] * ipoissonDoF);
	  }

	// Update the "norm" of the spectral function (1/DoF scale)
	poisson[i] += exp_rand();
	ipoissonDoF = R_pow(poisson[i], miDoF);
      }
    }
  }

  // Switch to unit Frechet margins
  for (int i=0;i<(*nObs * neffSite);i++)
    ans[i] = R_pow(ans[i], *DoF);

  PutRNGstate();
  free(covmat); free(scalemat); free(gp); free(mu); free(poisson);
  free(iter);

  free(evals); free(gp2); free(wkvec);

  return;
}

void rextremalskewtexact(double *coord, int *nObs, int *nSite, int *dim,
		     int *covmod, int *grid, double *nugget, double *range,
		     double *smooth, double *DoF, double *alpha, int *cholsky,
			 double *ans, int *ans2){

  /*
    This function generates random fields for the extremal-skew-t model
    using the exact procedure in Algorithm 1 of Dombry et al. (2106)
	"Exact simulation of max-stable processes" Biometrika


    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    covmod: the covariance model
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    the degree of freedom nu (NOT nu_d = nu + 1)
	cholsky: Cholsky (TRUE) or eigen (FALSE) decomposition
    ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double sill = 1 - *nugget,
    miDoF = - 1.0 / *DoF;

  double quadform, extval, ttsum, tt;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *scalemat = malloc(neffSite * neffSite * sizeof(double)),
	*gp = malloc(neffSite * sizeof(double)),
	*gp2 = malloc(neffSite * sizeof(double)),
	*mu = malloc(neffSite * sizeof(double)),
	*poisson = malloc(*nObs * sizeof(double)),
	*sdvec = malloc(neffSite * sizeof(double)),
	*mplus = malloc(neffSite * sizeof(double));

  int lwork = 10*neffSite + 2*neffSite*neffSite;
  double *evals = malloc(neffSite * sizeof(double)),
	*wkvec = malloc(lwork * sizeof(double));

  int *iter = malloc(*nObs * sizeof(double));
  for (int i=0; i<*nObs; i++) iter[i] = 0;

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);

  // calculate mplus vector for skew-t

  for(int j=0;j<neffSite;j++)
  {
	// numerator of alphastar
	int sval = 0;
	for(int m=0;m<neffSite;m++) {
		sval = sval + covmat[m + j * neffSite] * alpha[m];
	}

	// scaled matrix with zeros in jth row and column
	for (int l=0;l<neffSite;l++) {
      for (int m=l;m<neffSite;m++) {
	    scalemat[l + m * neffSite] = scalemat[m + l * neffSite] =
	      (covmat[l + m * neffSite] - covmat[j + l * neffSite] * covmat[j + m * neffSite]);
	  }
	}

	// quadratic form in denominator of alphastar
	R_qform(alpha, scalemat, &oneInt, &neffSite, &quadform);

	// alphastar
    mplus[j] = sval / sqrt(1 + quadform);

	// mplus
	// M_LN_SQRT_PI is log(pi)/2
	mplus[j] = exp((*DoF/2) * log(2) - M_LN_SQRT_PI + lgammafn((*DoF + 1)/2) +
      pt(mplus[j] * sqrt(*DoF + 1), *DoF + 1, 1, 1));

  }

  GetRNGstate();
  for (int j=0;j<neffSite;j++){

	//double mpj1 = R_pow(mplus[j], 1 / *DoF);
	//double mpj2 = R_pow(mplus[j], 2 / *DoF);

	double mpj1 = 1;
	double mpj2 = 1;

	// Get the mean vector of the extremal function
    for (int l=0;l<neffSite;l++)
      mu[l] = covmat[l + j * neffSite] * mpj1;

    // Get the extension parameter of the extremal function

	extval = 0;
    for (int l=0;l<neffSite;l++)
      extval += covmat[l + j * neffSite] * alpha[l];
    extval = sqrt(*DoF + 1) * extval;

    // Compute the scale matrix of the extremal function
    for (int l=0;l<neffSite;l++)
      for (int m=l;m<neffSite;m++)
				scalemat[l + m * neffSite] = scalemat[m + l * neffSite] =
				(covmat[l + m * neffSite] - covmat[j + l * neffSite] * covmat[j + m * neffSite]) / (1.0 + *DoF) * mpj2;

	// Normalize to its correlation and save the standard deviations
	// For shifted skew-t we need the scale matrix to be normalized
	for (int l=0;l<neffSite;l++)
	  sdvec[l] = sqrt(scalemat[l + l * neffSite]);

    // only need the upper triangle here but doing it all anyway
    for (int l=0;l<neffSite;l++)
      for (int m=0;m<neffSite;m++)
		 if(l != j && m != j) {
		scalemat[l + m * neffSite] =
	      scalemat[l + m * neffSite] / (sdvec[l] * sdvec[m]);
		 }

	/* Option One: Cholesky Decomposition */
	/* Option Two: Eigenvalue Decomposition */

    if(*cholsky == 1) {

      // By construction this matrix is singular since the j-th row is 0
      // so we regularize it
      scalemat[j + j * neffSite] = 1e-12;

      // Compute the Cholesky decomposition of this matrix
	  // Uses only the upper triangle
      int info = 0;
      F77_CALL(dpotrf)("U", &neffSite, scalemat, &neffSite, &info FCONE);

      if (info != 0)
        error("error code %d from Lapack routine '%s'", info, "dpotrf");

      // Reset it to the degenerate case
      scalemat[j + j * neffSite] = 0;
	} else {

	  // Compute the eigenvalue decomposition of this matrix
	  // Uses only the upper triangle of scalemat, but output uses
	  // whole of scalemat as the eigenvector matrix
	  int info = 0;
      F77_CALL(dsyev)("V", "U", &neffSite, scalemat, &neffSite, evals, wkvec, &lwork, &info
         FCONE FCONE);

      if (info != 0)
        error("error code %d from Lapack routine '%s'", info, "dsyev");

      for (int l=0;l<neffSite;l++) {
        for (int m=0;m<neffSite;m++) {
		  scalemat[l + m * neffSite] = scalemat[l + m * neffSite] * sqrt(evals[m]);
	    }
      }
	}

    for (int i=0; i<*nObs; i++){
      poisson[i] = exp_rand();
      double ipoissonDoF = R_pow(poisson[i] / mplus[j], miDoF);


      while (ans[j * lagj + i * lagi] < ipoissonDoF){
		  iter[i]++;
		  R_CheckUserInterrupt();

	do {

	  // Generate a proposal extremal function (log-scale)
	  for (int l=0;l<neffSite;l++)
	    gp[l] = norm_rand();

	  if(*cholsky == 1) {
	    F77_CALL(dtrmv)("U", "T", "N", &neffSite, scalemat, &neffSite, gp, &oneInt
        FCONE FCONE FCONE);
	  } else {
	    R_smult(scalemat, gp, &neffSite, &oneInt, &neffSite, gp2);
		for (int l=0;l<neffSite;l++) gp[l] = gp2[l];
	  }

	  double scaleDoF = sqrt((1 + *DoF) / rchisq(1 + *DoF));
	  for (int l=0;l<neffSite;l++)
		  gp[l] = gp[l] * scaleDoF;

	  // For extended skew-normal
	  ttsum = 0;
	  tt = rt(1 + *DoF);
	  for (int l=0;l<neffSite;l++)
	    ttsum += alpha[l] * gp[l];

	} while(tt >= sqrt(*DoF + 1.0) * ttsum + extval);
	// } while(tt >= ttsum + extval);

	for (int l=0;l<neffSite;l++) gp[l] = mu[l] + sdvec[l] * gp[l];


	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if ((ipoissonDoF * gp[l]) > ans[l * lagj + i * lagi]){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++) {
		if((gp[l] * ipoissonDoF) > ans[l * lagj + i * lagi]) ans2[l * lagj + i * lagi] = iter[i];
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], gp[l] * ipoissonDoF);
	  }

	// Update the "norm" of the spectral function (1/DoF scale)
	poisson[i] += exp_rand();
	ipoissonDoF = R_pow(poisson[i] / mplus[j], miDoF);

      }
    }
  }

  // Switch to unit Frechet margins
  for (int i=0;i<(*nObs * neffSite);i++) {
    ans[i] = R_pow(ans[i], *DoF);
  }
  for (int i=0; i<*nObs; i++){
	for (int j=0;j<neffSite;j++) {
	  ans[j * lagj + i * lagi] = ans[j * lagj + i * lagi] / mplus[j];
	}
  }

  PutRNGstate();
  free(covmat); free(scalemat); free(gp); free(mu); free(poisson);
  free(iter);

  free(mplus); free(sdvec);
  free(evals); free(gp2); free(wkvec);

  return;
}

void rschlatherexact(double *coord, int *nObs, int *nSite, int *dim,
		     int *covmod, int *grid, double *nugget, double *range,
		     double *smooth, double *ans, int *ans2){

  /*
    This function generates random fields for the Schlather model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    covmod: the covariance model
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double sill = 1 - *nugget;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *scalemat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
    *mu = malloc(neffSite * sizeof(double)),
    *poisson = malloc(*nObs * sizeof(double));

  int *iter = malloc(*nObs * sizeof(double));
  for (int i=0; i<*nObs; i++) iter[i] = 0;

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);

  GetRNGstate();
  for (int j=0;j<neffSite;j++){
    // Get the mean vector of the extremal function
    for (int l=0;l<neffSite;l++)
      mu[l] = covmat[l + j * neffSite];

    // Compute the scale matrix of the extremal function
    for (int l=0;l<neffSite;l++)
      for (int m=l;m<neffSite;m++)
	scalemat[l + m * neffSite] = scalemat[m + l * neffSite] =
	  (covmat[l + m * neffSite] - covmat[j + l * neffSite] * covmat[j + m * neffSite]) / 2.0;

    // By construction this matrix is singular since the j-th row is 0
    // so we regularize it
    scalemat[j + j * neffSite] = 1e-12;

    // Compute the Cholesky decomposition of this matrix
    int info = 0;
    F77_CALL(dpotrf)("U", &neffSite, scalemat, &neffSite, &info FCONE);

    if (info != 0)
      error("error code %d from Lapack routine '%s'", info, "dpotrf");

    // Reset it to the degenerate case
    scalemat[j + j * neffSite] = 0;

    for (int i=0; i<*nObs; i++){
      poisson[i] = exp_rand();

      while ((poisson[i] * ans[j * lagj + i * lagi]) < 1){
	    iter[i]++;
	R_CheckUserInterrupt();

	// Generate a proposal extremal function (log-scale)
	for (int l=0;l<neffSite;l++)
	  gp[l] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, scalemat, &neffSite, gp, &oneInt
     FCONE FCONE FCONE);

    double scale = sqrt(2.0 / rchisq(2.0));
	for (int l=0;l<neffSite;l++)
	  gp[l] = mu[l] + gp[l] * scale;

	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if (gp[l] > (poisson[i] * ans[l * lagj + i * lagi])){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++) {
		if((gp[l] / poisson[i]) > ans[l * lagj + i * lagi]) ans2[l * lagj + i * lagi] = iter[i];
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], gp[l] / poisson[i]);
	  }

	// Update the "norm" of the spectral function
	poisson[i] += exp_rand();
      }
    }
  }

  PutRNGstate();
  free(covmat); free(scalemat); free(gp); free(mu); free(poisson);
  free(iter);

  return;
}

/* Alec function for geometric; adapted from rbrownexact */
void rgeomexact(double *coord, int *nObs, int *nSite, int *dim,
		 int *covmod, int *grid, double *sigma2, double *nugget,
		 double *range, double *smooth, double *ans, int *ans2){

  /*
    This function generates random fields for the geom gaussian model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  const double one = 1, dzero = 0;
  double flag = 0;
  double sill = 1 - *nugget;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
    *vario = malloc(neffSite * sizeof(double)),
	*rho = malloc(neffSite * sizeof(double)),
    *shiftedCoord = malloc(*nSite * *dim * sizeof(double)),
    *orig = malloc(*dim * sizeof(double)),
    *poisson = malloc(*nObs * sizeof(double));

  int *iter = malloc(*nObs * sizeof(double));
  for (int i=0; i<*nObs; i++) iter[i] = 0;

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);

  /* Compute the Cholesky decomposition of the covariance matrix once for all */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info FCONE);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();
  for (int j=0;j<neffSite;j++){
    // Set the origin
    if (*grid){
      int idx1 = j / *nSite, idx2 = j % *nSite;//works only for 2d grid
      orig[0] = coord[idx1];
      orig[1] = coord[*nSite + idx2];
    }
    else {
      for (int d=0;d<*dim;d++)
	orig[d] = coord[j + d * *nSite];
    }

    // Compute the variogram gamma(s - origin)
    for (int l=0; l<*nSite;l++)
      for (int d=0; d<*dim; d++)
	shiftedCoord[d * *nSite + l] = coord[d * *nSite + l] - orig[d];

    distance2orig(shiftedCoord, *nSite, *dim, vario, *grid);


	switch (*covmod){
      case 1:
        flag = whittleMatern(vario, neffSite, dzero, one, *range, *smooth, rho);
        break;
      case 2:
        flag = cauchy(vario, neffSite, dzero, one, *range, *smooth, rho);
        break;
      case 3:
        flag = powerExp(vario, neffSite, dzero, one, *range, *smooth, rho);
        break;
      case 4:
        flag = bessel(vario, neffSite, *dim, dzero, one, *range, *smooth, rho);
        break;
	}
	if (flag != 0.0)
      error("The covariance parameters seem to be ill-defined. Please check\n");
	for (int l=0; l<neffSite; l++) {
	  vario[l] = *sigma2 * (1 - rho[l]);
	}

    for (int i=0; i<*nObs; i++){

      poisson[i] = exp_rand();
      double ipoisson = -log(poisson[i]);

      while (ans[j * lagj + i * lagi] < ipoisson){
		  iter[i]++;
	R_CheckUserInterrupt();

	// Generate a proposal extremal function
	for (int l=0;l<neffSite;l++)
	  gp[l] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt 
    FCONE FCONE FCONE);

	double dummy = gp[j];
	for (int l=0;l<neffSite;l++)
	  gp[l] -= dummy + vario[l];

	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if ((ipoisson + gp[l]) > ans[l * lagj + i * lagi]){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++) {
		if((ipoisson + gp[l]) > ans[l * lagj + i * lagi]) ans2[l * lagj + i * lagi] = iter[i];
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], ipoisson + gp[l]);
	  }

	// Update the "norm" of the spectral function (log-scale)
	poisson[i] += exp_rand();
	ipoisson = -log(poisson[i]);
      }
    }
  }

  // Switch to unit Frechet margins
  for (int i=0;i<(*nObs * neffSite);i++)
    ans[i] = exp(ans[i]);

  PutRNGstate();
  free(covmat); free(gp); free(vario); free(shiftedCoord); free(orig); free(poisson);
  free(iter); free(rho);

  return;
}
