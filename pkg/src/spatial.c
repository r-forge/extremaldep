#include "spatial.h"

void buildcovmat(int *nSite, int *grid, int *covmod, double *coord, int *dim,
		 double *nugget, double *sill, double *range,
		 double *smooth, double *covMat){

  int nPairs, effnSite = *nSite, zero = 0;
  const double one = 1, dzero = 0;
  double flag = 0;

  if (*grid)
    effnSite = R_pow_di(effnSite, *dim);

  nPairs = effnSite * (effnSite - 1) / 2;

  double *dist = malloc(nPairs * sizeof(double)),
    *rho = malloc(nPairs * sizeof(double)),
    *coordGrid = malloc(effnSite * *dim * sizeof(double));

  if (*grid){
    //Coord specify a grid
    for (int i = 0; i < *nSite; i++)
      for (int j = 0; j < *nSite; j++){
	coordGrid[j + i * *nSite] = coord[i];
	coordGrid[*nSite * (*nSite + i) + j] = coord[j];
      }

    distance(coordGrid, dim, &effnSite, &zero, dist);
  }

  else
    //Coord don't specify a grid
    distance(coord, dim, nSite, &zero, dist);

  switch (*covmod){
  case 1:
    flag = whittleMatern(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 2:
    flag = cauchy(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 3:
    flag = powerExp(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 4:
    flag = bessel(dist, nPairs, *dim, dzero, one, *range, *smooth, rho);
    break;
  case 6:
    if (*grid)
      flag = fbm(coordGrid, dist, *dim, effnSite, one, *range, *smooth, rho);

    else
      flag = fbm(coord, dist, *dim, effnSite, one, *range, *smooth, rho);

    break;
  }

  if (flag != 0.0)
    error("The covariance parameters seem to be ill-defined. Please check\n");

  //Fill the non-diagonal elements of the covariance matrix
  //#pragma omp parallel for
  for (int currentPair=0;currentPair<nPairs;currentPair++){
    int i = 0, j = 0;
    getSiteIndex(currentPair, effnSite, &i, &j);
    covMat[effnSite * i + j] = covMat[effnSite * j + i] = *sill * rho[currentPair];
  }

  //Fill the diagonal elements of the covariance matrix
  if (*covmod == 6){
    //Fractional brownian
    double irange2 = 1 / (*range * *range);

    if (*grid){
      for (int i = 0; i < effnSite;i++){
	covMat[i * (effnSite + 1)] = 0;

	for (int j= 0; j < *dim; j++)
	  covMat[i * (effnSite + 1)] += coordGrid[i + j * effnSite] * coordGrid[i + j * effnSite];

	covMat[i * (effnSite + 1)] = 2 * pow(covMat[i * (effnSite + 1)] * irange2, 0.5 * *smooth);
      }
    }

    else {
      for (int i = 0; i < effnSite; i++){
	covMat[i * (effnSite + 1)] = 0;

	for (int j = 0; j < *dim; j++)
	  covMat[i * (effnSite + 1)] += coord[i + j * effnSite] * coord[i + j * effnSite];

	covMat[i * (effnSite + 1)] = 2 * pow(covMat[i * (effnSite + 1)] * irange2, 0.5 * *smooth);
      }
    }
  }

  else
    for (int i = 0; i < effnSite; i++)
      covMat[i * (effnSite + 1)] = *sill + *nugget;


  free(dist); free(rho); free(coordGrid);
  return;
}

void distance(double *coord, int *nDim, int *nSite,
	      int *vec, double *dist){

  //This function computes either the euclidean distance or the
  //distance vector between each pair of locations
  const int nPair = *nSite * (*nSite - 1) / 2;

  if (*vec){
    // #pragma omp parallel for
    for (int pair=0;pair<nPair;pair++){
      int i, j;
      getSiteIndex(pair, *nSite, &i, &j);

      for (int k=0;k<*nDim;k++)
	dist[k * nPair + pair] = coord[k * *nSite + j] - coord[k * *nSite + i];
    }
  }

  else{
    //#pragma omp parallel for
    for (int pair=0;pair<nPair;pair++){
      int i, j;
      getSiteIndex(pair, *nSite, &i, &j);
      dist[pair] = 0;

      for (int k=0;k<*nDim;k++)
	dist[pair] += (coord[i + k * *nSite] -	coord[j + k * *nSite]) *
	  (coord[i + k * *nSite] - coord[j + k * *nSite]);

      dist[pair] = sqrt(dist[pair]);
    }
  }
}

void distance2orig(double *coord, int n, int dim, double *dist, int grid){
  //It computes the l_2 norm of points in R^d i.e. sqrt(x[1]^2 + ... + x[d]^2)
  if (grid){
    //Only works with two dimensional grids!!!
    int current = -1;
    double dummy;

    for (int i=0;i<n;i++){
      dummy = coord[i] * coord[i];

      for (int j=0;j<n;j++){
	current++;
	dist[current] = sqrt(dummy + coord[n + j] * coord[n + j]);
      }
    }
  }

  else {
    for (int i=0;i<n;i++){
      dist[i] = 0;
      for (int j=0;j<dim;j++)
	dist[i] += coord[i + j * n] * coord[i + j * n];

      dist[i] = sqrt(dist[i]);
    }
  }
}

void getSiteIndex(int currentPair, int nSite, int *site1, int *site2){
  int nFree = nSite - 2,
      cum = nSite - 2;

  *site1 = 0;

  while (currentPair > cum){
    (*site1)++;
    cum += nFree;
    nFree--;
  }

  *site2 = *site1 + currentPair - cum + nFree + 1;

  return;
}


/* LINEAR ALGEBRA ROUTINES: COLUMN-MAJOR ORDER */

/* Calculates A %*% B where A is m by p and B is p by n */
/* Stores result in m by n vector C */
/* Uses dgemm from BLAS */
void R_smult(double *A, double *B, int *m, int *n, int *p, double *C)
{

  double one = 1.0;
  double zero = 0.0;

  /*
  Rprintf("m = %d, n = %d, p = %d\n",*m,*n,*p);
	int i;
	for(i=0;i<(*m**p);i++){
		Rprintf("%f ",A[i]);
	}
	Rprintf("\n");
	for(i=0;i<(*p**n);i++){
		Rprintf("%f ",B[i]);
	}
	Rprintf("\n");
	for(i=0;i<(*m**n);i++){
		Rprintf("%f ",C[i]);
	}
	Rprintf("\n");
    */

	F77_CALL(dgemm)("N","N",m,n,p,&one,A,m,B,p,&zero,C,m FCONE FCONE);
}


/* Calculates A %*% B %*% t(A) where A is m by p and B is p by p */
/* Stores result in m by m vector C */
void R_qform(double *A, double *B, int *m, int *p, double *C)
{
  double *C2;
  C2 = (double *)Calloc(*m * *p * sizeof(double), double);

  R_smult(A, B, m, p, p, C2);
  R_xmult(C2, A, m, m, p, C);
  Free(C2);
}

/* Calculates A %*% t(B) where A is m by p and B is n by p */
/* Stores result in m by n vector C */
/* Uses dgemm from BLAS */
void R_xmult(double *A, double *B, int *m, int *n, int *p, double *C)
{

  double one = 1.0;
  double zero = 0.0;

  /*
  Rprintf("m = %d, n = %d, p = %d\n",*m,*n,*p);
  int i;
	for(i=0;i<(*p**m);i++){
		Rprintf("%f ",A[i]);
	}
	Rprintf("\n");
	for(i=0;i<(*p**n);i++){
		Rprintf("%f ",B[i]);
	}
	Rprintf("\n");
	for(i=0;i<(*m**n);i++){
		Rprintf("%f ",C[i]);
	}
	Rprintf("\n");
  */

	F77_CALL(dgemm)("N","T",m,n,p,&one,A,m,B,n,&zero,C,m FCONE FCONE);
}


double whittleMatern(double *dist, int n, double nugget, double sill, double range,
		     double smooth, double *rho){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When ans != 0.0, the whittle-matern parameters are ill-defined.

  const double cst = sill * R_pow(2, 1 - smooth) / gammafn(smooth),
    irange = 1 / range;

  //Some preliminary steps: Valid points?
  if (smooth < EPS)
    return (1 - smooth + EPS) * (1 - smooth + EPS) * MINF;

  else if (smooth > 100)
    /* Not really required but larger smooth parameters are unlikely
       to occur */
    return (smooth - 99) * (smooth - 99) * MINF;

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  //#pragma omp parallel for
  for (int i=0;i<n;i++){
    double cst2 = dist[i] * irange;

    if (cst2 == 0)
      rho[i] = sill + nugget;

    else
      rho[i] = cst * R_pow(cst2, smooth) * bessel_k(cst2, smooth, 1);
  }

  return 0.0;
}

double cauchy(double *dist, int n, double nugget, double sill, double range,
	      double smooth, double *rho){

  //This function computes the cauchy covariance function between each
  //pair of locations.
  //When ans != 0.0, the cauchy parameters are ill-defined.

  const double irange2 = 1 / (range * range);

  //Some preliminary steps: Valid points?
  if (smooth < 0)
    return (1 - smooth) * (1 - smooth) * MINF;

  else if (smooth > 100)
    return (smooth - 99) * (smooth - 99) * MINF;

  if (range <= 0.0)
    return (1 - range) * (1 - range)* MINF;

  if (sill <= 0.0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  //#pragma omp parallel for
  for (int i=0;i<n;i++){

    if (dist[i] == 0)
      rho[i] = nugget + sill;

    else
      rho[i] = sill * R_pow(1 + dist[i] * dist[i] * irange2, -smooth);
  }

  return 0.0;
}

double powerExp(double *dist, int n, double nugget, double sill, double range,
		double smooth, double *rho){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  const double irange = 1 / range;

  //Some preliminary steps: Valid points?
  if ((smooth < 0) || (smooth > 2))
    return (1 - smooth) * (1 - smooth) * MINF;

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  //#pragma omp parallel for
  for (int i=0;i<n;i++){
    if (dist[i] == 0)
      rho[i] = nugget + sill;

    else
      rho[i] = sill * exp(-R_pow(dist[i] * irange, smooth));
  }

  return 0.0;
}

double bessel(double *dist, int n, int dim, double nugget, double sill,
	      double range, double smooth, double *rho){
  //This function computes the bessel covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  const double irange = 1 / range, cst = sill * R_pow(2, smooth) * gammafn(smooth + 1);

  //Some preliminary steps: Valid points?
  if (smooth < (0.5 * (dim - 2)))
    return (1 + 0.5 * (dim - 2) - smooth) * (1 + 0.5 * (dim - 2) - smooth) * MINF;

  /* else if (smooth > 100)
    //Require as bessel_j will be numerically undefined
    return (smooth - 99) * (smooth - 99) * MINF; */

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  if (nugget < 0)
    return (1 - nugget) * (1 - nugget) * MINF;

  //#pragma omp parallel for
  for (int i=0;i<n;i++){
    double cst2 = dist[i] * irange;

    if (cst2 == 0)
      rho[i] = nugget + sill;

    else if (cst2 <= 1e5)
      rho[i] = cst * R_pow(cst2, -smooth) * bessel_j(cst2, smooth);

    else
      // approximation of the besselJ function for large x
      rho[i] = cst * R_pow(cst2, -smooth) * M_SQRT_2dPI *
	cos(cst2 - smooth * M_PI_2 - M_PI_4);

    /*if (!R_FINITE(rho[i]))
      return MINF;*/
  }

  return 0.0;
}

double fbm(double *coord, double *dist, int dim, int nSite, double sill, double range,
	   double smooth, double *rho){

  /* This function computes the covariance function related to a
  fractional Brownian motion.  When ans != 0.0, the parameters are
  ill-defined. */

  const int nPairs = nSite * (nSite - 1) / 2;
  const double irange = 1 / range;
  double *distOrig = malloc(nSite * sizeof(double));

  //Some preliminary steps: Valid points?
  if (smooth < EPS)
    return (1 - smooth + EPS) * (1 - smooth + EPS) * MINF;

  else if (smooth > 2)
    return (smooth - 1) * (smooth - 1) * MINF;

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  distance2orig(coord, nSite, dim, distOrig, 0);

  /* Rmk: 0.5 Var[Z(x)] = \gamma(||x||) as Z(o) = 0, where \gamma is
     the semi-variogram */
  //#pragma omp parallel for
  for (int i=0;i<nSite;i++)
    distOrig[i] = pow(distOrig[i] * irange, smooth);

  //#pragma omp parallel for
  for (int currentPair=0;currentPair<nPairs;currentPair++){
    int i, j;
    getSiteIndex(currentPair, nSite, &i, &j);

    rho[currentPair] = sill * (distOrig[i] + distOrig[j] -
			       pow(dist[currentPair] * irange, smooth));
  }

  free(distOrig);

  return 0;
}


void vandercorput(int *n, double *coord){
  //This function computes the Van der Corput sequence in R^3
  int i, k, r;
  double base, u, v;

  for (i=*n;i--;){
    //Binary decomposition
    k = i;
    u = 0;
    base = 2.;

    while (k){
      r = k % 2;
      u += r / base;
      base *= 2;
      k = k / 2;
    }

    //Ternary decomposition
    k = i;
    v = 0;
    base = 3;

    while (k){
      r = k % 3;
      v += r / base;
      base *= 3;
      k = k / 3;
    }

    coord[i] = cos(M_2PI * u) * sqrt(1 - v * v);
    coord[*n + i] = sin(M_2PI * u) * sqrt(1 - v * v);
    coord[2 * *n + i] = v;
  }

  return;
}

void rotation(double *coord, int *n, double *u, double *v,
	      double *w, double *angle){
  /* This function performs a rotation of coord w.r.t. to a rotation
     axis (u, v, w) and a rotation angle */

  int i;
  double r, px, py, pz, bx, by, bz, cx, cy, cz,
    cosAngle = cos(*angle), sinAngle = sin(*angle);


  for (i=*n;i--;){
    //Projection of the points onto the rotation axis
    r = coord[i] * *u + coord[*n + i] * *v + coord[2 * *n + i] * *w;

    px = r * *u;
    py = r * *v;
    pz = r * *w;

    r = sqrt((coord[i] - px) * (coord[i] - px) +
	     (coord[*n + i] - py) * (coord[*n + i] - py) +
	     (coord[2 * *n + i] - pz) * (coord[2 * *n + i] - pz));

    bx = (coord[i] - px) / r;
    by = (coord[*n + i] - py) / r;
    bz = (coord[2 * *n + i] - pz) / r;

    cx = *v * bz - *w * by;
    cy = *w * bx - *u * bz;
    cz = *u * by - *v * bx;

    coord[i] = px + r * cosAngle * bx + r * sinAngle * cx;
    coord[*n + i] = py + r * cosAngle * by + r * sinAngle * cy;
    coord[2 * *n + i] = pz + r * cosAngle * bz + r * sinAngle * cz;
  }

  return;
}
