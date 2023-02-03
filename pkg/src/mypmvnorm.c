#include <R.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

/* produces multivariate normal cdf */
void mypmvnorm(double *upper, int *d, double *chol, 
               int *Nmax, int *Nmin, double *eps, int *logeps, double *out)
  
{
  double delta,sval,neps;
  int i,j;
  int N=0;
  double intsum=0, varsum=0;
  
  double *y,*e,*f;
  y = (double *)R_alloc(*d, sizeof(double));
  e = (double *)R_alloc(*d, sizeof(double));
  f = (double *)R_alloc(*d, sizeof(double));
  
  e[0] = f[0] = pnorm(upper[0] / chol[0], 0, 1, 1, 0);
  
  RANDIN;
  
  do {
    for(i=1;i<*d;i++) 
    { 
      y[i-1] = qnorm(UNIF * e[i-1], 0, 1, 1, 0);
	    sval = 0;
	    for(j=0;j<i;j++) {
	      sval = sval + chol[i * *d + j] * y[j];
	    }
      e[i] = pnorm((upper[i] - sval) / chol[i * *d + i], 0, 1, 1, 0);
	    f[i] = e[i] * f[i-1];
    }
    N++;
	  delta = (f[*d - 1] - intsum)/N;
    intsum = intsum + delta;
    //Rprintf("%f ", delta);
    varsum = (N - 2) * varsum/N + delta*delta;
    if(*logeps) 
      neps = intsum * *eps; 
    else neps = *eps;
  } while ((varsum > (neps * neps/6.25) || N < *Nmin) && N != *Nmax);
  
  *out = intsum;
  
  RANDOUT;			  
}
			  
/* produces multivariate t cdf */
/* note that qt() and pt() is central t only and therefore has no ncp argument */
/* for non-central t the entry points are qnt() and pnt() */
void mypmvt(double *upper, int *d, double *chol, double *df,
            int *Nmax, int *Nmin, double *eps, int *logeps, double *out)
  
{
  double delta,sval,neps,sysqval,conval;
	int i,j;
	int N=0;
	double intsum=0, varsum=0;
			    
	double *y,*e,*f;
	y = (double *)R_alloc(*d, sizeof(double));
	e = (double *)R_alloc(*d, sizeof(double));
	f = (double *)R_alloc(*d, sizeof(double));
            
	e[0] = f[0] = pt(upper[0] / chol[0], *df, 1, 0);
			    
	RANDIN;
			    
	do {
	  sysqval = 0;
		for(i=1;i<*d;i++) 
		{ 
		  conval = sqrt((*df + sysqval) / (*df + i - 1));
			y[i-1] = qt(UNIF * e[i-1], *df + i - 1, 1, 0) * conval;
		  sysqval = sysqval + y[i-1]*y[i-1];
			sval = 0;
			for(j=0;j<i;j++) {
			   sval = sval + chol[i * *d + j] * y[j];
			}
			conval = sqrt((*df + i) / (*df + sysqval));
			e[i] = pt((upper[i] - sval) * conval / chol[i * *d + i], *df + i, 1, 0);
			f[i] = e[i] * f[i-1];
		}
		N++;
		delta = (f[*d - 1] - intsum)/N;
		intsum = intsum + delta;
		//Rprintf("%f ", delta);
		varsum = (N - 2) * varsum/N + delta*delta;
		if(*logeps) 
			neps = intsum * *eps; 
		else neps = *eps;
	} while ((varsum > (neps * neps/6.25) || N < *Nmin) && N != *Nmax);
			    
	*out = intsum;
			    
	RANDOUT;			  
}			  










