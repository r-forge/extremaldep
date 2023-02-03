#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

#include "cubature.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("mylib", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif




/* START ------ SKEW-NORMAL FUNCTIONS -------------- */

// univariate extended skew-normal probability density function:
double desn_int(double x, double mu, double omega, double alpha, double tau)
{
  double res=0, y=(x-mu)/omega;
  //computes the univariate extended skew-normal density:
  res=dnorm(x,mu,omega,0)*pnorm(alpha*y+tau,0,1,1,0)/pnorm(tau/sqrt(1+alpha*alpha),0,1,1,0);
  return res;
}
// interface for univariate extended sk-norm pdf:
void desn(double *x, double *mu, double *omega, double *alpha, double *tau, double *res)
{
  *res=desn_int(*x, *mu, *omega, *alpha, *tau);
  return;
}
// transformed univariate extended skew-norm pdf:
void desn_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  double *par, a=0, mu=0, omega=0, alpha=0, tau=0, y=0;
  //set variables:
  par=(double *) fdata;
  a=par[0];
  mu=par[1];
  omega=par[2];
  alpha=par[3];
  tau=par[4];
  // computed the function
  y=a+x[0]/(1+x[0]);
  fval[0]=desn_int(y,mu,omega,alpha,tau)*1/pow(1+x[0],2);
  return;
}
// univariate extended skew-normal cumulate distribution function
void pesn(double *xmin, double *xmax, double *par, double *val, double *err)
{
  adapt_integrate(1, desn_t, par, 1, xmin, xmax, 0, 0, 1e-8, val, err);
  return;
}
// bivariate normal probability density function:
double dmn_int(double x[2], double rho, double std)
{
  double omr2=0, res=0;
  //set variables:
  omr2=1-rho*rho;
  //computes the bivariate normal density:
  res=0.5*exp(-0.5*(x[0]*x[0]-2*rho*x[0]*x[1]+x[1]*x[1])/omr2)/M_PI/std/sqrt(omr2);
  return res;
}
// bivariate extended skew-normal probability density function:
double dmesn_int(double x[2], double mu[2], double omega[4], double alpha[2], double tau)
{
  double res=0, rho=0, std=0, z[2];
  //set variables:
  z[0]=(x[0]-mu[0])/sqrt(omega[0]);
  z[1]=(x[1]-mu[1])/sqrt(omega[3]);
  std=sqrt(omega[0]*omega[3]);
  rho=omega[1]/std;
  //computes the multivariate extended skew-normal density:
  res=dmn_int(z,rho,std)*pnorm((alpha[0]*z[0]+alpha[1]*z[1])+tau,0,1,1,0)/
    pnorm(tau/sqrt(1+alpha[0]*alpha[0]+2*rho*alpha[0]*alpha[1]+alpha[1]*alpha[1]),0,1,1,0);
  return res;
}
// interface for bivariate extended skew-norm pdf:
void dmesn(double *x, double *mu, double *omega, double *alpha, double *tau, double *res)
{
  *res=dmesn_int(x, mu, omega, alpha, *tau);
  return;
}
// transformed bivariate extended skew-norm pdf:
void dmesn_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  double *par, a[2], mu[2], omega[4], alpha[2], y[2], tau=0;
  //set variables:
  par=(double *) fdata;
  a[0]=par[0];a[1]=par[1];
  mu[0]=par[2];mu[1]=par[3];
  omega[0]=par[4];omega[1]=par[5];
  omega[2]=par[6];omega[3]=par[7];
  alpha[0]=par[8];alpha[1]=par[9];
  tau=par[10];
  // computed the function
  y[0]=a[0]+x[0]/(1+x[0]);
  y[1]=a[1]+x[1]/(1+x[1]);
  fval[0]=dmesn_int(y,mu,omega,alpha,tau)*1/pow(1+x[0],2)*1/pow(1+x[1],2);
  return;
}
// bivariate extended skew-normal cumulate distribution function
void pmesn(double *xmin, double *xmax, double *par, double *val, double *err)
{
  adapt_integrate(1, dmesn_t, par, 2, xmin, xmax, 0, 0, 1e-7, val, err);
  return;
}
// trivariate normal probability density function:
double dmn_int3(double x[3], double rho[3], double std[3])
{
    double det=0, x2[3], var[3], res=0;
    var[0]=std[0]*std[0]; var[1]=std[1]*std[1]; var[2]=std[2]*std[2];
    x2[0]=x[0]*x[0]; x2[1]=x[1]*x[1]; x2[2]=x[2]*x[2];
    //determinant:
    det=1-rho[0]*rho[0]-rho[1]*rho[1]-rho[2]*rho[2]+2*rho[0]*rho[1]*rho[2];
    //computes the bivariate normal density:
    res= exp(-0.5*(x2[0]/var[0]*(1-rho[2]*rho[2]) + x2[1]/var[1]*(1-rho[1]*rho[1]) + x2[2]/var[2]*(1-rho[0]*rho[0]) + 2*x[0]/std[0]*x[1]/std[1]*(rho[1]*rho[2]-rho[0]) + 2*x[0]/std[0]*x[2]/std[2]*(rho[0]*rho[2]-rho[1]) +2*x[1]/std[1]*x[2]/std[2]*(rho[0]*rho[1]-rho[2]))/det)/(pow(M_PI,1.5)*pow(2,1.5)*sqrt(fabs(det))*std[0]*std[1]*std[2]);
    return res;
}
// trivariate extended skew-normal probability density function:
double dmesn_int3(double x[3], double mu[3], double omega[9], double alpha[3], double tau)
{
    double res=0, rho[3], std[3], z[3];
    //set variables:
    std[0]=sqrt(omega[0]);
    std[1]=sqrt(omega[4]);
    std[2]=sqrt(omega[8]);
    rho[0]=omega[1]/(std[0]*std[1]);
    rho[1]=omega[2]/(std[0]*std[2]);
    rho[2]=omega[5]/(std[1]*std[2]);
    z[0]=(x[0]-mu[0])/std[0];
    z[1]=(x[1]-mu[1])/std[1];
    z[2]=(x[2]-mu[2])/std[2];
    //computes the multivariate extended skew-normal density:
    res=dmn_int3(x,rho,std)*pnorm((alpha[0]*z[0]+alpha[1]*z[1]+alpha[2]*z[2])+tau,0,1,1,0)/
    pnorm(tau/sqrt(1+alpha[0]*alpha[0]+alpha[1]*alpha[1]+alpha[2]*alpha[2]+2*rho[0]*alpha[0]*alpha[1]+2*rho[1]*alpha[0]*alpha[2]+2*rho[2]*alpha[1]*alpha[2]),0,1,1,0);
    return res;
}

// interface for trivariate extended skew-norm pdf:
void dmesn3(double *x, double *mu, double *omega, double *alpha, double *tau, double *res)
{
    *res=dmesn_int3(x, mu, omega, alpha, *tau);
    return;
}
// transformed trivariate extended skew-norm pdf:
void dmesn_t3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    double *par, a[3], mu[3], omega[9], alpha[3], y[3], tau=0;
    //set variables:
    par=(double *) fdata;
    a[0]=par[0];a[1]=par[1];a[2]=par[2];
    mu[0]=par[3];mu[1]=par[4];mu[2]=par[5];
    omega[0]=par[6];omega[1]=par[7];omega[2]=par[8];
    omega[3]=par[9];omega[4]=par[10];omega[5]=par[11];
    omega[6]=par[12];omega[7]=par[13];omega[8]=par[14];
    alpha[0]=par[15];alpha[1]=par[16];alpha[2]=par[17];
    tau=par[18];
    // computed the function
    y[0]=a[0]+x[0]/(1+x[0]);
    y[1]=a[1]+x[1]/(1+x[1]);
    y[2]=a[2]+x[2]/(1+x[2]);
    fval[0]=dmesn_int3(y,mu,omega,alpha,tau)*1/pow(1+x[0],2)*1/pow(1+x[1],2)*1/pow(1+x[2],2);
    return;
}
// check transformed 3d e.s.n pdf
void dmesn_T(unsigned *ndim, const double *x, void *fdata, unsigned *fdim, double *fval)
{
    dmesn_t3(*ndim, x, fdata, *fdim, fval);
}

// trivariate extended skew-normal cumulate distribution function
void pmesn3(double *xmin, double *xmax, double *par, double *val, double *err)
{
    adapt_integrate(1, dmesn_t3, par, 3, xmin, xmax, 0, 0, 1e-8, val, err);
    return;
}
// internal bivariate extended skew-normal cumulate density function:
double pmesn_int(double *par)
{
    double res=0;
    double xmin[2]={-1,-1}, xmax[2]={0,0}, *val=&res, *err;
    // set the result and error variables:
    /* val=(double *) malloc(1*sizeof(double)); */
    err=(double *) malloc(1*sizeof(double));
    // compute the cumulative distribution function:
    adapt_integrate(1, dmesn_t, par, 2, xmin, xmax, 0, 0, 1e-8, val, err);
    res=val[0];
    free(err);
    return res;
}
// internal trivariate extended skew-normal cumulate density function:
double pmesn_int3(double *par)
{
    double res=0;
    double xmin[3]={-1,-1,-1}, xmax[3]={0,0,0}, *val=&res, *err;
    // set the result and error variables:
    /* val=(double *) malloc(1*sizeof(double)); */
    err=(double *) malloc(1*sizeof(double));
    // compute the cumulative distribution function:
    adapt_integrate(1, dmesn_t3, par, 3, xmin, xmax, 0, 0, 1e-8, val, err);
    res=val[0];
    free(err);
    return res;
}


/* END ------ SKEW-NORMAL FUNCTIONS -------------- */



/* START ------ SKEW-T FUNCTIONS -------------- */
// univariate extended skew-t probability density function:
double dest_int(double x, double mu, double omega, double nu, double alpha, double tau)
{
  double res=0, y=(x-mu)/omega, nu1=nu+1;
  //computes the univariate extended skew-t density:
  res=dt(y,nu,0)/omega*pt((alpha*y+tau)*sqrt(nu1/(nu+y*y)),nu1,1,0)/pt(tau/sqrt(1+alpha*alpha),nu,1,0);
  return res;
}
// interface for univariate extended skew-t pdf:
void dest(double *x, double *mu, double *omega, double *nu, double *alpha, double *tau, double *res)
{
  *res=dest_int(*x, *mu, *omega, *nu, *alpha, *tau);
  return;
}
// transformed univariate extended skew-t pdf:
void dest_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  double *par, a=0, mu=0, omega=0, nu=0, alpha=0, tau=0, y=0;
  //set variables:
  par=(double *) fdata;
  a=par[0];
  mu=par[1];
  omega=par[2];
  nu=par[3];
  alpha=par[4];
  tau=par[5];
  // computed the function
  y=a+x[0]/(1+x[0]);
  fval[0]=dest_int(y,mu,omega,nu,alpha,tau)*1/pow(1+x[0],2);
  return;
}
// univariate extended skew-t cumulate distribution function
void pest(double *xmin, double *xmax, double *par, double *val, double *err)
{
  adapt_integrate(1, dest_t, par, 1, xmin, xmax, 0, 0, 1e-13, val, err);
  return;
}
// internal extended skew-t cumulate density function
double pest_int(double *par)
{
  double res=0;
  double *xmin, *xmax, *val = &res, *err;
  // set the result and error variables:
  //val=(double *) malloc(1*sizeof(double));
  err=(double *) malloc(1*sizeof(double));
  xmin=(double *) malloc(1*sizeof(double));
  xmax=(double *) malloc(1*sizeof(double));
  xmin[0]=-1; xmax[0]=0;
  // compute the cumulative distribution function:
  adapt_integrate(1, dest_t, par, 1, xmin, xmax, 0, 0, 1e-13, val, err);
  res=val[0];
  free(err);
  free(xmin);
  free(xmax);
  return res;
}
// bivariate t probability density function:
double dmt_int(double nu, double omr2, double Q, double std)
{
  double res=0, hnu=0.5*nu, hnu1=hnu+1;
  //computes the bivariate t density:
  res=gammafn(hnu1)/nu/M_PI/gammafn(hnu)/std/sqrt(omr2)*
    pow(1+Q/nu,-hnu1);
  return res;
}
// bivariate extended skew-t probability density function:
double dmest_int(double x[2], double mu[2], double omega[4], double nu, double alpha[2], double tau)
{
  double omr2=0, nu2=nu+2, Q=0, res=0, rho=0, std=0, z[2];
  //set variables:
  z[0]=(x[0]-mu[0])/sqrt(omega[0]);
  z[1]=(x[1]-mu[1])/sqrt(omega[3]);
  std=sqrt(omega[0]*omega[3]);
  rho=omega[1]/std;
  omr2=1-rho*rho;
  Q=(z[0]*z[0]-2*rho*z[0]*z[1]+z[1]*z[1])/omr2;
  //computes the multivariate extended skew-t density:
  res=dmt_int(nu,omr2,Q,std)*pt((alpha[0]*z[0]+alpha[1]*z[1]+tau)*sqrt(nu2/(nu+Q)),nu2,1,0)/
    pt(tau/sqrt(1+alpha[0]*alpha[0]+2*rho*alpha[0]*alpha[1]+alpha[1]*alpha[1]),nu,1,0);
  return res;
}
// interface for bivariate extended skew-t pdf:
void dmest(double *x, double *mu, double *omega, double *nu, double *alpha, double *tau, double *res)
{
  *res=dmest_int(x, mu, omega, *nu, alpha, *tau);
  return;
}
// transformed bivariate extended skew-t pdf:
void dmest_t(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  double *par, a[2], mu[2], omega[4], nu=0, alpha[2], y[2], tau=0;
  //set variables:
  par=(double *) fdata;
  a[0]=par[0];a[1]=par[1];
  mu[0]=par[2];mu[1]=par[3];
  omega[0]=par[4];omega[1]=par[5];
  omega[2]=par[6];omega[3]=par[7];
  nu=par[8];
  alpha[0]=par[9];alpha[1]=par[10];
  tau=par[11];
  // computed the function
  y[0]=a[0]+x[0]/(1+x[0]);
  y[1]=a[1]+x[1]/(1+x[1]);
  fval[0]=dmest_int(y,mu,omega,nu,alpha,tau)*1/pow(1+x[0],2)*1/pow(1+x[1],2);
  return;
}
// bivariate extended skew-t cumulate distribution function;
void pmest(double *xmin, double *xmax, double *par, double *val, double *err)
{
  adapt_integrate(1, dmest_t, par, 2, xmin, xmax, 0, 0, 1e-8, val, err);
  return;
}
// trivariate t probability density function:
double dmt_int3(double nu, double Q, double det)
{
  double res=0, hnu=0.5*nu, hnu1=hnu+1.5;
  //computes the bivariate t density:
  res=gammafn(hnu1)*pow(nu*M_PI,-1.5)/gammafn(hnu)/sqrt(det)*pow(1+Q/nu,-hnu1);
  return res;
}
// trivariate extended skew-t probability density function:
double dmest_int3(double x[3], double mu[3], double omega[9], double nu, double alpha[3], double tau)
{
  double det=0, nu3=nu+3, Q=0, res=0, z[3];
  double rho12=0, rho122=0, rho13=0, rho132=0, rho23=0, rho232=0;
  double so1=sqrt(omega[0]);
  double so2=sqrt(omega[4]);
  double so3=sqrt(omega[8]);
  //set variables:
  z[0]=(x[0]-mu[0])/so1;
  z[1]=(x[1]-mu[1])/so2;
  z[2]=(x[2]-mu[2])/so3;
  //define correlations:
  rho12=omega[1]/(so1*so2);
  rho13=omega[2]/(so1*so3);
  rho23=omega[5]/(so2*so3);
  rho122=rho12*rho12;
  rho132=rho13*rho13;
  rho232=rho23*rho23;
  // define determinant
  det=1+2*rho12*rho13*rho23-rho122-rho132-rho232;
  // compute the quadratic form
  Q=(z[0]*z[0]*(1-rho232)+z[1]*z[1]*(1-rho132)+z[2]*z[2]*(1-rho122)+
  2*z[0]*z[1]*(rho13*rho23-rho12)+2*z[0]*z[2]*(rho12*rho23-rho13)+2*z[1]*z[2]*(rho12*rho13-rho23))/det;
  //computes the multivariate extended skew-t density:
  res=dmt_int3(nu,Q,det*omega[0]*omega[4]*omega[8])*pt((alpha[0]*z[0]+alpha[1]*z[1]+alpha[2]*z[2]+tau)*sqrt(nu3/(nu+Q)),nu3,1,0)/
  pt(tau/sqrt(1+alpha[0]*alpha[0]+alpha[1]*alpha[1]+alpha[2]*alpha[2]+2*rho12*alpha[0]*alpha[1]+
  2*rho13*alpha[0]*alpha[2]+2*rho23*alpha[1]*alpha[2]),nu,1,0);
  return res;
}
// transformed trivariate extended skew-t pdf:
void dmest_t3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  double *par, a[3], mu[3], omega[9], nu=0, alpha[3], y[3], tau=0;
  //set variables:

  par=(double *) fdata;
  a[0]=par[0];a[1]=par[1];a[2]=par[2];
  mu[0]=par[3];mu[1]=par[4];mu[2]=par[5];
  omega[0]=par[6];omega[1]=par[7];
  omega[2]=par[8];omega[3]=par[9];
  omega[4]=par[10];omega[5]=par[11];
  omega[6]=par[12];omega[7]=par[13];
  omega[8]=par[14];
  nu=par[15];
  alpha[0]=par[16];alpha[1]=par[17];
  alpha[2]=par[18];
  tau=par[19];
  // computed the function
  y[0]=a[0]+x[0]/(1+x[0]);
  y[1]=a[1]+x[1]/(1+x[1]);
  y[2]=a[2]+x[2]/(1+x[2]);
  fval[0]=dmest_int3(y,mu,omega,nu,alpha,tau)*1/pow(1+x[0],2)*1/pow(1+x[1],2)*1/pow(1+x[2],2);
  return;
}
// check transformed 3d e.s.t pdf
void dmest_T(unsigned *ndim, const double *x, void *fdata, unsigned *fdim, double *fval)
{
    dmest_t3(*ndim, x, fdata, *fdim, fval);
}
// interface for trivariate extended skew-t pdf:
void dmest3(double *x, double *mu, double *omega, double *nu, double *alpha, double *tau, double *res)
{
  *res=dmest_int3(x, mu, omega, *nu, alpha, *tau);
  return;
}
// trivariate extended skew-t cumulate distribution function:
void pmest3(double *xmin, double *xmax, double *par, double *val, double *err)
{
  adapt_integrate(1, dmest_t3, par, 3, xmin, xmax, 0, 0, 1e-8, val, err);
  return;
}


// internal bivariate extended skew-t cumulate density function:
double pmest_int(double *par)
{
  double res=0;
  double xmin[2]={-1,-1}, xmax[2]={0,0}, *val=&res, *err;
  // set the result and error variables:
  /* val=(double *) malloc(1*sizeof(double)); */
  err=(double *) malloc(1*sizeof(double));
  // compute the cumulative distribution function:
  adapt_integrate(1, dmest_t, par, 2, xmin, xmax, 0, 0, 1e-8, val, err);
  res=val[0];
  free(err);
  return res;
}
// internal trivariate extended skew-t cumulate density function:
double pmest_int3(double *par)
{
  double res=0;
  double xmin[3]={-1,-1,-1}, xmax[3]={0,0,0}, *val=&res, *err;
  // set the result and error variables:
  /* val=(double *) malloc(1*sizeof(double)); */
  err=(double *) malloc(1*sizeof(double));
  // compute the cumulative distribution function:
  adapt_integrate(1, dmest_t3, par, 3, xmin, xmax, 0, 0, 1e-8, val, err);
  res=val[0];
  free(err);
  return res;
}

/* END ------ SKEW-T FUNCTIONS -------------- */

/* START ------ HUSLER-REISS FUNCTIONS -------------- */

// bivariate Husler-Reiss cumulative distribution function:
void pHuslerReiss(double *x, double *lambda, double *res)
{
  double y[2], lxy, a=*lambda;
  lxy=log(x[1]/x[0])/(2*a);
  y[0]=a+lxy;
  y[1]=a-lxy;
  
  res[0]=exp(-pnorm(y[0],0,1,1,0)/x[0]-pnorm(y[1],0,1,1,0)/x[1]);
    
  return;
}
// bivariate Husler-Reiss probability density function:
void dHuslerReiss(double *x, double *lambda, double *res)
{
  double y[2], a=*lambda;
  y[0]=x[0]; y[1]=x[1];
  *res=HuslerReiss(a,y);
    
  return;
}
/* END ------ HUSLER-REISS FUNCTIONS -------------- */

/* START ------ EXTREMAL SKEW-T FUNCTIONS -------------- */

// bivariate extremal skew-t cumulative distribution function:
void pmextst(double *x, double *omega, double *nu, double *alpha, double *res)
{
  double nu1=*nu+1, omo2=1-pow(omega[0],2), snu1=sqrt(nu1), nui=1/ *nu;
  double somo2=sqrt(omo2), alphab[2], alphat[2], alphas[2], taus[2];
  double xb[2], y[2], *par1, *par2;

  par1=(double *) malloc(6*sizeof(double));
  par2=(double *) malloc(6*sizeof(double));
  //defines the quantities alpha tilda:
  alphat[0]=alpha[0]+alpha[1]*omega[0];
  alphat[1]=alpha[1]+alpha[0]*omega[0];
  //defines the quantities alpha bar:
  alphab[0]=alphat[0]/sqrt(1+pow(alpha[1],2)*omo2);
  alphab[1]=alphat[1]/sqrt(1+pow(alpha[0],2)*omo2);
  //defines the quantities alpha star:
  alphas[0]=alpha[0]*somo2;
  alphas[1]=alpha[1]*somo2;
  // define conditional extended parameters:
  taus[0]=snu1*alphat[0];
  taus[1]=snu1*alphat[1];
  // define quantiles:
  xb[0]=pow((x[1]*pt(snu1*alphab[1],*nu,1,0))/(x[0]*pt(snu1*alphab[0],*nu,1,0)),nui);
  xb[1]=pow((x[0]*pt(snu1*alphab[0],*nu,1,0))/(x[1]*pt(snu1*alphab[1],*nu,1,0)),nui);
  y[0]=snu1*(xb[0]-omega[0])/somo2;
  y[1]=snu1*(xb[1]-omega[0])/somo2;
  // define the extended skew-t parameters:
  par1[0]=y[0];par1[1]=0;par1[2]=1;par1[3]=nu1;
  par1[4]=alphas[1];par1[5]=taus[0];
  par2[0]=y[1];par2[1]=0;par2[2]=1;par2[3]=nu1;
  par2[4]=alphas[0];par2[5]=taus[1];
  // computes the trivariate Pickand dependence function for the extremal skew-t model:
  res[0]=exp(-pest_int(par1)/x[0]-pest_int(par2)/x[1]);
    free(par1);
    free(par2);
  return;
}
// bivariate extremal skew-t probability density function:
void dmextst(double *x, double *omega, double *nu, double *alpha, double *res)
{
  double nu1=*nu+1, omo2=1-pow(omega[0],2), snu1=sqrt(nu1), nui=1/ *nu;
  double somo2=sqrt(omo2), alphab[2], alphat[2], alphas[2], taus[2];
  double beta=snu1/somo2, V=0.0, dx1V=0.0, dx2V=0.0, d2V=0.0, cdf1=0.0;
  double cdf2=0.0, pdf1=0.0, pdf2=0.0, xb[2], xs[2], y[2], ys[2];
  double nu2=nu1+1, pt1=0.0, pt2=0.0, u[2], t[2], v[2], w[2], *par1;
  double *par2, scdf1=0.0, scdf2=0.0, tpdf1=0.0, tpdf2=0.0, pw=0.5*(*nu+3);
  double sp=sqrt(M_PI*nu2), ga1=gammafn(pw), ga2=gammafn(0.5*nu2);

  par1=(double *) malloc(6*sizeof(double));
  par2=(double *) malloc(6*sizeof(double));
  xs[0]=pow(x[0],2);xs[1]=pow(x[1],2);
  //defines the quantities alpha tilda:
  alphat[0]=alpha[0]+alpha[1]*omega[0];
  alphat[1]=alpha[1]+alpha[0]*omega[0];
  //defines the quantities alpha bar:
  alphab[0]=alphat[0]/sqrt(1+pow(alpha[1],2)*omo2);
  alphab[1]=alphat[1]/sqrt(1+pow(alpha[0],2)*omo2);
  //defines the quantities alpha star:
  alphas[0]=alpha[0]*somo2;
  alphas[1]=alpha[1]*somo2;
  // define conditional extended parameters:
  taus[0]=snu1*alphat[0];
  taus[1]=snu1*alphat[1];
  // define quantiles:
  xb[0]=pow((x[1]*pt(snu1*alphab[1],*nu,1,0))/(x[0]*pt(snu1*alphab[0],*nu,1,0)),nui);
  xb[1]=pow((x[0]*pt(snu1*alphab[0],*nu,1,0))/(x[1]*pt(snu1*alphab[1],*nu,1,0)),nui);
  y[0]=beta*(xb[0]-omega[0]);
  y[1]=beta*(xb[1]-omega[0]);
  ys[0]=pow(y[0],2);
  ys[1]=pow(y[1],2);
  // define the extended skew-t parameters:
  par1[0]=y[0];par1[1]=0;par1[2]=1;par1[3]=nu1;
  par1[4]=alphas[1];par1[5]=taus[0];
  par2[0]=y[1];par2[1]=0;par2[2]=1;par2[3]=nu1;
  par2[4]=alphas[0];par2[5]=taus[1];
  // computes the univariate extended skew-t cdfs and pdfs:
  cdf1=pest_int(par1);
  cdf2=pest_int(par2);
  pdf1=dest_int(y[0],0,1,nu1,alphas[1],taus[0]);
  pdf2=dest_int(y[1],0,1,nu1,alphas[0],taus[1]);
  t[0]=nu1+ys[0];
  t[1]=nu1+ys[1];
  u[0]=alphas[1]*y[0]+taus[0];
  u[1]=alphas[0]*y[1]+taus[1];
  v[0]=1+ys[0]/nu1;
  v[1]=1+ys[1]/nu1;
  w[0]=sqrt(nu2/t[0]);
  w[1]=sqrt(nu2/t[1]);
  pt1=pt(u[0]*w[0],nu2,1,0);
  pt2=pt(u[1]*w[1],nu2,1,0);
  scdf1=cdf1/x[0];
  scdf2=cdf2/x[1];
  tpdf1=pdf1*xb[0];
  tpdf2=pdf2*xb[1];
  // computes the useful quantities for the extremal skew-t cdf:
  V=-scdf1-scdf2;
  dx1V=(scdf1+beta*(tpdf1/x[0]-tpdf2/x[1])/ *nu)/x[0];
  dx2V=(scdf2+beta*(tpdf2/x[1]-tpdf1/x[0])/ *nu)/x[1];

  d2V=tpdf1*beta/(xs[0]*x[1]* *nu)*
    (1+(1+beta*xb[0]*(ga1*pow(1+pow(u[0],2)/t[0],-pw)*w[0]*
  	      (alphas[1]-u[0]*y[0]/t[0])/(sp*ga2*pt1)-y[0]*nu2/(nu1*v[0]))) / *nu)+
    tpdf2*beta/(xs[1]*x[0]* *nu)*
    (1+(1+beta*xb[1]*(ga1*pow(1+pow(u[1],2)/t[1],-pw)*w[1]*
		      (alphas[0]-u[1]*y[1]/t[1])/(sp*ga2*pt2)-y[1]*nu2/(nu1*v[1]))) / *nu);
  res[0]=exp(V)*(dx1V*dx2V+d2V);
    free(par1);
    free(par2);
  return;
}
// bivariate extremal skew-t probability density function:
double dmextst_int(double *x, double *omega, double *nu, double *alpha)
{
  double nu1=*nu+1, omo2=1-pow(omega[0],2), snu1=sqrt(nu1), nui=1/ *nu, res=0.0;
  double somo2=sqrt(omo2), alphab[2], alphat[2], alphas[2], taus[2];
  double beta=snu1/somo2, V=0.0, dx1V=0.0, dx2V=0.0, d2V=0.0, cdf1=0.0;
  double cdf2=0.0, pdf1=0.0, pdf2=0.0, xb[2], xs[2], y[2], ys[2];
  double nu2=nu1+1, pt1=0.0, pt2=0.0, u[2], t[2], v[2], w[2], *par1;
  double *par2, scdf1=0.0, scdf2=0.0, tpdf1=0.0, tpdf2=0.0, pw=0.5*(*nu+3);
  double sp=sqrt(M_PI*nu2), ga1=gammafn(pw), ga2=gammafn(0.5*nu2);

  par1=(double *) malloc(6*sizeof(double));
  par2=(double *) malloc(6*sizeof(double));
  xs[0]=pow(x[0],2);xs[1]=pow(x[1],2);
  //defines the quantities alpha tilda:
  alphat[0]=alpha[0]+alpha[1]*omega[0];
  alphat[1]=alpha[1]+alpha[0]*omega[0];
  //defines the quantities alpha bar:
  alphab[0]=alphat[0]/sqrt(1+pow(alpha[1],2)*omo2);
  alphab[1]=alphat[1]/sqrt(1+pow(alpha[0],2)*omo2);
  //defines the quantities alpha star:
  alphas[0]=alpha[0]*somo2;
  alphas[1]=alpha[1]*somo2;
  // define conditional extended parameters:
  taus[0]=snu1*alphat[0];
  taus[1]=snu1*alphat[1];
  // define quantiles:
  xb[0]=pow((x[1]*pt(snu1*alphab[1],*nu,1,0))/(x[0]*pt(snu1*alphab[0],*nu,1,0)),nui);
  xb[1]=pow((x[0]*pt(snu1*alphab[0],*nu,1,0))/(x[1]*pt(snu1*alphab[1],*nu,1,0)),nui);
  y[0]=beta*(xb[0]-omega[0]);
  y[1]=beta*(xb[1]-omega[0]);
  ys[0]=pow(y[0],2);
  ys[1]=pow(y[1],2);
  // define the extended skew-t parameters:
  par1[0]=y[0];par1[1]=0;par1[2]=1;par1[3]=nu1;
  par1[4]=alphas[1];par1[5]=taus[0];
  par2[0]=y[1];par2[1]=0;par2[2]=1;par2[3]=nu1;
  par2[4]=alphas[0];par2[5]=taus[1];
  // computes the univariate extended skew-t cdfs and pdfs:
  cdf1=pest_int(par1);
  cdf2=pest_int(par2);
  pdf1=dest_int(y[0],0,1,nu1,alphas[1],taus[0]);
  pdf2=dest_int(y[1],0,1,nu1,alphas[0],taus[1]);
  t[0]=nu1+ys[0];
  t[1]=nu1+ys[1];
  u[0]=alphas[1]*y[0]+taus[0];
  u[1]=alphas[0]*y[1]+taus[1];
  v[0]=1+ys[0]/nu1;
  v[1]=1+ys[1]/nu1;
  w[0]=sqrt(nu2/t[0]);
  w[1]=sqrt(nu2/t[1]);
  pt1=pt(u[0]*w[0],nu2,1,0);
  pt2=pt(u[1]*w[1],nu2,1,0);
  scdf1=cdf1/x[0];
  scdf2=cdf2/x[1];
  tpdf1=pdf1*xb[0];
  tpdf2=pdf2*xb[1];
  // computes the useful quantities for the extremal skew-t cdf:
  V=-scdf1-scdf2;
  dx1V=(scdf1+beta*(tpdf1/x[0]-tpdf2/x[1])/ *nu)/x[0];
  dx2V=(scdf2+beta*(tpdf2/x[1]-tpdf1/x[0])/ *nu)/x[1];

  d2V=tpdf1*beta/(xs[0]*x[1]* *nu)*
    (1+(1+beta*xb[0]*(ga1*pow(1+pow(u[0],2)/t[0],-pw)*w[0]*
		      (alphas[1]-u[0]*y[0]/t[0])/(sp*ga2*pt1)-y[0]*nu2/(nu1*v[0]))) / *nu)+
    tpdf2*beta/(xs[1]*x[0]* *nu)*
    (1+(1+beta*xb[1]*(ga1*pow(1+pow(u[1],2)/t[1],-pw)*w[1]*
		      (alphas[0]-u[1]*y[1]/t[1])/(sp*ga2*pt2)-y[1]*nu2/(nu1*v[1]))) / *nu);
  res=exp(V)*(dx1V*dx2V+d2V);
  free(par1);
  free(par2);
  return res;
}

/* END ------ EXTREMAL SKEW-T FUNCTIONS -------------- */
