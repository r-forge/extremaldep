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




/* START ------ LIKELIHOOD FUNCTIONS -------------- */



// bivariate log-likelihood for the extremal skew-t model:
void llextst(double *x, int *n, double *omega, double *nu, double *alpha, double *res)
{
  int i=0;
  double ans=0.0, y[2];
  //checks the validity of the parameters:
  if(*omega<=-.999 || *omega>=.999 || *nu<=1 || *nu>100){
    *res= -1.0e8;
    return;}
  for(i=0;i<*n;i++){
    y[0]=x[i];y[1]=x[i+*n];
    ans=dmextst_int(y, omega, nu, alpha);
    *res+=log(ans);}
  //check if the log-likelihood is finite:
  if(!R_FINITE(*res))
    *res= -1.0e8;
  return;
  }
  
  // Pairwise likelihood for the Husler-Reiss model:
double HuslerReiss(double a, double *xx)
{
  //double ao2=0.0, ax2=0.0, axy=0.0, ay2=0.0, d2V=0.0;
  double ax2=0.0, axy=0.0, ay2=0.0, d2V=0.0;
  double dx=0.0, dxV=0.0, dy=0.0, dyV=0.0, lyx=0.0, px=0.0;
  double x=xx[0],y=xx[1];
  double py=0.0, V=0.0, x2=0.0, y2=0.0, w=0.0, z=0.0;

  //ao2=0.5*a;
  axy=2*a*x*y;
  x2=pow(x,2);
  y2=pow(y,2);
  ax2=2*a*x2;
  ay2=2*a*y2;
  lyx=log(y/x)/(2*a);
  z=a+lyx;
  w=a-lyx;
  px=pnorm(z,0,1,1,0);
  py=pnorm(w,0,1,1,0);
  dx=dnorm(z,0,1,0);
  dy=dnorm(w,0,1,0);
  // Defines the likelihood components:
  V=-px/x-py/y;
  dxV=px/x2+dx/ax2-dy/axy;
  dyV=py/y2+dy/ay2-dx/axy;
  d2V=(w*dx*y+z*dy*x)/(ax2*ay2);
  return exp(V)*(dxV*dyV+d2V);
}

// log-likelihood function of the HR model
void llHRmax(double *x, double *lambda, int *n, double *res)
{
  double y[2], a=*lambda;
  int i;
  for(i=0;i<*n;i++){
    y[0]=x[i];y[1]=x[i+*n];
    *res+=log(HuslerReiss(a, y));
  }
  return;
}

double d1x_dt(double x, double df)
{
  double res=0.0;
  res=-dt(x,df,0)*(df+1)*x/
    (1+pow(x,2)/df)/df;
  return res;
}
// log-likelihood for the extremal-t model (one observation):
double ExtremalT(double *data, double df, double rho)
{
  double a=0.0, ac=0.0, aci=0.0, c=0.0, df1=df+1;
  double d1tx=0.0, d1ty=0.0, dtx=0.0, dty=0.0, d2V=0.0, dxV=0.0, dyV=0.0;
  double opdf=1+1/df, ptx=0.0, pty=0.0, res=0.0, x=data[0], x2=0.0;
  double x2d=0.0, xyd=0.0, y=data[1], y2=0.0, y2d=0.0, V=0.0, w=0.0, z=0.0;
  // Computes the log-likelihood:
  a=sqrt(df1/(1-pow(rho, 2)));
  c=pow(y/x,1/df);
  ac=a*c;
  aci=a/c;
  z=(c-rho)*a;
  w=(1/c-rho)*a;
  x2=pow(x,2);
  y2=pow(y,2);
  x2d=x2*df;
  y2d=y2*df;
  xyd=x*y*df;
  ptx=pt(z,df1,1,0);
  pty=pt(w,df1,1,0);
  dtx=dt(z,df1,0);
  dty=dt(w,df1,0);
  d1tx=d1x_dt(z,df1);
  d1ty=d1x_dt(w,df1);
  //defines the log-likelihood components:
  V=-ptx/x-pty/y;
  dxV=ptx/x2+dtx*ac/x2d-dty*aci/xyd;
  dyV=pty/y2+dty*aci/y2d-dtx*ac/xyd;
  d2V=ac*(dtx*opdf+d1tx*ac/df)/x2d/y+aci*(dty*opdf+d1ty*aci/df)/y2d/x;
  res=exp(V)*(d2V+dxV*dyV);//log-likelihood
  return res;
}

// log-likelihood for the extremal-t model (n observation):
void llETmax(double *x, double *par, int *n, double *res)
{
  double ans=0.0, lans=0.0, data[2], df=par[0], rho=par[1];
  int i=0;
  // Checks the validity of the df and correlation parameters:
  if(df<=0 || rho< -1 || rho>1){*res=LOW; return;}
  for(i=0;i<*n;i++){
    data[0]=x[i];data[1]=x[i+*n];
    ans=ExtremalT(data, df, rho);
    if(!isnan(ans)){
      lans=log(ans);
      if(R_FINITE(lans)) *res+=log(ans);
      else *res+=-1.0e5;}}
  if(!R_FINITE(*res)) *res=LOW;
  return;
}


/* END ------ LIKELIHOOD FUNCTIONS -------------- */
