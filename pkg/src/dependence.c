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




/* START ------ EXTREMAL DEPENDENCE FUNCTIONS -------------- */


// Upper tail dependence coefficient for the extremal skew-t model:
void chistup(double *omega, double *nu, double *alpha, double *res)
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
  xb[0]=pt(snu1*alphab[1],nu1,1,0)/pt(snu1*alphab[0],nu1,1,0);
  xb[1]=pt(snu1*alphab[0],nu1,1,0)/pt(snu1*alphab[1],nu1,1,0);
  y[0]=snu1*(pow(xb[0],nui)-omega[0])/somo2;
  y[1]=snu1*(pow(xb[1],nui)-omega[0])/somo2;
  // define the extended skew-t parameters:
  par1[0]=y[0];par1[1]=0;par1[2]=1;par1[3]=nu1;
  par1[4]=alphas[1];par1[5]=taus[0];
  par2[0]=y[1];par2[1]=0;par2[2]=1;par2[3]=nu1;
  par2[4]=alphas[0];par2[5]=taus[1];
  // computes the upper tail dependence coefficient for the extremal skew-t model:
  res[0]=2-pest_int(par1)-pest_int(par2);
    free(par1);
    free(par2);
  return;
}
// Lower tail dependence coefficient for the extremal skew-t model:
void chistlo(double *omega, double *nu, double *alpha, double *res)
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
  xb[0]=pt(-snu1*alphab[1],nu1,1,0)/pt(-snu1*alphab[0],nu1,1,0);
  xb[1]=pt(-snu1*alphab[0],nu1,1,0)/pt(-snu1*alphab[1],nu1,1,0);
  y[0]=snu1*(-pow(xb[0],nui)+omega[0])/somo2;
  y[1]=snu1*(-pow(xb[1],nui)+omega[0])/somo2;
  // define the extended skew-t parameters:
  par1[0]=y[0];par1[1]=0;par1[2]=1;par1[3]=nu1;
  par1[4]=alphas[1];par1[5]=-taus[0];
  par2[0]=y[1];par2[1]=0;par2[2]=1;par2[3]=nu1;
  par2[4]=alphas[0];par2[5]=-taus[1];
  // computes the tail dependence coefficient for the extremal skew-t model:
  res[0]=pest_int(par1)+pest_int(par2);
    free(par1);
    free(par2);
  return;
}
// Bivariate Pickand dependence function for the extremal skew-t model:
void bivpkst(double *x, double *omega, double *nu, double *alpha, double *res)
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
  xb[0]=(x[0]*pt(snu1*alphab[1],nu1,1,0))/((1-x[0])*pt(snu1*alphab[0],nu1,1,0));
  xb[1]=((1-x[0])*pt(snu1*alphab[0],nu1,1,0))/(x[0]*pt(snu1*alphab[1],nu1,1,0));
  y[0]=snu1*(pow(xb[0],nui)-omega[0])/somo2;
  y[1]=snu1*(pow(xb[1],nui)-omega[0])/somo2;
  // define the extended skew-t parameters:
  par1[0]=y[0];par1[1]=0;par1[2]=1;par1[3]=nu1;
  par1[4]=alphas[1];par1[5]=taus[0];
  par2[0]=y[1];par2[1]=0;par2[2]=1;par2[3]=nu1;
  par2[4]=alphas[0];par2[5]=taus[1];
  // computes the trivariate Pickand dependence function for the extremal skew-t model:
  res[0]=x[0]*pest_int(par1)+(1-x[0])*pest_int(par2);
    free(par1);
    free(par2);
  return;
}
// Trivariate Pickand dependence function for the extremal skew-t model:
void trivpkst(double *x, double *omega, double *nu, double *alpha, double *res)
{
  double nu1=*nu+1, omo2[3]={1-pow(omega[1],2),1-pow(omega[2],2),1-pow(omega[5],2)};
  double snu1=sqrt(nu1), somo2[3], alphab[3], alphat[3], alphas[6], taus[3];
  double nui=1/ *nu, psi=0, xb[6], y[6], psi1[4], psi2[4], psi3[4];
  double omegac[3], *par1, *par2, *par3;

  par1=(double *) malloc(12*sizeof(double));
  par2=(double *) malloc(12*sizeof(double));
  par3=(double *) malloc(12*sizeof(double));
  //matrix of sqrts:
  somo2[0]=sqrt(omo2[0]);
  somo2[1]=sqrt(omo2[1]);
  somo2[2]=sqrt(omo2[2]);
  //matrix of conditional variances:
  omegac[0]=omega[5]-omega[1]*omega[2];
  omegac[1]=omega[2]-omega[1]*omega[5];
  omegac[2]=omega[1]-omega[2]*omega[5];
  //defines the quantities alpha tilda:
  alphat[0]=alpha[0]+alpha[1]*omega[1]+alpha[2]*omega[2];
  alphat[1]=alpha[1]+alpha[0]*omega[1]+alpha[2]*omega[5];
  alphat[2]=alpha[2]+alpha[0]*omega[2]+alpha[1]*omega[5];
  //defines the quantities alpha bar:
  alphab[0]=alphat[0]/sqrt(1+pow(alpha[1],2)*omo2[0]+2*alpha[1]*alpha[2]*omegac[0]+pow(alpha[2],2)*omo2[1]);
  alphab[1]=alphat[1]/sqrt(1+pow(alpha[0],2)*omo2[0]+2*alpha[0]*alpha[2]*omegac[1]+pow(alpha[2],2)*omo2[2]);
  alphab[2]=alphat[2]/sqrt(1+pow(alpha[0],2)*omo2[1]+2*alpha[0]*alpha[1]*omegac[2]+pow(alpha[1],2)*omo2[2]);
  //defines the quantities alpha star:
  alphas[0]=alpha[1]*somo2[0];
  alphas[1]=alpha[2]*somo2[1];
  alphas[2]=alpha[0]*somo2[0];
  alphas[3]=alpha[2]*somo2[2];
  alphas[4]=alpha[0]*somo2[1];
  alphas[5]=alpha[1]*somo2[2];
  // define conditional extended parameters:
  taus[0]=snu1*alphat[0];
  taus[1]=snu1*alphat[1];
  taus[2]=snu1*alphat[2];
  // define quantiles:
  xb[0]=(x[0]*pt(snu1*alphab[1],nu1,1,0))/(x[1]*pt(snu1*alphab[0],nu1,1,0));
  xb[1]=(x[0]*pt(snu1*alphab[2],nu1,1,0))/(x[2]*pt(snu1*alphab[0],nu1,1,0));
  y[0]=snu1*(pow(xb[0],nui)-omega[1])/somo2[0];
  y[1]=snu1*(pow(xb[1],nui)-omega[2])/somo2[1];

  xb[2]=(x[1]*pt(snu1*alphab[0],nu1,1,0))/(x[0]*pt(snu1*alphab[1],nu1,1,0));
  xb[3]=(x[1]*pt(snu1*alphab[2],nu1,1,0))/(x[2]*pt(snu1*alphab[1],nu1,1,0));
  y[2]=snu1*(pow(xb[2],nui)-omega[1])/somo2[0];
  y[3]=snu1*(pow(xb[3],nui)-omega[5])/somo2[2];

  xb[4]=(x[2]*pt(snu1*alphab[0],nu1,1,0))/(x[0]*pt(snu1*alphab[2],nu1,1,0));
  xb[5]=(x[2]*pt(snu1*alphab[1],nu1,1,0))/(x[1]*pt(snu1*alphab[2],nu1,1,0));
  y[4]=snu1*(pow(xb[4],nui)-omega[2])/somo2[1];
  y[5]=snu1*(pow(xb[5],nui)-omega[5])/somo2[2];
  // define the conditional Psi dispersion matrices of the extended skew-t models:
  psi=omegac[0]/(somo2[0]*somo2[1]);
  psi1[0]=1;psi1[1]=psi;psi1[2]=psi;psi1[3]=1;
  psi=omegac[1]/(somo2[0]*somo2[2]);
  psi2[0]=1;psi2[1]=psi;psi2[2]=psi;psi2[3]=1;
  psi=omegac[2]/(somo2[1]*somo2[2]);
  psi3[0]=1;psi3[1]=psi;psi3[2]=psi;psi3[3]=1;
  // define the extended skew-t parameters:
  par1[0]=y[0];par1[1]=y[1];par1[2]=0;par1[3]=0;
  par1[4]=psi1[0];par1[5]=psi1[1];par1[6]=psi1[2];par1[7]=psi1[3];
  par1[8]=nu1;par1[9]=alphas[0];par1[10]=alphas[1];par1[11]=taus[0];
  par2[0]=y[2];par2[1]=y[3];par2[2]=0;par2[3]=0;
  par2[4]=psi2[0];par2[5]=psi2[1];par2[6]=psi2[2];par2[7]=psi2[3];
  par2[8]=nu1;par2[9]=alphas[2];par2[10]=alphas[3];par2[11]=taus[1];
  par3[0]=y[4];par3[1]=y[5];par3[2]=0;par3[3]=0;
  par3[4]=psi3[0];par3[5]=psi3[1];par3[6]=psi3[2];par3[7]=psi3[3];
  par3[8]=nu1;par3[9]=alphas[4];par3[10]=alphas[5];par3[11]=taus[2];
  // computes the trivariate Pickand dependence function for the extremal skew-t model:
  res[0]=x[0]*pmest_int(par1)+x[1]*pmest_int(par2)+x[2]*pmest_int(par3);
    free(par1);
    free(par2);
    free(par3);
  return;
}


/* END ------ EXTREMAL DEPENDENCE FUNCTIONS -------------- */
