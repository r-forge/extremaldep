#include "cubature.h"

void dGEV(double *x, double *loc, double *scale, double *shape, double *res)
{
    res[0] = dgev(x[0], loc[0], scale[0], shape[0]);
}

double dgev(double x, double loc, double scale, double shape)
{
  double y=0.0, res=0.0;

  y = (x - loc) / scale;

  if(shape==0)
    res = exp(-exp(-y) - y) / scale;
  else
    res = exp(-pow(fmax(1 + shape * y, 0), - 1 / shape)) *
      pow(fmax(1 + shape * y, 0), - 1 / shape - 1) / scale;

  return res;
}

void GevLogLik(double *data, int *ndata, double *par, double *res)
{
  int n=0;

  if(par[1] <= 0)
    {
      *res = LOW;
      return;
    }

  for(n = 0; n < *ndata; n++)
    *res += log(dgev(data[n], par[0], par[1], par[2]));

  if(!R_FINITE(*res))
    *res = LOW;

  return;
}

double pgev(double x, double loc, double scale, double shape)
{
  double y=0.0, result=0.0;

  y = (x - loc) / scale;

  if(shape==0)
    result = exp(-exp(-y));
  else
    result = exp(-pow(fmax(1 + shape * y, 0), - 1 / shape));

  return result;
}

void pGEV(double *x, double *loc, double *scale, double *shape, double *res)
{
    res[0] = pgev(x[0], loc[0], scale[0], shape[0]);
}


double qgev(double x, double loc, double scale, double shape)
{
  double res=0.0;

  if(shape==0)
    res = loc - scale * log(-log(x));
  else
    res = loc + scale * (pow(-log(x), -shape) - 1) / shape;

  return res;
}

void qGEV(double *x, double *loc, double *scale, double *shape, double *res)
{
    res[0] = qgev(x[0], loc[0], scale[0], shape[0]);
}
