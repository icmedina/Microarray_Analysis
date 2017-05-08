#include "mex.h"

/* mex -v -f /usr/tuelocal/matlab-5.2/bin/mexopts.sh dwt.c */

void idwt(double *Win, double *Vin, int *M, int *L, double *h, double *g, 
         double *Xout)
{

  int i, j, l, t, u;
  int m = -2, n = -1;

  for(t = 0; t < *M; t++) {
    m += 2;
    n += 2;
    u = t;
    i = 1;
    j = 0;
    Xout[m] = h[i] * Win[u] + g[i] * Vin[u];
    Xout[n] = h[j] * Win[u] + g[j] * Vin[u];
    if(*L > 2) {
      for(l = 1; l < *L/2; l++) {
       u += 1;
       if(u > *M) u = 0;
       i += 2;
       j += 2;
       Xout[m] += h[i] * Win[u] + g[i] * Vin[u];
       Xout[n] += h[j] * Win[u] + g[j] * Vin[u];
      }
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], 
                int nrhs, const mxArray *prhs[])
{
  int M, L;
  int m, n;
  double *Xout;
  double *Win, *Vin, *h, *g;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 4) {
    mexErrMsgTxt("DWT requires four input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("DWT requires one output argument.");
  }
  
  Win = mxGetPr(prhs[0]);
  Vin = mxGetPr(prhs[1]);
  h = mxGetPr(prhs[2]);
  g = mxGetPr(prhs[3]);

  m = mxGetM(prhs[0]);
  mexPrintf("m = %d\n", m);
  n = mxGetN(prhs[0]);
  mexPrintf("n = %d\n", n);
  
  /* Create matrices for the return arguments */
  
  plhs[0] = mxCreateDoubleMatrix(m, 2*n, mxREAL);
  
  /* Assign pointers to the various parameters */
  
  Xout = mxGetPr(plhs[0]);
  
  M = mxGetNumberOfElements(prhs[0]);
  L = mxGetNumberOfElements(prhs[2]);
  
  /* Do the actual computations in a subroutine */
  
  idwt(Win, Vin, &M, &L, h, g, Xout);
  return;
}

