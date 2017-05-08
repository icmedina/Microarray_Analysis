#include "mex.h"

/* mex -v -f /usr/tuelocal/matlab-5.2/bin/mexopts.sh dwt.c */

void dwt(double *Vin, int *M, int *L, double *h, double *g, 
         double *Wout, double *Vout)
{

  long n, t, u;

  for(t = 0; t < *M/2; t++) {
    u = 2 * t + 1;
    Wout[t] = h[0] * Vin[u];
    Vout[t] = g[0] * Vin[u];
    for(n = 1; n < *L; n++) {
      u -= 1;
      if(u < 0) u = *M - 1;
      Wout[t] += h[n] * Vin[u];
      Vout[t] += g[n] * Vin[u];
    } 
  }
}

void mexFunction(int nlhs, mxArray *plhs[], 
                int nrhs, const mxArray *prhs[])
{
  int M, L;
  int m, n;
  double *Wout, *Vout;
  double *Vin, *h, *g;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 3) {
    mexErrMsgTxt("DWT requires three input arguments.");
  } else if (nlhs > 2) {
    mexErrMsgTxt("DWT requires two output arguments.");
  }
  
  Vin = mxGetPr(prhs[0]);
  h = mxGetPr(prhs[1]);
  g = mxGetPr(prhs[2]);

  m = mxGetM(prhs[0]);
  mexPrintf("m = %d\n", m);
  n = mxGetN(prhs[0]);
  mexPrintf("n = %d\n", n);
  
  /* Create matrices for the return arguments */
  
  plhs[0] = mxCreateDoubleMatrix(m, n/2, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m, n/2, mxREAL);
  
  /* Assign pointers to the various parameters */
  
  Wout = mxGetPr(plhs[0]);
  Vout = mxGetPr(plhs[1]);
  
  M = mxGetNumberOfElements(prhs[0]);
  L = mxGetNumberOfElements(prhs[1]);
  
  /* Do the actual computations in a subroutine */
  
  dwt(Vin, &M, &L, h, g, Wout, Vout);
  return;
}

