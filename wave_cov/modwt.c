#include "mex.h"
#include "math.h"

/* mex -v -f /usr/tuelocal/matlab-5.2/bin/mexopts.sh modwt.c */

void modwt(double *Vin, int *N, int *j, int *L, double *ht, double *gt, 
          double *Wout, double *Vout)
{

  int k, n, t;

  for(t = 0; t < *N; t++) {
    k = t;
    Wout[t] = ht[0] * Vin[k];
    Vout[t] = gt[0] * Vin[k];
    for(n = 1; n < *L; n++) {
      k -= (int) pow(2.0, (double) *j - 1.0);
      if(k < 0) k += *N;
      Wout[t] += ht[n] * Vin[k];
      Vout[t] += gt[n] * Vin[k];
    }
  }

}

void mexFunction(int nlhs, mxArray *plhs[], 
                int nrhs, const mxArray *prhs[])
{
  int M, J, L;
  int m, n;
  double *Wout, *Vout;
  double *Vin, *ht, *gt;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 4) {
    mexErrMsgTxt("DWT requires four input arguments.");
  } else if (nlhs > 2) {
    mexErrMsgTxt("DWT requires two output arguments.");
  }
  
  Vin = mxGetPr(prhs[0]);
  ht = mxGetPr(prhs[1]);
  gt = mxGetPr(prhs[2]);
  J = (int) mxGetScalar(prhs[3]);
  /* mexPrintf("J = %d\n", J); */

  m = mxGetM(prhs[0]);
  /* mexPrintf("m = %d\n", m); */
  n = mxGetN(prhs[0]);
  /* mexPrintf("n = %d\n", n); */
  
  /* Create matrices for the return arguments */
  
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
  
  /* Assign pointers to the various parameters */
  
  Wout = mxGetPr(plhs[0]);
  Vout = mxGetPr(plhs[1]);
  
  M = mxGetNumberOfElements(prhs[0]);
  L = mxGetNumberOfElements(prhs[1]);
  
  /* Do the actual computations in a subroutine */
  
  modwt(Vin, &M, &J, &L, ht, gt, Wout, Vout);
  return;
}
