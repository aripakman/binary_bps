/*==========================================================
 * Author: ari pakman
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include <iostream>
#include <math.h>
#include <stdint.h>
#include "mex.h"

#include "BouncyBinarySampler.h"
#include "BinaryDistribution.h"
#include "MRF.h"

using namespace std;

/****************************/

typedef pair< vector<int> ,double> ST_Pair;    

static
void bouncy_binary(

  mxArray * plhs[],
  int Utype,     //0: exponential, else: gaussian
  double * M,
  double * r,
  double max_time,
  int d, 
  double * last_Y,
  int * seed_ptr ){


  vector < ST_Pair >  STs;
  vector<double> log_likes; 

  vector <double>  Y(d);
  for (int i=0; i< d; i++)
    Y[i] = last_Y[i];


  
  MRF * mrf  = new MRF(d, M, r);

  BouncyBinarySampler sampler = BouncyBinarySampler(mrf);

  if (Utype == 0) {
    sampler.runSamplerExponential(max_time, STs, log_likes, seed_ptr, Y );
  }     
  else{
    sampler.runSamplerGaussian(max_time, STs, log_likes, seed_ptr, Y );
  }


  for (int i=0; i< d; i++)
    last_Y[i] = Y[i];


  int L = STs.size();

  // Create variables to return
  plhs[0] = mxCreateNumericMatrix(d, L, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, L, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericMatrix(1, L, mxDOUBLE_CLASS, mxREAL);


  double * samples  = mxGetPr(plhs[0]);
  double * times = mxGetPr(plhs[1]);
  double * loglikes = mxGetPr(plhs[2]);

  //copy the samples, times and likes to the variables to be returned
  for(int i=0; i < L; i++) {

    loglikes[i] = log_likes[i];
    times[i] = STs[i].second;
    
    for (int j=0; j < d; j++) {
      
      samples[d*i+j ] =  STs[i].first[j];      
    }
  }


  return;
}


void mexFunction(
		 int          nlhs,
		 mxArray     *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
 
  
    // Input variables
    int Utype;
    double *M;                /* NxN input matrix */
    double *r;                /* 1xN input matrix */
    double max_time;
    double *last_Y;
    int * seed_ptr = new int();

  /* Check for proper number of arguments */

  if (nrhs != 5 && nrhs != 6) {
    cout << nrhs << endl;
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "MEXCPP requires four or five input arguments.");
  } else if (nlhs > 3) {
    mexErrMsgIdAndTxt("MATLAB:bouncy_binary:nargout",
            "hbouncy_binary yields up to three output arguments.");
  }

  Utype = (int) mxGetScalar(prhs[0]);
  M = (double *) mxGetPr(prhs[1]);
  r = (double *) mxGetPr(prhs[2]);
  max_time = (double ) mxGetScalar(prhs[3]);
  last_Y = (double *) mxGetPr(prhs[4]);

  if (nrhs == 6)
     *seed_ptr =  (int ) mxGetScalar(prhs[5]);
   else 
      *seed_ptr = rand();


  /* check that number of cols in r is 1 */


  if(mxGetN(prhs[2])!=1) {
     mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input r must be a column vector.");
  }

  int d = mxGetM(prhs[1]);

  if(mxGetM(prhs[1])!=d || mxGetN(prhs[1])!=d) {
     mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input M must be a square matrix with same dim as r.");
  }


  bouncy_binary(plhs, Utype, M,r, max_time, d, last_Y, seed_ptr);

  delete seed_ptr;

  return;
}
