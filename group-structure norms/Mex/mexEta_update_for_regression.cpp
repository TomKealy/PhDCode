#include "mex.h"
#include <cmath>

  // Argument 1: w,    whose dim is p x 1
  //          2: G,    whose dim is p x ng
  //          3: D,    whose dim is p x ng
  //          4: epsilon
  //   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  bool   *G;
  double *w, *D, *Eta, *Zeta;

  const double epsilon = *mxGetPr(prhs[3]);

  int p = mxGetM(prhs[0]), ng = mxGetN(prhs[1]);

  w       = (double*)mxGetPr(prhs[0]);
  G       = (bool*)  mxGetPr(prhs[1]);
  D       = (double*)mxGetPr(prhs[2]);

  plhs[0]  = mxCreateDoubleMatrix(p, 1, mxREAL);
  Zeta     = (double*)mxGetPr( plhs[0] );

  plhs[1]  = mxCreateDoubleMatrix(ng, 1, mxREAL);
  Eta      = (double*)mxGetPr( plhs[1] );

  for (int j = 0; j < p; j++)
    Zeta[j] = 0.0;	

  for (int g = 0; g < ng; g++)
    Eta[g] = 0.0;

  //--------------------------------------------------------------------------------

  double *w_sq          = new double[p];
  double *D_sq          = new double[p*ng];

  for (int j = 0; j < p; j++)
  	w_sq[j] = pow( w[j],  2 );  

  for (int i = 0; i < p*ng; i ++)
	if ( G[i] )
		D_sq[i] = pow( D[i], 2 );
    else
      	D_sq[i] = 0.0; 

  //--------------------------------------------------------------------------------
  
  for (int g = 0; g < ng; g++)
      {

	for (int row = 0; row < p; row++)
	{
	  if ( G[p*g+row] )
	    Eta[g] += w_sq[row] * D_sq[p*g+row];
	}


	Eta[g] = pow( Eta[g], 0.5 );

       }

  //--------------------------------------------------------------------------------
  double inv_eta;

  for (int g = 0; g < ng; g++)
	{

	inv_eta = 1.0 / ( Eta[g] + epsilon );

	for (int row = 0; row < p; row++)
  		if ( G[p*g+row] )
    			Zeta[row] += D_sq[p*g+row] * inv_eta;
	}

  for (int j = 0; j < p; j++)
      Zeta[j] = 1.0 / Zeta[j];


  delete [] w_sq;
  delete [] D_sq; 
}
