#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include <string>
#include <cmath>

//   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.


  //Equivalent matlab code:
  // D  = diag(sqrt(Zeta))
  // X_ = X * D;
  // w  = D * ( ( X_'*X_ + n*lambda*I) \ X_'*Y );

  // Argument 1: X,     whose dim is n x p
  //          2: Y,     whose dim is n x 1
  //          3: Zeta,  whose dim is p x 1
  //          4: lambda

/* Prototype for routines */

//extern void dcopy_(int*, double*, int*, double*, int*);
//extern void dscal_(int*, double*, double*, int*);
//extern void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*, double*, int*);
//extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
//extern void dposv_(char*, int*, int*, double*, int*, double*, int*, int*);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *X, *Y, *Zeta, *w;
  double one_db = 1.0, zero_db = 0.0;

  int one_int = 1;
  int n  = mxGetM(prhs[0]);
  int p  = mxGetN(prhs[0]);
  int np = n * p;
  int info;
          
  const double lambda_n =  (*mxGetPr(prhs[3])) * n;

  std::string NoTranspose = "N", Transpose = "T", Lower = "L";

  X    = (double*)mxGetPr(prhs[0]);
  Y    = (double*)mxGetPr(prhs[1]);
  Zeta = (double*)mxGetPr(prhs[2]);

  //if (p != mxGetM(prhs[1])) {
  //  mexErrMsgTxt("Inner dimensions of matrix multiply do not match");
  //}

  //--------------------------------------------------------------------------------
  //X_ is initialized as a copy of X

  double *X_ = new double[n*p];
  
  dcopy_( &np, X, &one_int, X_, &one_int );

 
  //--------------------------------------------------------------------------------
  // scale the columns of X_ with sqrt_zeta

  double *sqrt_Zeta = new double[p];

  for ( int j = 0; j < p; j++ )
    {
      sqrt_Zeta[j] = sqrt( Zeta[j] );

      dscal_ ( &n, sqrt_Zeta+j, X_+j*n, &one_int );
    }

  //--------------------------------------------------------------------------------

  
  // We fill only the lower part of X_tX_plus_lambda_Id
  
  double *X_tX_plus_lambda_Id = new double[p*p];

  for (int col = 0; col < p; col++ )
  {  
      
    X_tX_plus_lambda_Id[col*p + col] = lambda_n;
    
	for ( int row = col+1; row < p; row++ )
        
        X_tX_plus_lambda_Id[col*p + row] = 0.0;
    
  }
  
  
  // X_tX_plus_lambda_Id <- X_^T * X_ + X_tX_plus_lambda_Id

  dsyrk_ ( &Lower[0], &Transpose[0], &p, &n, &one_db, X_, &n, &one_db,  X_tX_plus_lambda_Id, &p );
  
  //--------------------------------------------------------------------------------

  
  plhs[0]  = mxCreateDoubleMatrix(p, 1, mxREAL);
  w        = (double*)mxGetPr( plhs[0] );
  
  // w <- X_' * Y + 0*w

  dgemv_ (&Transpose[0], &n, &p, &one_db, X_, &n, Y, &one_int, &zero_db, w, &one_int );
  
  //--------------------------------------------------------------------------------
  
  // w <- X_tX_plus_lambda_Id^(-1) X_' * Y

  dposv_ ( &Lower[0], &p, &one_int, X_tX_plus_lambda_Id, &p, w, &p, &info );
  
  //--------------------------------------------------------------------------------

  if ( info != 0 )
	mexErrMsgTxt( "\n Error while computing the matrix inversion ! \n " );


  for (int j = 0; j < p; j++)
	w[j] *= sqrt_Zeta[j];

  
  delete [] X_;
  delete [] X_tX_plus_lambda_Id;
  delete [] sqrt_Zeta;

}
