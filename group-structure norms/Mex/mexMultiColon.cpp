#include <mex.h>
#include <iostream>

//   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

void mexFunction(int nlhs, mxArray * plhs[],int nrhs, const mxArray * prhs[]) 
{

// Given I and J with |I|=|J|, the function outputs [I(1):J(1), I(2):J(2), ...]

	if ( nrhs < 2 )
	{
       		mexErrMsgTxt("Two arguments needed.");
       		return;
   	}

   	if ( nlhs != 1 )
   	{
       		mexErrMsgTxt("One output needed.");
       		return;
   	}

	const mwSize * dimsI  = mxGetDimensions(prhs[0]);
	const mwSize * dimsJ  = mxGetDimensions(prhs[1]);
	
   	if ( dimsI[0] != dimsJ[0] )
   	{
       		mexErrMsgTxt("Vectors must have the same length.");
       		return;
   	}	
	
	double * pI  = (double *)mxGetPr(prhs[0]);
	double * pJ  = (double *)mxGetPr(prhs[1]);
	
	int OutputSize = 0;
	
	for (int j = 0; j < dimsI[0]; j++) 
	  if ( pJ[j] >= pI[j] )
        OutputSize = OutputSize + int( pJ[j] - pI[j] + 1 );

	
	plhs[0] = mxCreateDoubleMatrix( OutputSize, 1, mxREAL);

	double * Output = (double *)mxGetPr(plhs[0]);
	
	for (int j = 0; j < dimsI[0]; j++)
	  {
	    for (int i = int(pI[j]); i <= int(pJ[j]); i++)
	      {
		*Output = i;
		Output++;
	      }
	  }
	  
	return;
}
