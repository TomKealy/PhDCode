#include <mex.h>

void mexFunction(int nlhs, mxArray * plhs[],int nrhs, const mxArray * prhs[]) 
{
   if ( nrhs != 2 )
   {
       mexErrMsgTxt("Two arguments needed.");
       return;
   }
   if ( nlhs != 1 )
   {
       mexErrMsgTxt("One output needed.");
       return;
   }   
   const mwSize * dimsX  = mxGetDimensions(prhs[0]);
   const mwSize * dimsM  = mxGetDimensions(prhs[1]);
   if ( dimsX[0] != dimsM[0])
   {
       mexErrMsgTxt("The column and the matrix must have same height.");
       return;
   }
   if ( mxIsLogical(prhs[0]) && mxIsLogical(prhs[1]) )
   {
       bool * C   = (bool*)mxGetPr(prhs[0]);
       bool * M   = (bool*)mxGetPr(prhs[1]);
       plhs[0] = mxCreateNumericArray(2,dimsM,mxLOGICAL_CLASS,mxREAL);
       bool * Y   = (bool*)mxGetPr(plhs[0]);
       for(int col = 0 ; col < dimsM[1] ; col++)
           for(int r = 0 ; r < dimsM[0] ; r++)
               *Y++ = *M++ && C[r];
       return;
   }
   mexErrMsgTxt("Arguments must both be 'logical'.");
}
