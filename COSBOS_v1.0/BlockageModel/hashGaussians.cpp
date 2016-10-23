/**
 * Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>,
 * Signal Analysis and Machine Perception Laboratory,
 * Department of Electrical, Computer, and Systems Engineering,
 * Rensselaer Polytechnic Institute, Troy, NY 12180, USA
 */

/** 
 * This is the C++/MEX code for hasing Gaussians and distances
 *
 * compile: 
 *     mex hashGaussians.cpp
 *
 * usage:
 *     H=hashGaussians(sensors,lights,dim,sigma)
 *         sensors: 3D spatial coordinates of sensors
 *         lights: 3D spatial coordinates of lights
 *         dim: 3D dimension of the room
 *         sigma: the standard deviation of the Gaussian kernel
 *         H: the data to be hashed
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/* return two numbers, d and on
 * d is distance
 * on is an indicator for whether projection is on line segment
 */
double *pointToLineDistance(double x, double y, double z, double x1, double y1, double z1, double x2, double y2, double z2)
{
    double *result=new double[2];
    
    double x3=x2-x1;
    double y3=y2-y1;
    double z3=z2-z1;
    
    double alpha = ( x3*(x-x1) + y3*(y-y1) + z3*(z-z1) ) / ( x3*x3 + y3*y3 + z3*z3 );
    
    double xv = x1 - x + alpha * x3;
    double yv = y1 - y + alpha * y3;
    double zv = z1 - z + alpha * z3;
    
    result[0]=sqrt(xv*xv+yv*yv+zv*zv);
    
    if(0<=alpha && alpha<=1)
    {
        result[1]=1;
    }
    else
    {
        result[1]=0;
    }
        
    return result;
}




void generateGaussian(double *H, double *sensors, double *lights, int *dim, double sigma, int ns, int nl)
{
    long i,j; // important, must be long
    long a; // a is a copy of i
    int x,y,z; // positions in V
    double x1,y1,z1; // sensor position
    double x2,y2,z2; // light position
    int sc,lc,s,l; // sensor/light channel, and sensor/light
    double *result; // the result of calling pointToLineDistance() 
    double d; // distance from point to line
    int on; // whether point projection is on line (1/0)
    double gaussian;

    cout<<"Rendering volume H ..."<<endl;
    for(i=0;i<dim[0]*dim[1]*dim[2];i++)
    {
        a=i;
        x=a%dim[0];
        a/=dim[0];
        y=a%dim[1];
        z=a/dim[1];
        
        for(j=0;j<ns*nl;j++)
        {
            s=j%ns;
            l=j/ns;
            
            x1=sensors[s];
            y1=sensors[s+ns];
            z1=sensors[s+ns*2];
            x2=lights[l];
            y2=lights[l+nl];
            z2=lights[l+nl*2];
            
            result=pointToLineDistance((double)x, (double)y, (double)z, x1, y1, z1, x2, y2, z2);
            d=result[0];
            on=(int)result[1];
            
            gaussian=exp(-d*d/2/sigma/sigma);
            H[i+dim[0]*dim[1]*dim[2]*j]=gaussian;

            // free result; 
            delete[] result;
        }
    }
}



/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *sensors;
    double *lights;
    int ns; // number of sensors, 12
    int nl; // number of lights, 12
    int *dim;
    double *dim2;
    double sigma;
    double *H;
    
    /*  check for proper number of arguments */
    if(nrhs!=4)
    {
        mexErrMsgIdAndTxt( "MATLAB:volumevolumeHashing:invalidNumInputs",
                "Four inputs required.");
    }
    if(nlhs>1)
    {
        mexErrMsgIdAndTxt( "MATLAB:volumeHashing:invalidNumOutputs",
                "One output required.");
    }
    
    /* check to make sure sigma is a scalar */
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
            mxGetN(prhs[3])*mxGetM(prhs[3])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:volumeHashing:sigmaNotScalar",
                "Input sigma must be a scalar.");
    }
    
    /*  get the scalar input sigma */
    sigma = (int) mxGetScalar(prhs[3]);
    

    /*  create a pointer to the input matrix sensors */
    sensors = mxGetPr(prhs[0]);
    
    /*  create a pointer to the input matrix lights */
    lights = mxGetPr(prhs[1]);
    
    /*  create a pointer to the input matrix dim */
    dim2 = mxGetPr(prhs[2]);
    dim = new int[3];
    dim[0]=(int)dim2[0];
    dim[1]=(int)dim2[1];
    dim[2]=(int)dim2[2];
    
    
    
    /*  get the number of sensors */
    ns = (int) mxGetM(prhs[0]);
    cout<<"Number of sensors: "<<ns<<endl;
    
    /*  get the number of lights */
    nl = (int) mxGetM(prhs[1]);
    cout<<"Number of lights: "<<nl<<endl;
            
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericMatrix((long)dim[0]*dim[1]*dim[2]*ns*nl, 1, 
         mxDOUBLE_CLASS, mxREAL);
    
    /*  create a C++ pointer to a copy of the output matrix */
    H = mxGetPr(plhs[0]);
    
    /*  call the C subroutine */
    generateGaussian(H,sensors,lights,dim,sigma,ns,nl);
    
    return;
    
}
