/**
 * Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>,
 * Signal Analysis and Machine Perception Laboratory,
 * Department of Electrical, Computer, and Systems Engineering,
 * Rensselaer Polytechnic Institute, Troy, NY 12180, USA
 */

/** 
 * This is the C++/MEX code for rendering the volume from hashed Gaussians
 *
 * compile: 
 *     mex volumeFromHashing.cpp
 *
 * usage:
 *     V=volumeFromHashing(sensors,lights,dim,H,E)
 *         sensors: 3D spatial coordinates of sensors
 *         lights: 3D spatial coordinates of lights
 *         dim: 3D dimension of the room
 *         H: the data that has been hashed
 *         E: the difference matrix, E=A0-A 
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;



void generateVolume(double *V, int *dim, double *H, double *E, int ns, int nl)
{
    double **L=new double *[ns]; // aggregation of E
    long i,j; // important, must be long
    int sc,lc,s,l; // sensor/light channel, and sensor/light
    double normalizor,gaussian;
    long dimProd=(long)dim[0]*dim[1]*dim[2];
    
    for(i=0;i<ns;i++)
    {
        L[i]=new double[nl];
        for(j=0;j<nl;j++)
        {
            L[i][j]=0;
        }
    }
    
    cout<<"Constructing matrix L ..."<<endl;
    for(i=0;i<4*ns*3*nl;i++)
    {
        sc=i%(4*ns);
        lc=(i-sc)/(4*ns);
        s=sc/4;
        l=lc/3;
        if(sc%4 == lc%3)
        {
            L[s][l]=L[s][l]+E[i];
        }
    }

    cout<<"Rendering volume V ..."<<endl;
    for(i=0;i<dimProd;i++)
    {
        normalizor=0;
        V[i]=0;
        for(j=0;j<ns*nl;j++)
        {
            s=j%ns;
            l=j/ns;
            
            gaussian=H[i+dimProd*j];
            V[i]+=L[s][l]*gaussian;
            normalizor+=gaussian;
        }
        V[i]/=normalizor;
    }
    
    // free
    for(i=0;i<11;i++)
    {
        delete[] L[i];
    }
    delete[] L;
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
    double *H;
    double *E;
    double *V;
    
    /*  check for proper number of arguments */
    if(nrhs!=5)
    {
        mexErrMsgIdAndTxt( "MATLAB:volumeFromHashing:invalidNumInputs",
                "Five inputs required.");
    }
    if(nlhs>1)
    {
        mexErrMsgIdAndTxt( "MATLAB:volumeFromHashing:invalidNumOutputs",
                "One output required.");
    }

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
    
    /*  create a pointer to the input matrix A */
    H = mxGetPr(prhs[3]);
    E = mxGetPr(prhs[4]);
    
    /*  get the number of sensors */
    ns = (int) mxGetM(prhs[0]);
    cout<<"Number of sensors: "<<ns<<endl;
    
    /*  get the number of lights */
    nl = (int) mxGetM(prhs[1]);
    cout<<"Number of lights: "<<nl<<endl;
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericArray(3, dim, 
         mxDOUBLE_CLASS, mxREAL);
    
    /*  create a C++ pointer to a copy of the output matrix */
    V = mxGetPr(plhs[0]);
    
    /*  call the C subroutine */
    generateVolume(V,dim,H,E,ns,nl);
    
    return;
    
}
