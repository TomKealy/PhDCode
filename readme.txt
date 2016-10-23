DPSS Approximation and Recovery Toolbox (DART) version 1.0


Compressive sensing (CS) has recently emerged as a framework for efficiently capturing signals that are sparse or compressible in an appropriate basis. While often motivated as an alternative to  Nyquist-rate sampling, there remains a gap between the discrete, finite-dimensional CS framework and the problem of acquiring a continuous-time signal. 

This software package provides a set of tools for bridging this gap through the use of Discrete Prolate Spheroidal Sequences (DPSS's), a collection of functions that trace back to the seminal work by Slepian, Landau, and Pollack on the effects of time-limiting and bandlimiting operations. DPSS's form a highly efficient basis for sampled bandlimited functions; by modulating and merging DPSS bases, we obtain a dictionary that offers high-quality sparse approximations for most sampled multiband signals. This multiband modulated DPSS dictionary can be readily incorporated into the CS framework. 

For further details, see the paper "Compressive Sensing of Analog Signals Using Discrete Prolate Spheroidal Sequences," by Mark A. Davenport and Michael B. Wakin (included in the distribution.) 

DART contains all of the software necessary to reproduce the results presented in this paper (see the "experiments" directory for the specific results presented).

To begin experimenting with DPSS's for CS, first look at the demo.m file, which shows the core routines of DART in action.  The "tests" directory provides a number of files which illustrate the performance of each of the core components of DART.

Included in the distribution are the files csvd.m and lsqi.m, which are part of the regtools package. See the included license file for further details.

Please e-mail markad@stanford.edu if you find any bugs or have any questions. 