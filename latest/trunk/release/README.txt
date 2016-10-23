Thank you for downloading "Generalized Approximate Message
Passing (GAMP) in MRI" MATLAB package. If you have any questions or
difficulties, please contact Sundeep Rangan at srangan@poly.edu.

Installation
------------

1. Download the GAMP_MRI package from the sourceForge web page for GAMP,
   "https://gampmatlab.svn.sourceforge.net/svnroot/gampmatlab"

2. Download the Wavelab package from  (http://www-stat.stanford.edu/~wavelab/) to run the wavelet transforms.

    2.1. Download the Wavelab850.zip file.

    2.2. Extract the .zip file into the folder \matlab\toolbox\
         Note that if your Matlab root directory is named differently than \matlab then use its correct name instead of \matlab.

    2.3. In Matlab set your current path to matlab\toolbox\WaveLab850
         Note that for continuing you need to have a C compiler.

    2.4. Run WavePath

    2.5. If you have a 64-bit machine download the modified version of InstallMEX from
         "https://gampmatlab.svn.sourceforge.net/svnroot/gampmatlab/trunk/release/" and replace it with the old one.

    2.6. Run InstallMEX
         Note that if you don't have permission to create files in the MATLAB toolbox directory,
         you can built the mex files in another directory. Then copy the entire directory to the destination.
         Remark: if you type "which FWT_PBS" command you should get FWT_PBS.mexw64.

3. The MATLAB files are in a subdirectory, you should add the directory to your MATLAB path.
   You can use the MATLAB addpath command.

4. If you enter code/test and run estimTest. This will run the GAMP algorithm for a Bernoulli-Gaussian sparse vector
   and compare the performance against a linear least-squares estimator.
   Also, you can enter code/MRI and run test program. There you can also find examples of
   Professor Otazo's compressed sensing reconstruction of Cartesian MRI data acquired with multiple coils code.
   if you want to run it, the data files are available upon your request.


Reading further
----------------
For more information, consult the User's Guide at:

http://gampmatlab.sourceforge.net/wiki/index.php/Users%27_Guide

The User's Guide has very little right now. There is a
reference to one example, but I am hoping to find more.

Also, due to a server error we are trying to debug, the webpage
may time several minutes to open the first time. So, please be
patient...