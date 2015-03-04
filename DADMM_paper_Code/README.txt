----------------------------------------------------------------------------------
D-ADMM: A Communication-Efficient Distributed Algorithm For Separable Optimization
Copyright (C) 2012-2013 João Mota, João Xavier, Pedro Aguiar, Markus Püschel.
----------------------------------------------------------------------------------

Thank you for downloading this software. 

The purpose of this software package is to make available the source code of the
experiments whose results are presented in the following paper.

J. Mota, J. Xavier, P. Aguiar, M. Püschel,
"D-ADMM: A Communication-Efficient Distributed Algorithm For Separable Optimization"
IEEE Transactions on Signal Processing, 2013

The implementation of the proposed algorithm, D-ADMM, is in the folder /D-ADMM. 
In the same folder, there are several application examples on how to use the code.

All code is implemented in Matlab. Although all algorithms are distributed, all 
code here is centralized and runs sequentially. The purpose is to simulate
distributed algorithms.


-------------------------------------------------------------------------------
This software is distributed under the GNU General Public License; see gpl.txt

This package is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
-------------------------------------------------------------------------------



Description of each folder:


  /CompareAlgs: contains the code to run the experiments presented in Fig.3 of
                the paper mentioned above. In each directory, there is a file 
                named RunExperimentsNodes.m. To replicate the plots of Fig.3,
                run RunExperimentsNodes(50). WARNING: these experiments take
                a long time to execute because each algorithm is run for several
                values of rho. The results of the experiments that were used in
                the paper are stored in the folder /Results.


  /D-ADMM: contains the code of the algorithm proposed in the paper mentioned 
           above. There are four other folders that contain a script to run 
           D-ADMM for a particular application, solvers for that application,
           and functions that measure the relative error.


  /D-Lasso: contains the code for the algorithm in 
           
            H. Zhu, G. Giannakis, and A. Cano, "Distributed In-Network Channel
            Decoding," IEEE Trans. Sig. Proc., Vol. 57, No. 10, 2009

            The name D-Lasso, which we decided to use here, is due to 

            J. Bazerque and G. Giannakis, "Distributed Spectrum Sensing for 
            Cognitive Radio Networks by exploiting Sparsity," IEEE Trans. Sig. 
            Proc., Vol. 58, No. 3, 2010

            This is our implementation of the algorithm.


  /D-MLE: contains the code for the algorithm in

          I. Schizas, A. Ribeiro, and G. Giannakis, "Consensus in Ad Hoc WSNs
          With Noisy Links - Part I: Distributed Estimation of Deterministic 
          Signals," IEEE Trans. Sig. Proc., Vol. 56, No. 1, 2008


          This is our implementation of the algorithm.


  /GenerateData: contains the all the data. This includes networks and problem
                 data; each one has its own folder.


  /SpecificProblemSolvers: contains the implementation of algorithms solve 
                           a particular case of a separable problem, i.e.,
                           they cannot be applied to all problem instances.
                           These algorithms are from the papers

                           B. Oreshkin, M. Coates, M. Rabbat, "Optimization and 
                           Analysis of Distributed Averaging With Short Node 
                           Memory," IEEE Trans. Sig. Proc., Vol. 58, No. 5, 2010

                           G. Mateos, J. Bazerque, G. Giannakis, "Distributed 
                           Sparse Linear Linear Regression," IEEE Trans. Sig. Proc., 
                           Vol. 58, No. 10, 2010








