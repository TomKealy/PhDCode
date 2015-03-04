SPARCO: A Matlab Toolbox for Sparse Reconstruction
--------------------------------------------------

1. Introduction
===============

Thank you for downloading the SPARCO toolbox!  SPARCO is an
environment for testing and benchmarking algorithms for sparse
reconstruction.  It includes a collection of sparse reconstruction
problems that have appeared in the compressed sensing literature.  The
toolbox is also a framework for implementing new test problems, and
includes a library of linear operators that appears in this context.

At the core of the sparse recovery problem is the linear system

    A*x + r = b,

where A is an m-by-n linear operator and the m-vector b is the
observed signal.  The goal is to find a sparse n-vector x such that r
is small in norm.

Please visit the Sparco website to download the latest version, and to
see information on the test problems and operators:

  http://www.cs.ubc.ca/labs/scl/sparco

2. (Optional) external packages
===============================

SPARCO is prepackaged with all its required dependencies, but a few of
the test problems may rely on external packages.  These packages only
need to be installed if you wish to use these particular problems.

The curevelet-based test problems 50-51 rely on the CurveLab toolbox:

    CurveLab     http://www.curvelet.org/software.html

3. Installation and Setup
=========================

Start Matlab and make sure that the working directory is set to the
directory that contains the SPARCO source files.  At the Matlab
prompt, run

  >> sparcoSetup

This script adds various directories to your Matlab path.  The script
will try to permanently add these directories to your path (in
pathdefs.m), but may fail if that file is read-only.  In that case,
please copy and paste to your startup.m file the 'addpath' commands
printed to the screen.

To verify your installation, run

  >> checkProblems

This script instantiates each of the test problems and does some quick
testing to be sure all is in order.

4. Quick Start
==============

The main interface to the test problems is through the main 'sparco'
script found at the top-level directory.  To instantiate a particular
test problem, simply call the main 'sparco' script with the problems
ID number or name.  For example, to instantiate problem 5, do

   >> prob = generateProblem(5);

or

   >> prob = generateProblem('gcosspike');

This create a problem structure ('prob' in this case) which contains
all the information needed to access this problem.  This structure
contains many bits and pieces, including function handles and data
vectors.  The most important components are

   prob.A      a function handle to the operator A.
   prob.b      the right-hand-side vector.
   prob.sizeA  is a tuple with the number of rows and columns in A.

The function handle prob.A behaves as follows:

   z = prob.A( x, mode )  gives  z = A *x    if  mode == 1
                                 z = A'*x    if  mode == 2
                                 z = size(A) if  mode == 0.

It is also possible to use the classOp tool to instantiate the
operator as an overloaded object, e.g.,

   >> C = classOp(prob.A);  b = prob.b;
   >> g = C'*b;  # Equivalent to g = prob.A(b,2);
   >> y = C *g;  # Equivalent to y = prob.A(g,1);

The best way to get started is to browse the examples in the EXAMPLES
subdirectory.

To get a full list of problem numbers, do

   >> plist = generateProblem('list');

The conversion between problem name or number is done using for example

   >> name  = generateProblem(5,'getname');
   >> index = generateProblem('gcosspike','lookup');

5. The Test Problems
====================

An up-to-date list of test problems is maintained at

  http://www.cs.ubc.ca/labs/scl/sparco/

6. License
==========

SPARCO is open-source code released under the GNU Public License.  We
are delighted to acknowledge several other open-source packages that
we could build on:

- Rice Wavelet Toolbox:  http://www.dsp.rice.edu/software/rwt.shtml
- Nonuniform FFT Toolbox:   http://www.eecs.umich.edu/~fessler/code

Several individuals have given kind permission to include section of
code from their packages:

- Michael Lustig, SparseMRI, http://www.stanford.edu/~mlustig/SparseMRI.html
- Stephen J. Wright, GPSR, http://www.lx.it.pt/~mtf/GPSR/

7. The Authors
==============

We hope that SPARCO proves useful in your experiments.  If you have
any bug reports or comments, please feel free to email one of the
toolbox authors:

  Ewout van den Berg <ewout78@cs.ubc.ca>
  Michael P. Friedlander <mpf@cs.ubc.ca>

Enjoy!
Ewout and Michael
27 October 2007

$Id: README 787 2008-02-16 00:55:44Z ewout78 $
