% Sparco Toolbox.
% Testing environment for sparse reconstruction.
% Version 1.2  25-June-2008
% by   Ewout van den Berg (ewout78@cs.ubc.ca)
% and  Michael P. Friedlander (mpf@cs.ubc.ca)
% Copyright 2008, University of British Columbia
%
% SPARCO
%   README           - About the SPARCO toolbox.
%   sparco           - Interface to individual problems.
%   sparcoSetup      - Setup the SPARCO toolbox.
%   literature.bib   - BibTeX file with cites to problem sources.
%
% EXAMPLES
%   ex1              - General demo.
%   exGPSR           - Example that calls GPSR (GPSR not included).
%   exOp             - Operator demo.
%   exSPGL1          - Example that calls SPGL1 (SPGL1 not included).
%   exTV             - Total variation demo (requires SparseMRI).
%   
% OPERATORS
%   opBinary         - Binary (0/1) ensemble.
%   opBlockDiag      - Diagonal operator.
%   opBlur           - Two-dimensional blurring operator.
%   opColumnRestrict - Restriction operator on matrix columns.
%   opConvolve1d     - One-dimensional convolution operator.
%   opCurvelet2d     - Two-dimensional curvelet operator.
%   opDCT            - One-dimensional discrete cosine transform (DCT).
%   opDiag           - Scaling operator (i.e., diagonal matrix).
%   opDictionary     - Dictionary of concatenated operators.
%   opDirac          - Identity operator.
%   opFFT            - One-dimensional fast Fourier transform (FFT).
%   opFFT2C          - Centralized two-dimensional fast Fourier transform (FFT).
%   opFFT2d          - Two-dimensional fast Fourier transform (FFT).
%   opFoG            - Concatenate sequence of operators into a single operator.
%   opGaussian       - Gaussian ensemble.
%   opHaar           - One-dimensional Haar wavelet.
%   opHaar2d         - Two-dimensional Haar wavelet.
%   opHadamard       - Hadamard matrix.
%   opHeaviside      - Heaviside operator.
%   opMask           - Selection mask.
%   opMatrix         - Apply arbitrary matrix as operator.
%   opPadding        - Padding operator.
%   opReal           - Discard imaginary components.
%   opRestriction    - Selection and restriction operator.
%   opSign           - Sign-ensemble operator.
%   opSplitComplex   - Split a complex vector into real and imaginary parts.
%   opTranspose      - Transpose operator.
%
% PROBLEMS
%   prob001          - FFT/Heaviside dictionary, HeaviSine singal
%   prob002
%   ...
% 
% TOOLS
%

% $Id: Contents.m 1040 2008-06-26 20:29:02Z ewout78 $
