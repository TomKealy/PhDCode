% SparseLab Demo Utilities Directory, Version 100
%
% This is the utilities directory of SparseLab. Some .m files are
% part of WaveLab ver 802, http://www-stat.stanford.edu/~wavelab
%
%              .m files in this directory
%
%   Contents.m           -   This file
%   aconv.m              -   Convolution Tool for Two-Scale Transform
%   AutoImage.m          -   Automatic Scaling for Image Display
%   DownDyadHi.m         -   Hi-Pass Downsampling operator (periodized)
%   DownDyadLo.m         -   Lo-Pass Downsampling operator (periodized)
%   dyad.m               -   Index entire j-th dyad of 1-d wavelet xform
%   dyadlength.m         -   Find length and dyadic length of array
%   LockAxes.m           -   Version-independent axis command
%   lshift.m             -   Circular left shift of 1-d signal
%   MakeONFilter.m       -   Generate Orthonormal QMF Filter for Wavelet Transform
%   MirrorFilt.m         -   Apply (-1)^t modulation
%   Noisemaker.m         -   Add Noise to Signal
%   NormNoise.m          -   Estimates noise level, Normalize signal to noise level 1
%   packet.m             -   Packet table indexing
%   PlotSpikes.m         -   Plot 1-d signal as baseline with series of spikes
%   PlotWaveCoeff.m      -   Spike-plot display of wavelet coefficients
%   RegisterPlot.m       -   Add legend with file name, date, flag
%   reverse.m            -   Reverse order of elements in 1-d signal
%   rshift.m             -   Circular right shift of 1-d signal
%   ShapeAsRow.m         -   Reshape 1d vector as row
%   ShapeLike.m          -   Reshape first argument like second argument
%   TIDenoise.m          -   Translation invariant denoising of a 1-D signal
%   twonorm.m            -   Computes ||v||_2
%   UnlockAxes.m         -   Version-independent axis command
%   UpDyadHi.m           -   Hi-Pass Upsampling operator; periodized
%   UpDyadLo.m           -   Lo-Pass Upsampling operator; periodized
%   UpSample.m           -   Upsampling operator
%
%              Subdirectories
%
%   \BuildDatasets
%   \Transforms

    
%
% Copyright (c) 2007. Victoria Stodden
%  

%
% Part of SparseLab Version:200
% Created Tuesday March 24, 2007
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%