function op = opHaar(n,levels)
%OPHAAR   One-dimensional Haar wavelet
%
%   OPHAAR(N,LEVELS) creates a Haar wavelet operator for one
%   dimensional signal of length N. LEVELS indicates the number of
%   scales to use, default is 5.
%
%   See also opWavelet

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opHaar.m 1040 2008-06-26 20:29:02Z ewout78 $

if (nargin < 2), levels = 5; end

op = opWavelet(n,1,'haar',0, levels,'min');
