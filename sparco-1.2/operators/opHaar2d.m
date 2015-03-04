function op = opHaar2d(m,n,levels)
%OPHAAR2D   Two-dimensional Haar wavelet
%
%   OPHAAR2D(N,LEVELS) creates a Haar wavelet operator for two
%   dimensional signal of size M by N. LEVELS indicates the number of
%   scales to use, default is 5.
%
%   See also opWavelet

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opHaar2d.m 1040 2008-06-26 20:29:02Z ewout78 $

if (nargin < 3), levels = 5; end

op = opWavelet(m,n,'haar',0, levels,'min');
