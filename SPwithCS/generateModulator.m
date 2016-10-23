function modulator = generateModulator(N,B)
%% Generates a modulator matrix
% 
% Usage: modulator = generateModulator(N,B)
% Inputs: N - signal length
%         B - bandwidth of each band 
%             (B must be <= .5 and 1/B must be an integer)
% Outputs: modulator - matrix of modulation tones
%
% Most recent change - 8/27/2011
%
% Copyright 2011, Mark Davenport, Michael Wakin
%
% This file is part of DPSS Approximation and Recovery Toolbox version 1.0.
%
%    DART is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DART is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DART.  If not, see <http://www.gnu.org/licenses/>.

%% Construct a vector of frequencies indicating center frequencies 
% for bands of width B with edges at -.5, -.5+B, ... , .5-B, .5.  
% This translates to center frequencies of -.5+B/2, ... , .5-B/2    
f_vec = [-.5+B/2:B:.5-B/2];

%% Construct time index
t_vec = [0:1:N-1];

%% Generate modulator matrix
modulator = exp(2j*pi*t_vec'*f_vec);
