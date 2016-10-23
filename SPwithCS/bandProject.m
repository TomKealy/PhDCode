function x_hat = bandProject(x,supp,Psi,modulator)
%% Projects a signal x onto subspace indicated by supp
% 
% Usage: x_hat = bandProject(x,supp,Psi,modulator)
% Inputs: x - signal to be projected
%         supp - vector indicating occupied frequency bands
%         Psi - baseband DPSS basis used to capture each band
%         modulator - matrix with columns corresponding to tones for each 
%                     possible center frequency 
% Outputs: x_hat - projected signal
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

%% Initialize parameters
M = size(Psi,2);
K = length(supp);

%% Form complete basis for identified space 
basis = repmat(Psi,1,K).*kron(modulator(:,supp),ones(1,M));

%% Recover optimal x_hat
basis = orth(basis);
x_hat = basis*(ctranspose(basis)*x);
