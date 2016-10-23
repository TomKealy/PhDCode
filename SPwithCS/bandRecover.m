function x = bandRecover(y,senseOpts,supp,Psi,modulator,normBound)
%% Given measurements y and supp(x), recover the signal x
% 
% Usage: x = bandRecover(y,senseOpts,supp,Psi,modulator,normBound)
% Inputs: y - signal to be projected
%         senseOpts - options specifying sensing architecture
%         supp - vector indicating occupied frequency bands
%         Psi - baseband DPSS basis used to capture each band
%         modulator - matrix with columns corresponding to tones for each 
%                     possible center frequency
%         normBound - an upper bound on the norm of the signal to be
%                     recovered for use in Tikhonov regularization
%
% Outputs: x - recovered signal
%
% Most recent change - 9/8/2011
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

%% Initilize parameters
[junk,S] = size(Psi);
M = length(y);
K = length(supp);

%% Form reduced orthogonal basis for identified space
basis = repmat(Psi,1,K).*kron(modulator(:,supp),ones(1,S));
basis = orth(basis);

%% Form Phi*Psi, i.e., the projected version of the basis
projBasis = zeros(M,size(basis,2));
for jj=1:size(basis,2),
    projBasis(:,jj) = sense(basis(:,jj),0,senseOpts); 
end

%% Recover alpha via Tikhonov regularized least-squares
[U,s,V] = csvd(projBasis);
[alpha,lambda] = lsqi(U,s,V,y,normBound);

%% From alpha, generate x
x = basis*alpha;