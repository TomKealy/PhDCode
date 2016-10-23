function [supp varargout] = bandDetect(x,K,PsiDetect,PsiRecover,modulator,opts)
%% Find the frequency bands containing the most energy in x
%  Also computes the projection onto the set of K-band-sparse signals
% 
% Usage: [supp x_hat] = bandDetect(x,K,PsiDetect,PsiRecover,modulator,opts)
% Inputs: x - signal to be projected
%         K - signal number of active bands
%         PsiDetect - baseband DPSS basis for detection
%         PsiRecover - baseband DPSS basis for recovery
%         modulator - matrix with columns corresponding to tones for each possible center frequency 
%         opts - argument setting algorithm parameters 
%           opts.alg - selects detection algorithm
%                      current options: 'threshold', 'blockOMP'
% Outputs: supp - frequency band indices
%          x_hat - projected signal (optional)
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


%% Transpose DPSS bases and threshold
if strcmp(opts.alg,'threshold')
    % Modulate x to each possible frequency
    x_mod = repmat(x,1,size(modulator,2)).*conj(modulator);

    % Project each modulated vector into DPSS basis, which is done simply using
    % the matrix transpose since Psi is an ONB
    r = ctranspose(PsiDetect)*x_mod;

    % Find the K columns of r that have the most energy
    energy = sum(r.*conj(r),1);
    [junk,index] = sort(energy,'descend');
    supp = index(1:K);

    % Optionally compute the projection of x onto the identified bands
    if nargout==2
        varargout{1} = bandProject(x,supp,PsiRecover,modulator);
    end
end

%% Block OMP
if strcmp(opts.alg,'blockOMP')
    
    supp = [];
    x_hat = zeros(size(x));    
    for jj=1:K,
        residual = x-x_hat;
        % Modulate residual to each possible frequency
        x_mod = repmat(residual,1,size(modulator,2)).*conj(modulator);
        
        % Project each modulated vector into DPSS basis, which is done simply using
        % the matrix transpose since Psi is an ONB
        r = ctranspose(PsiDetect)*x_mod;
        
        % Find the column of r (not yet selected) that has the most energy
        energy = sum(r.*conj(r),1);
        energy(supp) = 0;
        [junk,index] = sort(energy,'descend');
        supp = union(supp,index(1));
        
        % Project x onto identified subspace to form x_hat
        x_hat = bandProject(x,supp,PsiRecover,modulator);
    end
    
    % Optionally output the projection of x onto the identified bands
    if nargout==2
        varargout{1} = x_hat;
    end
end