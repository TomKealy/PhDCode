function [s,idxNZ] = sparseComplexSigGen4plusNoise(N,NZ,sigma_off,fixedActiveValue)
% - Generates an Nx1 random sparse vector with "N-NZ" number of its components 
%   (randomly) chosen to be nearly zero (inactive components). That is NZ is
%   the number of Non-Zero (active) components.
% - The rest of the compoents will have the fixed value "fixedActiveValue" 
%   active components). If not supplied with the "fixedActiveValue", the
%   function generates a gaussian noise with unit varinace over the active
%   components.
% - The "sigma_off" parameter controls the variance of (gaussian) noise over the
%   inactive components.
% - The function also returns the indices of the Non-zero (i.e., active)
%   components.
if nargin < 2
    NZ = 0.1*N;
end

if nargin < 3
    sigma_off = 0.01;
end

idxNZ = zeros(N,1);
rndPerm = randperm(N);
idxNZ(rndPerm(1:NZ)) = 1;
idxNZ = logical(idxNZ);

s = zeros(N,1);
if nargin < 4
    s(idxNZ) = 1/sqrt(2)*randn(NZ,1) + sqrt(-1)/sqrt(2)*randn(NZ,1);
else
    s(idxNZ) = fixedActiveValue;
end
s(~idxNZ) = sigma_off/sqrt(2)*randn(N-NZ,1) + sigma_off*sqrt(-1)/sqrt(2)*randn(N-NZ,1);

% s = filter(1,[1,-0.2],s);