function A = MakeMatrix(n, p, type)
% MakeMatrix: creates various types of test matrices
% Usage
%	A = MakeMatrix(n, p, type)
%   A = MakeMatrix(d, n, type)
% Input
%	n           number of rows
%	p           number of columns
%	type        one of: USE (Uniform Spherical Ensemble), fourier, hadamard,
%	            randomortho, randomshperical, uniformrandomprojection.
%                       
% Outputs
%	 A          nxp test matrix
%

if nargin < 3
    type = 'USE';
end

if strcmp(type, 'USE')
    A = UniformSphericalMatrix(n,p);
elseif strcmp(type, 'fourier')
    A = PartialFourierMatrix(n,p);
elseif strcmp(type, 'hadamard')
    A = PartialHadamardMatrix(n,p);
elseif strcmp(type, 'randomortho')
    A = RandomOrthoMatrix(n,p);
elseif strcmp(type, 'randomspherical')
    A = RandomSphericalMatrix(n,p);
elseif strcmp(type, 'uniformrandomprojection')
    A = UniformRandomProjectionMatrix(n,p);
end
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
