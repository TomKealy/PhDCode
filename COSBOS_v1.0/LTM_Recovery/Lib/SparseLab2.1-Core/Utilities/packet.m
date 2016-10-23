function p=packet(d,b,n)
% packet -- Packet table indexing
%  Usage
%    p = packet(d,b,n)
%  Inputs
%    d     depth of splitting in packet decomposition
%    b     block index among 2^d possibilities at depth d
%    n     length of signal
%  Outputs
%    p     linear indices of all coeff's in that block
%

npack = 2^d;
p =  ( (b * (n/npack) + 1) : ((b+1)*n/npack ) ) ;

    
%
% Copyright (c) 2006. David Donoho
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
