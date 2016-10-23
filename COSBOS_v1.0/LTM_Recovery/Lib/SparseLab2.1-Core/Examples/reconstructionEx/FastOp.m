function y = FastOp(mode, p, x)

global I;

if mode == 1
    
    z = TransformSparsify(mode, p, x, 'Difference');
    y = z(I);

elseif mode == 2

    z = zeros(p,1);
    z(I) = x;
    y = TransformSparsify(mode, p, z, 'Difference');
    
end

%
% Copyright (c) 2006. Iddo Drori
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
