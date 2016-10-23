function y = FastOpLS(mode, m, n, x, I, dim)

if (mode == 1)
    
    u = zeros(dim, 1);
    u(I) = x;
    y = FastOp(1,dim,u);
    
elseif (mode == 2)
    
    x2 = FastOp(2,dim,x);
    y = x2(I);
    
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
