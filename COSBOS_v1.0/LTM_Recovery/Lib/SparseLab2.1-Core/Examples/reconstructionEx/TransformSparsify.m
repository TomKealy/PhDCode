function y = TransformSparsify(mode, n, x, transform)

if nargin < 3
    transform = 'Wavelet';
end

switch transform
    case 'Wavelet'    
        qmf = MakeONFilter('Haar');
        s = 3;
        X = reshape(x,sqrt(n),sqrt(n));
        if mode == 1
            Y = IWT2_PO(X,s,qmf);
        elseif mode == 2
            Y = FWT2_PO(X,s,qmf);
        end
        y = reshape(Y,n,1);
    case 'Laplacian'
        L = sparse([1:n,2:n,1:n-1], [1:n,1:n-1,2:n], [4*ones(1,n) -1*ones(1,2*(n-1))]);
        if mode == 1
            y = L*x;
        elseif mode == 2
            y = L'*x;
        end      
    case 'Difference'
        if mode == 1
            for k = 1:n % lower tri
                y(k) = sum(x(1:k));
            end
        elseif mode == 2
            for k = 1:n % upper tri
                y(k) = sum(x(k:end));
            end
        end      
        y = y(:);
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
