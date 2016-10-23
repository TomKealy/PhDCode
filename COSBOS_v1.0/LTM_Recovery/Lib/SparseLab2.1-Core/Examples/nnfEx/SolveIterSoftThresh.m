function [xhat, iter] = SolveIterSoftThresh(A, y, p, maxiters, fullpath, opttol, gamma, trace)
% IterSoftThreshOp -- Apply Iterative Soft Thresholding
%  Usage 
%    [xhat, iter] = IterSoftThreshOp(A, y, p, maxiters, fullpath, opttol, gamma, trace)
%  Input
%    y          nx1 observation vector
%    p          solution dimension
%    maxiters   maximum number of iterations
%    fullpath   0 returns the solution, 1 returns the full path
%    opttol     error tolerance, default is 1e-4
%  Outputs 
%    xhat       solution path
%    iter       number of iterations
%  

if nargin < 4
    maxiters = 100;
end

if nargin < 5
    fullpath = 0;
end

if nargin < 6
    opttol = 1e-6;
end

if nargin < 7
    gamma = 1/2;
end

if nargin < 8
    trace = 1;
end

xk = zeros(p,1);
rho = 1; 

if fullpath
    xhat = [];
end

fastop = (ischar(A) || isa(A, 'function_handle'));

iter = 0;
done = 0;
while ~done
    
    if (fastop)
        rk = y - feval(A,1,p,xk);
        ck = feval(A,2,p,rk); 
    else
        rk = y - A*xk;
        ck = A'*rk;
    end
    
    nrk = norm(rk);
    nck = norm(ck);
    
    if iter == 0
        nc0 = nck;
    end
    
    % threshold
    if iter == 0
        lambda_k = max(abs(ck));
    else
        lambda_k = lambda_k./(1+gamma);
    end

    % adjustment
    ak = SoftThresh(ck, lambda_k);

    % update
    xk = xk + ak;
    
    iter = iter + 1;

    if fullpath
        xhat = [xhat xk(:)];
    end

    % check stopping condition
    if (nck/nc0 < opttol)
        done = 1;
    end

    if iter > maxiters
        done = 1;
    end

end

if ~fullpath
    xhat = xk(:);
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
