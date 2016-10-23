function [xhat,errornorm] = SolveISTBlock(A, y, x, L)
%function xhat = SolveISTBlock(A, y)
% IterSoftThresh -- Apply Iterative Soft Thresholding 
%                   Block variant with least squares projection
%  Usage 
%    xhat = SolveISTBlock(A, y)
%  Inputs 
%    A     dxn matrix
%    y     dx1 observation vector
%  Outputs 
%    xhat         Solution
%    xhat_b       Solution of block variant with least squares projection

[d,n] = size(A);
At = A';

if nargin < 4
    L = 20;
end
s = n/L;
Bstarj = zeros(L,s,d);
r = randperm(n);
for j = 1:L
    B = A(:,r(1+(j-1)*s:j*s));
    Bstar = inv(B'*B)*B';
    Bstarj(j,:,:) = Bstar;
end

xt = zeros(n,1);
xt_b = zeros(n,1);
xtj_b = zeros(L,s);
iters = 15;
M = sqrt(2*log(n)/d);

for t = 1:iters

    % residual
    rt = y - A*xt;
    
    % correlations
    ct = At*rt;
    
    % threshold
    if t == 1
        sct = sort(abs(ct),'descend');
        lambdat = sct(10);
    else
        lambdat = lambdat/2;
    end

    % adjustment
    at = SoftThresh(ct,lambdat);

    % update
    rho = 1;
    xt = xt + rho*at;

    % block variant
    errornorm(t) = norm(x - xt_b);
    for j = 1:L
        Bstar = squeeze(Bstarj(j,:,:));
        u = xtj_b(j,:)' + rho*SoftThresh((Bstar*(y-A*xt_b)),lambdat);
        xtj_b(j,:) = u';
        xt_b(r(1+(j-1)*s:j*s)) = xtj_b(j,:);
    end

end

xhat = xt_b(:);
    
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
