function S = RSTMat(N);
% This function generates an NxN RST matrix S, s.t. writing
% x_hat = S*x is equivalent to x_hat = RST(x)

S = eye(N);

for jj = 1:N
    S(:,jj) = RST(S(:,jj));
end
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
