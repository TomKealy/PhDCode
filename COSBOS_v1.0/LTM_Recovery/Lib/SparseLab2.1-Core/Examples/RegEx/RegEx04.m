function RegEx04(X,y,b,zn)
% Plots results for Orthogonial Matching Pursuit Algorithm against the true
% underlying model.

%Find sparse solutions, with true model unknown

[bhatOMP, iters activationHist] = SolveOMP(X, y, 100, norm(zn), 0);

% Plot Results

figure;
plot(1:200,bhatOMP,'-');
hold on
plot(1:200,b,'.r')
title('Orthogonal Matching Pursuit Regression Solution')
xlabel('Variable Number')
ylabel('Coefficient Value')
legend('Estimated Coefficient Values','True Coefficient Values')%

%
% Copyright (c) 2006. Victoria Stodden
%

% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
