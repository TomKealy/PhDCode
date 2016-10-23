function RegEx03(X,y,b,zn)
% Plots results for Matching Pursuit Algorithm against the true
% underlying model.

%Find sparse solutions, with true model unknown

[bhatMP, iters activationHist] = SolveMP(X, y, 100, norm(zn), 0);

% Plot Results

figure;
plot(1:200,bhatMP,'-');
hold on
plot(1:200,b,'.r')
title('Matching Pursuit Regression Solution')
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
