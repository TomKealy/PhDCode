% This example creates a noisy version of a true underlying model, for
% SparseLab tools to recover it.

p = 200; %number of variables
n = 50;  %number of observations
k = 10;  %number of variables in the true model

z = randn(n,1);
zn = z*2; %Noise level drawn from N(0,4)

X = MakeMatrix(n,p,'USE');

% Choose random locations for the signal variables

pos=randperm(200); 
pos=pos(1:k); 
b=zeros(p,1);

%Build Model

b(pos)=10^2*rand(k,1); %true variables amplified 10^2
y = X*b+zn;%

%
% Copyright (c) 2006. Victoria Stodden
%

% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
