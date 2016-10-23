%
% Script for Reconstruction Example
%
% solves the problem
% min ||D(x)||_1 s.t. y = S(x),
%
% where D is successive differences 
% and S random sampling
%
clear all;

% Signal
p = 2^9;
sname = 'HeaviSine';
x = MakeSignal(sname,p); 
x = x(:);
a = min(x)-1; b = max(x)+1;

% random sampling
r = randperm(p);
n = p/4;
global I;
I = sort(r(1:n));
y = x(I);
    
% Solve
tic
[xhat_tilda, iters] = SolveIST('FastOp',y,p);
xhat = TransformSparsify(1,p,xhat_tilda,'Difference');
toc

figure; 
plot(x,'-.r','LineWidth',2); hold on;
plot(I,y,'g.');
plot(xhat);
axis([1 p a b]);
legend('signal','samples','reconstruction');
title('IST reconstruction');

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
