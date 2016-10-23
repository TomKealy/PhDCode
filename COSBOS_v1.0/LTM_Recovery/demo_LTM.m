% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

%% demo 1: overdetermined, full rank recovery with pseudo inverse

clear;clc;close all;
load 'Data/U_18565.mat';
X=bsxfun(@minus,TestLight(2:end,:),TestLight(1,:));
Y=bsxfun(@minus,cdata(2:end,:),cdata(1,:));
X=X';
Y=Y';
A_fullrank=solve_A_fullrank(X,Y);

%% demo 2-4: underdetermined, find A0 first

load 'Data/0_30876.mat';
X=bsxfun(@minus,TestLight(2:end,:),TestLight(1,:));
Y=bsxfun(@minus,cdata(2:end,:),cdata(1,:));
X=X';
Y=Y';
A0=solve_A_fullrank(X,Y);

N2=20; % assume only N2 perturbation patters 

%% demo 2: underdetermined, low rank recovery

load 'Data/U_18565.mat';
X=bsxfun(@minus,TestLight(2:end,:),TestLight(1,:));
Y=bsxfun(@minus,cdata(2:end,:),cdata(1,:));
X=X(1:N2,:)';
Y=Y(1:N2,:)';
Z=A0*X-Y;

E=solve_A_Fnorm(X,Z);
A_Fnorm=A0-E;

%% demo 3: underdetermined, L0 recovery

E=solve_A_0norm(X,Z);
A_0norm=A0-E;

%% demo 4: underdetermined, L1 recovery

E=solve_A_1norm(X,Z);
A_1norm=A0-E;

