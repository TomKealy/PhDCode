% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

%% This is the demo for the light reflection model with ceiling-mounted sensors

clear;clc;close all;

addpath('../LTM_Recovery');

%% get A0
load 'Data/0_24323.mat';
X=bsxfun(@minus,TestLight(2:end,:),TestLight(1,:));
Y=bsxfun(@minus,cdata(2:end,:),cdata(1,:));
X=X';
Y=Y';
A0=solve_A_fullrank(X,Y);

%% get A
load 'Data/A_17327.mat';
X=bsxfun(@minus,TestLight(2:end,:),TestLight(1,:));
Y=bsxfun(@minus,cdata(2:end,:),cdata(1,:));
X=X';
Y=Y';
A=solve_A_fullrank(X,Y);

%% generate all kernels
coordinates;
para=1;
K=generateAllKernels(lights,sensors,dim,para); % slow, better store K in a mat file

%% get floor-plane occupancy map
E=A0-A;
E(E<0)=0;

C=zeros(dim(1),dim(2)); % floor-plane occupancy map
sumK=zeros(dim(1),dim(2)); % for normalization

lambda1=1; % see Eq. (16) in [1]
lambda2=1; % see Eq. (16) in [1]

for s=1:size(sensors,1)
    for l=1:size(lights,1)
        a=E(4*s-3,3*l-2)+E(4*s-2,3*l-1)+E(4*s-1,3*l);
        a=a^lambda1;
        C=C+a*K{s,l};
        sumK=sumK+K{s,l};
    end
end

%% visualize the floor-plane occupancy
C=C./(sumK.^lambda2);
imagesc(C);
axis equal off;
title('floor-plane occupancy');

