% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

%% This is the demo for the light blockage model with wall-mounted sensors

clear;clc;close all;

compile; % compile the two cpp files (you need MATLAB compiler)
addpath('../LTM_Recovery');

%% get A0
load 'Data/0_30876.mat';
X=bsxfun(@minus,TestLight(2:end,:),TestLight(1,:));
Y=bsxfun(@minus,cdata(2:end,:),cdata(1,:));
X=X';
Y=Y';
A0=solve_A_fullrank(X,Y);

%% get A
load 'Data/U_85164.mat';
X=bsxfun(@minus,TestLight(2:end,:),TestLight(1,:));
Y=bsxfun(@minus,cdata(2:end,:),cdata(1,:));
X=X';
Y=Y';
A=solve_A_fullrank(X,Y);

%% get volume
E=A0-A;
E(E<0)=0;
sigma=20;
coordinates;
H=hashGaussians(sensors,lights,dim,sigma); % once computed, H can be used many times
V=volumeFromHashing(sensors,lights,dim,H,E);
V=V(:,end:-1:1,:); % to be consistent with occupancy scenarios in the papers

%% visualize the floor-plane occupancy
FloorPlane=sum(V,3);
imagesc(FloorPlane);
axis equal off;
title('floor-plane occupancy');

%% save result to TIFF 3D image (can be viewed with Amira, ImageJ, etc.)
writeTiff(V,'V.tif');

