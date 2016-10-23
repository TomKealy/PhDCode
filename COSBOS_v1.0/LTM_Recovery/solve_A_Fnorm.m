% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function A=solve_A_Fnorm(X,Y)
%% solve Y=AX for A, while minimizing rank(A), or equivalently, minimizing F-norm of A 

threshold=0.01; % threshold for singular values

l=size(Y,1);
m=size(X,1);
[U,S,V] = svd(X); % X = U*S*V'
nn=sum(diag(S)>threshold); % remove small singular values
S=S(:,1:nn); % remove small singular values
V=V(:,1:nn); % remove small singular values
Y=Y*V;
Sinv=diag(1./diag(S));
Z1=Y*Sinv;
Z=[Z1 zeros(l,m-nn)];
A=Z*U';
