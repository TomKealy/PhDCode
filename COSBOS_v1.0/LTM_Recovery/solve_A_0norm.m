% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function A=solve_A_0norm(X,Y)
%% solve Y=AX for A, while minimizing L0-norm of vec(A)

addpath('Lib/SparseLab2.1-Core/Solvers');

m1=size(X,1);
m2=size(Y,1);

XX=kron(X',eye(m2,m2));
YY=Y(:);
AA = SolveOMP(XX, YY, m1*m2);
A=reshape(AA,[m2 m1]);
