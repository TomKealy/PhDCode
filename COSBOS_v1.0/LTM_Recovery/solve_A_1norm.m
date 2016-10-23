% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function A=solve_A_1norm(X,Y)
%% solve Y=AX for A, while minimizing L1-norm of vec(A) 

addpath('Lib/l1magic/Optimization');

m1=size(X,1);
m2=size(Y,1);

XX=kron(X',eye(m2,m2));
YY=Y(:);
evalc('AA = l1eq_pd(zeros(m2*m1,1), XX, [], YY)'); % avoid echo 
A=reshape(AA,[m2 m1]);
