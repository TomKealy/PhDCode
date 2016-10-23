% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function A=solve_A_fullrank(X,Y)
%% solve Y=AX for A, while system is overdetermined

if size(X,2)>=size(X,1)
    A=Y/X;
else
    error('Error: Rank is not full. Cannot use full rank recovery of LTM.');
end