% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function Iq=lightDistribution(theta)
%% The luminous intensity distribution of our Vivia 7DR3-RGB fixture (from Renaissance Lighting)
%    theta: angle to normal direction
%    Iq: luminous intensity

a=0:5:90;
I=[509 505 456 398 333 269 203 142 91 49 15 1 0 0 0 0 0 0 0];

Iq = interp1(a,I,theta/pi*180,'pchip');
