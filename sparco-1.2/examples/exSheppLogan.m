function exSheppLogan(n)
% Example illustrating the use of the ELLIPSOIDS function.
%
% See also: ELLIPSOIDS

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: exSheppLogan.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin < 1, n  = 128; end

figure(1); clf;
[P1,x,y,z] = ellipsoids(n,[0,0,-0.25],[0,0,0]);
surf(x,y,z,P1), shading flat; hold on;

[P2,x,y,z] = ellipsoids(n,[0,0,-0.25],[0,90,0]);
surf(x,y,z,P2), shading flat;

[P3,x,y,z] = ellipsoids(n,[0,0,-0.25],[90,0,0]);
surf(x,y,z,P3), shading flat;

[P4,x,y,z] = ellipsoids(n,[0.0,0.3,-0.25],[35,0,20]);
h = surf(x,y,z,P4); shading flat;
set(h,'FaceAlpha',0.2);
colormap gray; hold off;


figure(2); clf;
subplot(2,2,1); imagesc(P1), colormap gray;
subplot(2,2,2); imagesc(P2), colormap gray;
subplot(2,2,3); imagesc(P3), colormap gray;
subplot(2,2,4); imagesc(P4), colormap gray;


