function s = scaleImage(x,m,n);
% SCALEIMAGE  Basic image rescaling
%
%    SCALEIMAGE(X,M,N) resamples image X based on a 2M by 2N
%    grid. The scaled image is formed by taking the average between
%    pairs of grid points.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: scaleImage.m 1040 2008-06-26 20:29:02Z ewout78 $

idx = round(linspace(1,size(x,1),2*m));
s   = (x(idx(1:2:2*m),:) + x(idx(2:2:2*m),:)) / 2;
idx = round(linspace(1,size(x,2),2*n));
s   = (s(:,idx(1:2:2*n)) + s(:,idx(2:2:2*n))) / 2;
