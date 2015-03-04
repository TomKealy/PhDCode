function s = scalarToRGB(x,colors)
% input values are assumed to lie between 0 and 1

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: scalarToRGB.m 1040 2008-06-26 20:29:02Z ewout78 $

l  = size(colors,1);
m  = size(x,1);
n  = size(x,2);
s  = zeros(m,n,3);

for i=1:m
  for j=1:n
     idx = max(1,min(l,1+floor((l-1) * x(i,j))));
     s(i,j,:) = colors(idx,:);
  end
end
