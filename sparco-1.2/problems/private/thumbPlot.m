function P = thumbPlot(P,x,y,color)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: thumbPlot.m 1040 2008-06-26 20:29:02Z ewout78 $

m = size(P,1);
n = size(P,2);
if (size(P,3) == 0) & (length(color) == 3)
  % Convert to gray-scale
  color = 0.30*color(1) + 0.59*color(2) + 0.11*color(3);
end

mnx = min(x);  % Minimum x
mxx = max(x);  % Maximum x
mny = min(y);  % Minimum y
mxy = max(y);  % Maximum y
dy  = (mxy - mny) * 0.1;   % Offset on vertical axis
sx  = (mxx - mnx) * 1.0;   % Scale of horizontal axis
sy  = (mxy - mny) * 1.2;   % Scale of vertical axis

if (sx < 1e-6), sx = 1; end
if (sy < 1e-6), sy = 1; end

for i=1:length(x)-1
   x0 = floor(1 + (n-1) * (x(i  ) - mnx) / sx);
   x1 = floor(1 + (n-1) * (x(i+1) - mnx) / sx);
   y0 = floor(    (n-1) * (y(i  ) - mny + dy) / sy);
   y1 = floor(    (n-1) * (y(i+1) - mny + dy) / sy);
   
   samples = 1+2*max(abs(x1-x0)+1,abs(y1-y0)+1);
   c       = linspace(0,1,samples);
   idx     = round((1-c)*x0 + c*x1);
   idy     = n - round((1-c)*y0 + c*y1);
   for j=1:samples
      P(idy(j),idx(j),:) = color;
   end
end
