function mask = scratchMask(n,m,nscratch);
%SCRATCHMASK  Generation of a scratch mask
%
%   SCRATCHMASK(N,M,NSCRATCH) creates an N by M binary mask with
%   NSCRATCH linear scratches originating from one of the four
%   sides of the image and ending in one of the other three.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: scratchMask.m 1040 2008-06-26 20:29:02Z ewout78 $

%p =  n / (n+m);
p = 0.5;

mask = zeros(n,m);

for i=1:nscratch

   while (1),
      if rand(1,1) < p
         x1 = rand(1,1) * 2 - 1;
         y1 = sign(rand(1,1)-0.5);
      else
         y1 = rand(1,1) * 2 - 1;
         x1 = sign(rand(1,1)-0.5);
      end

      if rand(1,1) < p
         x2 = rand(1,1) * 2 - 1;
         y2 = sign(rand(1,1)-0.5);
      else
         y2 = rand(1,1) * 2 - 1;
         x2 = sign(rand(1,1)-0.5);
      end

      if x1~=x2 & y1~=y2, break; end
   end

   x1 = round((m-1) * (1+x1) / 2);
   y1 = round((n-1) * (1+y1) / 2);
   x2 = round((m-1) * (1+x2) / 2);
   y2 = round((n-1) * (1+y2) / 2);
   
   if abs(y1-y2) > abs(x2-x1)
      % Mostly vertical line
      w = linspace(0,1,abs(y2-y1)+1);
      idx = 1 + (y1:sign(y2-y1):y2) + n * round((1-w) * x1 + w * x2);
   else
      % Mostly horizontal line
      w = linspace(0,1,abs(x2-x1)+1);
      idx = 1 + round((1-w) * y1 + w * y2) + n * (x1:sign(x2-x1):x2);
   end
   mask(idx) = 1; % Use mask(idx) = mask(idx) + 1; to determine line distribution
end
