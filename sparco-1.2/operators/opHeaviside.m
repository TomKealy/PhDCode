function op = opHeaviside(n,s)
% OPHEAVISIDE  Heaviside operator
%
%    OPHEAVISIDE(N,S) creates an operator for multiplication by an
%    N by N Heaviside matrix. These matrices have ones below and on
%    the diagonal and zeros elsewhere. S is a flag indicating
%    whether the columns should be scaled to unit Euclidean
%    norm. By default the columns are unnormalized.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opHeaviside.m 1040 2008-06-26 20:29:02Z ewout78 $

if (nargin < 2), s = 0; end;

op = @(x,mode) opHeaviside_intrnl(n,s,x,mode);


function y = opHeaviside_intrnl(n,s,x,mode)
checkDimensions(n,n,x,mode);
if mode == 0
   y = {n,n,[0,1,0,1],{'Heaviside',n}};
elseif mode == 1
   % Scale if normalized columns requested
   if s ~= 0, x = [1./sqrt(n:-1:1)]'.*x(:); end

   y = cumsum(x);
else
   y        = cumsum(x);
   ym       = y(end);
   y(2:end) = ym - y(1:end-1);
   y(1)     = ym;
   
   % Scale if normalized columns requested
   if s ~= 0, y = [1./sqrt(n:-1:1)]'.*y(:); end
end
