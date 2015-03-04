function op = opPadding(sizeLHS, sizeRHS)
% OPPADDING   Padding operator
%
%   OPPADDING(SIZELHS,SIZERHS) create an operator that pads and
%   unpads matrices equally around each size. SIZELHS and SIZERHS
%   give the matrix size before and after padding. SIZERHS must
%   be greater than or equal to SIZELHS, and both must have an even
%   number of rows and columns.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opPadding.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opPadding_intrnl(sizeLHS,sizeRHS,x,mode);


function y = opPadding_intrnl(sizeLHS,sizeRHS,x,mode)
m = prod(sizeLHS);
n = prod(sizeRHS);

checkDimensions(m,n,x,mode);
if mode == 0
   y = {m,n,[0,1,0,1],{'Padding'}};
elseif mode == 1
   z    = reshape(x,sizeRHS);
   idx1 = (sizeRHS(1)/2 - (sizeLHS(1)/2 - 1)):(sizeRHS(1)/2 + sizeLHS(1)/2);
   idx2 = (sizeRHS(2)/2 - (sizeLHS(2)/2 - 1)):(sizeRHS(2)/2 + sizeLHS(2)/2);
   y    = reshape(z(idx1,idx2),[prod(sizeLHS),1]);
else
   z    = zeros(sizeRHS);
   idx1 = (sizeRHS(1)/2 - (sizeLHS(1)/2 - 1)):(sizeRHS(1)/2 + sizeLHS(1)/2);
   idx2 = (sizeRHS(2)/2 - (sizeLHS(2)/2 - 1)):(sizeRHS(2)/2 + sizeLHS(2)/2);
   z(idx1,idx2) = reshape(x,sizeLHS);
   y            = reshape(z,[prod(sizeRHS),1]);
end
