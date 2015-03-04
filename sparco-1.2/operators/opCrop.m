function op = opCrop(opIn, m,n, mOffset, nOffset, scale)
% OPCROP   Crop rectangular region from given operator
%
%    OPCROP(OP,M,N,MOFFSET,NOFFSET,SCALE) creates an M by N
%    operator that restricts OP to the M by N rectangular region
%    starting at (zero-based) row MOFFSET and column
%    NOFFSET. Parameter SCALE specifies a scalar scaling factor for
%    the resulting operator. Note that the offsets can be negative
%    and that the specified region can extend beyond the size of
%    the given operator.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opCrop.m 1040 2008-06-26 20:29:02Z ewout78 $

[mOp,nOp] = opsize(opIn);
if nargin < 2, m = mOp; end
if nargin < 3, n = nOp; end;
if nargin < 4, mOffset = 0; end;
if nargin < 5, nOffset = 0; end;
if nargin < 6, scale   = 1; end;

% Check arguments
if (m <= 0) || (n <= 0)
    error('Size of crop operator must be positive');
end

op = @(x,mode) opCrop_intrnl(opIn,m,n,mOffset,nOffset,mOp,nOp, ...
                             scale, x,mode);


function y = opCrop_intrnl(opIn,m,n,mOffset,nOffset,mOp,nOp, scale,x,mode)
checkDimensions(m,n,x,mode);
if mode == 0
   info = opIn([],0);
   if ~isreal(scale), info{3} = ones(1,4); end;
   y = {m,n,info{3},{'Crop',opIn,mOffset,nOffset,scale}};
   return;
end

% Compute padding sizes and offsets
mPrepad  = max(0,mOffset);
mPostpad = max(0,mOp - (m + mOffset));
mOverlap = max(0,mOp - mPrepad - mPostpad);
nPrepad  = max(0,nOffset);
nPostpad = max(0,nOp - (n + nOffset));
nOverlap = max(0,nOp - nPrepad - nPostpad);

if ((mOffset + m <= 0) || (nOffset + n <= 0) || ...
       (mOffset >= mOp) || (nOffset >= nOp))
       % Crop has no overlap with operator, result y will be 0
   if mode == 1, y = zeros(m,1); else y = zeros(n,1); end;
   return
end

if mode == 1
   % Set initial solution
   y = zeros(m,1);
   
   % Set modified x
   x = [zeros(nPrepad,1);
        x(max(0,-nOffset) + (1:nOverlap)); ...  % Overlap
        zeros(nPostpad,1)];

   % Apply operator to x and scale       
   z = scale * opIn(x,1);
       
   % Copy the right part of z into y 
   y(max(0,-mOffset) + (1:mOverlap)) = z(max(0,mOffset)+(1:mOverlap));
else
   % Set initial solution
   y = zeros(n,1);
   
   % Set modified x
   x = [zeros(mPrepad,1);
        x(max(0,-mOffset) + (1:mOverlap)); ...  % Overlap
        zeros(mPostpad,1)];

   % Apply operator to x and scale       
   z = conj(scale) * opIn(x,2);
       
   % Copy the right part of z into y 
   y(max(0,-nOffset) + (1:nOverlap)) = z(max(0,nOffset)+(1:nOverlap));
end
