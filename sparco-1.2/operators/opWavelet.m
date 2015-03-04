function op = opWavelet(m,n, family, filter, levels, type)
% OPWAVELET  Wavelet operator
%
%    OPWAVELET(M,N,FAMILY,FILTER,LEVELS,TYPE) creates a wavelet
%    operator of given FAMILY, for M by N matrices. The wavelet
%    transformation is computed using the Rice Wavelet Toolbox.
%
%    The remaining three parameters are optional. FILTER = 8
%    specifies the filter length and must be even. LEVELS = 5 gives
%    the number of levels in the transformation. Both M and N must
%    be divisible by 2^LEVELS. TYPE = 'min' indicates what type of
%    solution is desired; 'min' for minimum phase, 'max' for
%    maximum phase, and 'mid' for mid-phase solutions. 

%   Copyright 2008, Rayan Saab, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opWavelet.m 1040 2008-06-26 20:29:02Z ewout78 $

if (nargin < 3), family = 'Daubechies'; end;
if (nargin < 4), filter =     8;        end;
if (nargin < 5), type   = 'min';        end;
if (nargin < 6), levels =     5;        end;

family = lower(family);

switch family
 case {'daubechies'}
    h = daubcqf(filter);
    
 case {'haar'}
    h = daubcqf(0);
    
 otherwise
    error(sprintf('Wavelet family %s is unknown', family));
end

op = @(x,mode) opWavelet_intrnl(m,n,family,filter,levels,type,h,x,mode);


function y = opWavelet_intrnl(m,n,family,filter,levels,type,h,x,mode)
checkDimensions(n*m,n*m,x,mode);
if mode == 0
   y = {n*m,n*m,[0,1,0,1],{'Wavelet',family,filter,levels,type}};
elseif mode == 1
   if isreal(x)
     [y,l] = midwt(reshape(x,[m,n]), h, levels);
   else
     [y1,l] = midwt(reshape(real(x),[m,n]), h, levels);
     [y2,l] = midwt(reshape(imag(x),[m,n]), h, levels);
     y = y1 + sqrt(-1) * y2;    
   end
   y = reshape(y,[m*n,1]);
else
   if isreal(x)
      [y,l] = mdwt(reshape(x,[m,n]), h, levels);
   else
      [y1,l] = mdwt(reshape(real(x),[m,n]), h, levels);
      [y2,l] = mdwt(reshape(imag(x),[m,n]), h, levels);
     y = y1 + sqrt(-1) * y2;    
   end   
   y = reshape(y,[m*n,1]);
end
