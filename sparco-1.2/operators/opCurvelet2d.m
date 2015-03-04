function op = opCurvelet2d(m,n,nbscales,nbangles,ttype)
% OPCURVELET2D  Two-dimensional curvelet operator
%
%    OPCURVELET2D(M,N,NBSCALES,NBANGLES,TTYPE) creates a
%    two-dimensional curvelet operator for M by N matrices. The
%    curvelet transform is computed using the Curvelab code.
%
%    The remaining three parameters are optional; NBSCALES gives
%    the number of scales and is set to ceil(log2(min(M,N)) - 3)
%    by default, as suggested by Curvelab. NBANGLES gives the number
%    of angles at the second coarsest level which must be a
%    multiple of four with a minimum of 8. By default NBANGLES is
%    set to 16. TTYPE determines the type of transformation and is
%    set to 'WRAP' by default.
%
%    See also CURVELAB.

%   Copyright 2008, Gilles Hennenfent, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opCurvelet2d.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin < 3, nbscales = ceil(log2(min(m,n)) - 3); end;
if nargin < 4, nbangles = 16;                       end;
if nargin < 5, ttype = 'WRAP';                      end; % Wrapping

finest  = 1;
is_real = 1;

% Compute length of curvelet coefficient vector
if strcmp(ttype,'ME')
  C = mefcv2(randn(m,n),m,n,nbscales,nbangles);

  hdr{1}{1} = size(C{1}{1});
  cn = prod(hdr{1}{1});
  for i = 2:nbscales
    nw = length(C{i});
    hdr{i}{1} = size(C{i}{1});
    hdr{i}{2} = size(C{i}{nw/2+1});
    cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
  end
else
  C = fdct_wrapping_mex(m,n,nbscales,nbangles,finest,randn(m,n));
  
  hdr{1}{1} = size(C{1}{1});
  cn = prod(hdr{1}{1});  
  for i = 2:nbscales
    nw = length(C{i});
    hdr{i}{1} = size(C{i}{1});
    hdr{i}{2} = size(C{i}{nw/4+1});
    cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
  end
end

parms = {m,n,cn,hdr,finest,nbscales,nbangles,is_real,ttype};
op = @(x,mode) opCurvelet2d_intrnl(parms{:},x,mode);


function y = opCurvelet2d_intrnl(m,n,cn,hdr,ac,nbs,nba,is_real,ttype,x,mode)
checkDimensions(cn,m*n,x,mode);
if mode == 0
   cmplx = ~is_real;
   y = {cn,m*n,[cmplx,cmplx,cmplx,cmplx],{'Curvelet2d'}};
elseif mode == 1
   % Analysis mode
   if strcmp(ttype,'ME')
      y = mefcv2(reshape(x,m,n),m,n,nbs,nba);
   else
      y = fdct_wrapping_mex(m,n,nbs,nba,ac,reshape(x,m,n));
      y = fdct_wrapping_c2r(y);
   end
   y = fdct_c2v(y,cn);
else
   % Synthesis mode  
   if strcmp(ttype,'ME')
      x = mefdct_v2c(x,hdr,nba);
      y = meicv2(x,m,n,nbs,nba);
   else
      x = fdct_v2c(x,hdr,ac,nba);
      x = fdct_wrapping_r2c(x);
      y = ifdct_wrapping_mex(m,n,nbs,nba,ac,x);
   end
   if is_real
      y = real(y);
   end
   y = y(:);
end
