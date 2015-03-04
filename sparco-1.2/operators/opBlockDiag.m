function op = opBlockDiag(varargin)
% OPBLOCKDIAG  Operator-diagonal operator
%
%    OPBLOCKDIAG(WEIGHTS, OP1, OP2, ...) creates a compound
%    operator with the given weighted operators on the diagonal;
%    DIAG([WEIGHT1*OP1, WEIGHT2*OP2, ...]). The WEIGHT parameter is
%    optional and can either be a scalar or a vector. In case of a
%    scalar, all operators are weighted equally. If only a single
%    operator is given it is replicated as many times as there are
%    weights.
%
%    See also opFoG, opKron, opWindowedOp, opDictionary.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opBlockDiag.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin > 0 && isnumeric(varargin{1})
  weight = varargin{1};
  opList = varargin(2:end);
else
  weight = 1;
  opList = varargin;
end


if length(opList) == 0
   error('At least one operator must be specified');
else
   % Check operator space and compute size
   m = 0; n = 0; c = [0 0 0 0];
   for i=1:max(length(weight),length(opList))
      opIdx = i; if length(opList) == 1, opIdx = 1; end;
      info = opList{opIdx}([],0);

      m = m + info{1};            % Total number of rows
      n = n + info{2};            % Total number of columns
      c = c | info{3};            % Complex or real operator
   end
   
   % Check weight vector
   weight = weight(:);
   if (~isscalar(weight)) ...
      && (length(opList) > 1) ...
      && (length(weight) ~= length(opList))
     error(['Length of weight vector must be equal to the number ' ...
            'of operators']);
   end
   if any(~isreal(weight)), c = ones(1,4); end;

   op = @(x,mode) opBlockDiag_intrnl(m,n,c,weight,opList,x,mode);
end


function y = opBlockDiag_intrnl(m,n,c,weight,opList, x, mode)
checkDimensions(m,n,x,mode);
if mode == 0
   y = {m,n,c,{'BlockDiag',weight,opList}};
elseif mode == 1
   y  = zeros(m,1);
   kx = 0;
   ky = 0;
   for i=1:max(length(weight),length(opList))
      opIdx = i; if length(opList) == 1, opIdx = 1; end;
      wIdx  = i; if length(weight) == 1, wIdx  = 1; end;
      
      s               = opList{opIdx}([],0);
      y(ky+1:ky+s{1}) = weight(wIdx)*opList{opIdx}(x(kx+1:kx+s{2}),mode);
      kx              = kx + s{2};
      ky              = ky + s{1};
   end;
else   
   y  = zeros(n,1);
   kx = 0;
   ky = 0;
   for i=1:max(length(weight),length(opList))
      opIdx = i; if length(opList) == 1, opIdx = 1; end;
      wIdx  = i; if length(weight) == 1, wIdx  = 1; end;
      
      s               = opList{opIdx}([],0);
      y(ky+1:ky+s{2}) = conj(weight(wIdx))*(opList{opIdx}(x(kx+1:kx+s{1}),mode));
      kx              = kx + s{1};
      ky              = ky + s{2};
   end
end
