function op = opWindowedOp(m,operator,window,overlap)
% OPWINDOWEDOP  Overcomplete windowed operator
%
%    OPWINDOWEDOP(M,OPERATOR,WINDOW,OVERLAP) creates an
%    overcomplete windowed operator with overlap:
%
%      /  +------------+-------------------------+
%     |   |            |                         |
%     |   | Diag(W)*OP |                         |
%     |   |            +------------+            | --
%  M  /   +------------+            |            | -- OVERLAP
%     \   |            | Diag(W)*OP |            |
%     |   |            |            +------------+
%     |   |            +------------+            |
%     |   |                         | Diag(W)*OP |
%      \  +-------------------------+------------+
%    
%    with W for WINDOW and OP for OPERATOR.
%
%    Parameter OVERLAP must be strictly less than the length of
%    WINDOW, which in turn must coincide with the number of rows
%    in the OPERATOR. M indicates the signal length and prescribes
%    the number of rows in the resulting operator. The number of
%    columns is always an integer multiple of the number of columns
%    in OPERATOR and depends on the required number of blocks.

%   Copyright 2008, Rayan Saab, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opWindowedOp.m 1040 2008-06-26 20:29:02Z ewout78 $

info   = operator([],0);
window = window(:);

% Window length should coincide with number of rows in the operator
if length(window) ~= info{1}
  error('Window length and operator size do not agree');
end

if overlap >= info{1}
  error('Overlap must be strictly less than the window length');
end

% Compute the number of blocks
offset = length(window) - overlap;
blocks = max(1,ceil((m - length(window)) / offset)+1);

% Create the operator
op = @(x,mode) opWindowedOp_intrnl(m,operator,window,offset,blocks,info,x,mode);


function y = opWindowedOp_intrnl(m,operator,window,offset,blocks,info,x,mode)
checkDimensions(m,blocks*info{2},x,mode);

if mode == 0
    type = info{3};
    if isreal(window), c = 0; else c = 1; end;
    c = [c | type(1), c | type(2), type(3+c), type(4)];
    y = {m, blocks*info{2}, c, {'WindowedOp',operator}};
elseif mode == 1
    y = zeros(info{1} + (blocks-1) * offset, 1);
    for i=0:blocks-1
      idx     = i*offset + (1:info{1});
      xWindow = x(i*info{2}+(1:info{2}));
      y(idx)  = y(idx) + window.*operator(xWindow,1);
    end
    y = y(1:m);
else
    x = [x; zeros(info{1} + (blocks-1) * offset - m, 1)];
    y = zeros(info{2} * blocks, 1);
    for i=0:blocks-1
      idx     = i*info{2}+(1:info{2});
      xWindow = x(i*offset + (1:info{1}));
      y(idx)  = operator(window.*xWindow,2);
    end
end
