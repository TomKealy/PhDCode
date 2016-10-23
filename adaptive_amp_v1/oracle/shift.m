function ys = shift(y,e,bc)

if ~exist('bc','var')
    bc = 'zero';
end;

ys = zeros(size(y));

switch bc

case 'circular'
if e > 0	% shift right
    ys(e+1:end) = y(1:end-e);
    ys(1:e) = y(end-e+1:end);
else		% shift left
    ys(1:end+e) = y(1-e:end);
    ys(end+e+1:end) = y(1:1-e-1);
end;

case 'zero'
if e > 0	% shift right
    ys(e+1:end) = y(1:end-e);
else		% shift left
    ys(1:end+e) = y(1-e:end);
end;

% mirror for the right shift, adjoint for the left
case 'causal-mirror'
if e > 0	% shift right
    ys(e+1:end) = y(1:end-e);
    ys(1:e) = y(e:-1:1);
else		% shift left
    ys(1:end+e) = y(1-e:end);
    ys(1:-e) = ys(1:-e) + y(-e:-1:1);
end;

% mirror for the left shift, adjoint for the right
case 'anticausal-mirror'
if e > 0	% shift right
    ys(e+1:end) = y(1:end-e);
    ys(end-e+1:end) = ys(end-e+1:end) + y(end:-1:end-e+1);
else		% shift left
    ys(1:end+e) = y(1-e:end);
    ys(end+e+1:end) = y(end:-1:end+e+1);
end;

end;
