function op = opKron(op1,op2)
% OPKRON  Kronecker tensor product
%
%    OPKRON(OP1,OP2) creates an operator that is the Kronecker
%    tensor product of OP1 and OP2.

%   Copyright 2008, Rayan Saab, Ewout van den Berg and
%                   Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opKron.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opKron_intrnl(op1,op2,x,mode);


function y = opKron_intrnl(op1,op2,x,mode)
info1 = op1([],0);
info2 = op2([],0);

checkDimensions(info1{1}*info2{1},info1{2}*info2{2},x,mode);

if mode == 0
   c1 = (info1{3} ~= 0); % Ensure values are 0 or 1
   c2 = (info2{3} ~= 0); % Ensure values are 0 or 1
   ck = [c1(1+c2(1)), c1(1+c2(2)), c1(3+c2(3)), c1(3+c2(4))];
   y  = {info1{1}*info2{1}, info1{2}*info2{2}, ck, {'Kron',op1,op2}};
elseif mode == 1
   % Compute y = (op1 @ op2) x
   x = reshape(x,info2{2},info1{2});
   z = zeros(info2{1},info1{2});
   for i=1:info1{2}
      z(:,i) = op2(x(:,i),1);
   end
   z = z';

   y = zeros(info1{1},info2{1});
   for i=1:info2{1}
      y(:,i) = op1(z(:,i),1);
   end
   y = y';
   y = y(:);
else
   % Compute y = (op1 @ op2)' * x
   x = reshape(x,info2{1},info1{1});
   z = zeros(info2{2},info1{1});
   for i=1:info1{1}
      z(:,i) = op2(x(:,i),2);
   end
   z = z';
   
   y = zeros(info1{2},info2{2});
   for i=1:info2{2}
      y(:,i) = op1(z(:,i),2);
   end
   y = y';
   y = y(:);
end
