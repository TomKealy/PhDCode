function y=deci2bin1(x,l)
% Converts a given decimal number x into a binary number of l bits
if x==0, y=0;
 else y=[];
      while x>=1,  y=[rem(x,2) y]; x=floor(x/2);  end
end
if nargin>1, y=[zeros(size(x,1),l-size(y,2)) y]; end