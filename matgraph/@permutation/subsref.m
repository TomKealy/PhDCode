function j = subsref(p,S)
% p(k) returns the result of p applied to k. 
% If k is out of bounds, 0 is returned.

n  = length(p);
args = S.subs;
a = args(1);
i = a{1};
if (i<1) || (i>n)
    j = 0;
else
    j = p.array(i);
end
