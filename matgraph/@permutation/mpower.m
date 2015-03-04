function q = mpower(p,t)
% p^t is p multiplied by itself t times. May have t <= 0.
if t==0
    q = permutation(length(p));
    return
end

if t<0
    q = (inv(p))^(-t);
    return
end

if t == 1
    q = p;
	return
end

% recurse when t>1

if mod(t,2) == 0
    q = p^(t/2);
    q = q*q;
else
    q = p^((t-1)/2);
    q = p*q*q;
end
