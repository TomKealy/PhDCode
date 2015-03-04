function s = sign(p)
% sign of a permutation (1 for even, -1 for odd)
n = length(p);
if (n==0)
    s = 1;
    return
end

c = cycles(p);
s = 0;
for k=1:length(c)
    s = s + length(c{k}) + 1;
end

s = mod(s,2);
s = (-1)^s;
