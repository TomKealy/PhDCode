function trim(g,d)
% trim(g,d) --- iteratively remove vertices of degree at most d

if nargin==1
    d = 0;
end

ds = deg(g);
x = find(ds <= d);
x=x';
while (length(x) > 0)
    if (nv(g) == 0)
        return
    end
    delete(g,x);
    ds = deg(g);
    x = find(ds <= d);
    x=x';
end
