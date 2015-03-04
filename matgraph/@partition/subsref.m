function vec = subsref(p,S)
% subsref --- implements the notation p(v) and p(v,w)
% p(v) returns the part containing v
% p(v,w) returns 1 if v and w are in the same part and 0 otherwise


args = S.subs;
[m,n] = size(p.array);

% single subscript case

if length(args)==1
    v = args(1);
    v = v{1};
    if (v<1) || (v>n)
        vec = [];
    else
        r = find(p.array(:,v));
        vec = find(p.array(r,:));
    end
    return
end

% two subscript case

v = args(1);
v = v{1};
w = args(2);
w = w{1};

if (v<1) || (v>n) || (w<1) || (w>n) 
    vec = false;
    return
end

r1 = find(p.array(:,v));
r2 = find(p.array(:,w));
vec = r1==r2;
