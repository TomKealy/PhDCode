function yn = has_path(g,u,v)
% has_path(g,u,v) --- determine if there is a path from u to v in g
% this is implemented a bit differently from find_path

global GRAPH_MAGIC

n = nv(g);

if (u<1) || (u>n) || (v<1) || (v>n)
    yn = 0;
    return
end

if u==v
    yn = 1;
    return
end

x = zeros(n,1);
x(u) = 1;

last_sum = 0;

while last_sum < sum(x)
    last_sum = sum(x);
    x = GRAPH_MAGIC.graphs{g.idx}.array * x + x;
    x = double(x>0);
    if (x(v))
        yn = 1;
        return
    end
end

yn = 0;
