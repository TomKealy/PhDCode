function vtcs = component(g,v)
% component(g,v) -- find vertices in v's component of g

n = nv(g);
if (v < 1) || (v>n)
    error('Vertex out of range')
end

vec = zeros(n,1);
vec(v) = 1;
A = matrix(g);
while(true)
    w = double(A*vec +vec > 0);
    if nnz(w) == nnz(vec)
        break
    end
    vec = w;
end
vtcs = find(vec);