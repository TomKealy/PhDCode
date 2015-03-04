function tf = eq(p,q)
% test if p == q
np = length(p);
nq = length(q);

if np ~= nq
    tf = false;
    return
else
    tf = all(p.array == q.array);
end
