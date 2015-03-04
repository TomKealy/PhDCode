function yn = eq(p,q)
% eq(p,q) --- p==q checks if two partitions are equal

if ~all(size(p)==size(q))
    yn = false;
    return
end

yn = all(all(p.array==q.array));
yn = full(yn);