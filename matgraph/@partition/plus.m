function r = plus(p,q)
% plus --- p+q is the join of p and q

n1 = nv(p);
n2 = nv(q);

if n1 ~= n2
    error('+ only defined for partitions on the same ground set');
end

r = q;

n = n1;

for i=1:n-1
    idx = find(p.array(:,i));
    irow = p.array(idx,:);
    for j=find(irow)
        r=merge(r,i,j);
    end
end
r.array = sortrows(r.array);