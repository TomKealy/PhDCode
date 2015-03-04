function r = mtimes(p,q)
% mtimes(p,q) --- p*q is the meet of the partitions p and q

n1 = nv(p);
n2 = nv(q);

if n1 ~= n2
    error('* only defined for partitions on the same ground set')
end

m1 = np(p);
m2 = np(q);

mat = sparse(0,0);

for i=1:m1
    for j=1:m2
        a = p.array(i,:) & q.array(j,:);
        if sum(a) > 0
            mat = [mat;a];
        end
    end
end

mat = sortrows(mat);
mat = unique(mat,'rows');
r.array = mat;

r = class(r,'partition');