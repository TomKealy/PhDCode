function q = merge(p,x,y)
% merge(p,x,y) --- form a partition from p by combining x and y's parts.
% If x or y are out of range, p is returned.

q = p;
n = nv(p);

if (x<1) || (x>n) || (y<1) || (y>n) || (x==y)
    return
end

rx = find(p.array(:,x));
ry = find(p.array(:,y));

if rx == ry
    return
end

% merge row ry into row rx
q.array(rx,:) = q.array(rx,:) | q.array(ry,:);

% remove row ry

idx = [1:ry-1,ry+1:np(p)];
q.array = q.array(idx,:);
q.array = sortrows(q.array);
