function cartesian(g,h1,h2)
% cartesian(g,h1,h2) --- overwrite g with the product of h1 and h2

n1 = nv(h1);
n2 = nv(h2);
n = n1 * n2;

resize(g,n);
clear_edges(g);

for u=1:n-1
    for v=u+1:n
        [u1,u2] = splitout(u,n1,n2);
        [v1,v2] = splitout(v,n1,n2);
        
        if (u1==v1) && has(h2,u2,v2)
            add(g,u,v);
        end
        if has(h1,u1,v1) && (u2==v2)
            add(g,u,v);
        end
    end
end


if hasxy(h1) && hasxy(h2)
    xy1 = getxy(h1);
    xy2 = getxy(h2);
    
    xy = zeros(n,4);
    
    for u=1:n
        [u1,u2] = splitout(u,n1,n2);
        xy(u,:) = [xy1(u1,:),xy2(u2,:)];
    end
    m = randn(4,2);
    [q,r] = qr(m);
    p = q(:,1:2);
    xy = xy*p;
    embed(g,xy);
    
end


if is_labeled(h1) && is_labeled(h2)
    for u=1:n
        [u1,u2] = splitout(u,n1,n2);
        lab1 = get_label(h1,u1);
        lab2 = get_label(h2,u2);
        label(g,u, ['(', lab1, ',', lab2, ')']);
    end
end





function [i,j] = splitout(k, n1, n2)

j = mod(k-1,n2);
i = (k-1-j)/n2;
%i = floor((k-1)/n2);

i=i+1;
j=j+1;
