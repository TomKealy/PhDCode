function interval_graph(g,ilist)
% interval_graph(g,ilist) --- create an interval graph
% ilist is an n-by-2 list of intervals.

[n,x] = size(ilist);

for i=1:n
    ilist(i,:) = sort(ilist(i,:));
end

[ilist,idx] = sortrows(ilist);

resize(g,n);
clear_edges(g);
rmxy(g);

for i=1:n-1
    a = ilist(i,:);
    for j=i+1:n
        b = ilist(j,:);
        if (meets(a,b))
            add(g,i,j);
        end
        if a(2) < b(1)
            break
        end
    end
end

renumber(g,idx);


function yn = meets(a,b)

if (a(2) < b(1)) || (b(2) < a(1))
    yn = 0;
else
    yn = 1;
end
