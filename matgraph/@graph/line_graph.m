function line_graph(g,h)
% line_graph(g,h) --- set g to be the line graph of h
% The line graph of h is the intersection graph of its edges.

elist = sortrows(edges(h));
[m,c] = size(elist);

resize(g,m);
rmxy(g);
clear_edges(g);

for i=1:m-1;
    a = elist(i,:);
    for j=i+1:m
        b = elist(j,:);
        if common_endpoint(a,b)
            add(g,i,j);
        end
        if (a(2) < b(1))
            break
        end
    end
end


if hasxy(h)
    xy = getxy(h);
    gxy = zeros(m,2);
    for i=1:m
        a = elist(i,1);
        b = elist(i,2);
        gxy(i,:) = (xy(a,:) + xy(b,:))/2;
    end
    embed(g,gxy);
end


function yn = common_endpoint(e1,e2)
yn = ( e1(1) == e2(1) ) | ( e1(1) == e2(2) ) | ( e1(2) == e2(1) ) | ...
    (e1(2) == e2(2));
    
