function mycielski(g,h)
% mycielski(g,n) --- create n'th order Mycielski graph (n > 0)
% mycielski(g,h) --- apply Mycielski step to h to give g
% e.g., mycielski(g,3) overwrites G with the Grotztsch graph
% Thanks: Kim Tucker


% if 2nd argument is a graph, apply construction one time

if (isa(h,'graph'))
    B = matrix(h);
    n = length(B);

    A = [ B, B,  zeros(n,1);
          B,   zeros(n),  ones(n,1);
          zeros(1,n), ones(1,n), 0 ];
    fast_set_matrix(g,A);
    rmxy(g);
    return
end


% otherwise, 2nd argument is an integer

h = floor(h);

if (h<=0)
    complete(g,1);
    return
end

if (h == 1)
    complete(g,2);
    return
end

mycielski(g,h-1);
mycielski(g,g);