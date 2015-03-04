function [yn,p] = iso(g,h)
% [yn,p] = iso(g,h) --- is g isomorphic to h?
% Given graphs g and h, the function attempts to determine if g and h are
% isomorphic by applying various heurisitcs. If it can determine the graphs
% are not isomorphic, yn is set to 0. If it can determine the graphs are
% isomorphic, yn is set to 1 and the isomophism (as a permutation object)
% is returned in p. Otherwise, if the heuristics are inconclusive, yn is
% set equal to 0.5.
%
% When a permutation is found, then renumber(g,p) will be equal to h. If
% the graphs are not isomorphic or if no isomorphism can be found, an empty
% permutation is returned.
%
% The heuristic approach here is not state of the art; actually, it's 
% really quite bad. If the graphs are isomorphic with trivial automorphism
% group, then this function is likely to find the isomorphism. 
%
% Users are strongly encouraged to help create a better routine.

%  basic tests 

% (0) check if the graphs are simply equal
if (g==h)
    p = permutation(nv(g));
    yn = 1;
    return
end


% (1) number of vertices & edges

if any(size(g) ~= size(h))
    yn = 0;
    return
end

% (2) degreee sequences
dg = sort(deg(g));
dh = sort(deg(h));
if any(dg ~= dh) 
    yn = 0;
    return
end

% (3a) same number of connected components 
gc = components(g);
hc = components(h);

if np(gc) ~= np(hc)
    yn = 0;
    return;
end

% (3b) components must all have the same sizes
gcp = parts(gc);
hcp = parts(hc);

gx = zeros(1,np(gc));
gy = gx;

for k=1:np(gc)
    gx(k) = length(gcp{k});
    gy(k) = length(hcp{k});
end
gx = sort(gx);
gy = sort(gy);
if any(gx ~= gy)
    yn = 0;
    return
end


% we now know that g and h have the same number of vertices and edges, have
% the same degree sequence, the same number of components, and the sizes of
% those components agree.

% we now work on copies of g and h

G = graph; copy(G,g);
H = graph; copy(H,h);

% if the graphs are not connected, work on their complements.

if (np(gc) > 1)
    complement(G);
    complement(H);
end

Lg = full(laplacian(G));
Lh = full(laplacian(H));


[Vg, Dg] = eig(Lg);
[Vh, Dh] = eig(Lh);

delta = sort(diag(Dg)) - sort(diag(Dh));
if (norm(delta) > 1000*eps)
    yn = 0;
    free(G);
    free(H);
    return;
end

Vg = Vg(:,2:end);
Vh = Vh(:,2:end);

[X,idx1] = sortrows(Vg);
[X,idx2] = sortrows(Vh);
p1 = permutation(idx1);
p2 = permutation(idx2);

yn = 0.5;

p = p2*inv(p1);

renumber(G,p)
if (G==H)
    free(G)
    free(H)
    yn = 1;
    return
end

Vh = -Vh;


[X,idx1] = sortrows(Vg);
[X,idx2] = sortrows(Vh);
p1 = permutation(idx1);
p2 = permutation(idx2);

p = p2*inv(p1);

copy(G,g)
renumber(G,p)
if (G==H)
    yn = 1;
else
    yn = 0.5;
    p = permutation(0);
end
free(G)
free(H)





