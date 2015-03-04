%% Checking graph isormorphism in Matgraph
% We show how to determine if two graphs are isomorphic and, if so, how to
% find the permutation that transforms one into the other.

%% Cycle graph
% Create two copies of a cycle on 10 vertices, but randomly renumber the
% second. 

g = graph;
h = graph;
cycle(g,10);
cycle(h,10);
renumber(h,random(permutation(nv(h))));
figure(1); clf; ndraw(g);
figure(2); clf; ndraw(h);

[yn,p] = iso(g,h)

renumber(g,p)
if g==h
	disp('Graphs are the same')
else
	disp('Graphs are different')
end

%% Large, vertex regular graph
% The bucky graph is a vertex transitive graph with 60 vertices. We make
% two copies (one that is randomly renumbered) and check they are the same.

bucky(g);
bucky(h);
renumber(h,random(permutation(nv(h))));

[yn,p] = iso(g,h)

renumber(g,p)
if g==h
	disp('Graphs are the same')
else
	disp('Graphs are different')
end

%% Standard random graph
% We generate a standard Erdos-Renyi random graph and a randomly renumbered
% copy, and then verify they are isomorphic.

random(g,100,0.5)
copy(h,g)
renumber(h,random(permutation(nv(h))));

[yn,p] = iso(g,h);

renumber(g,p)
if g==h
	disp('Graphs are the same')
else
	disp('Graphs are different')
end


%% Random regular graph
% We generate a random 3-regular graph and a randomly renumbered
% copy, and then verify they are isomorphic.

random_regular(g,100,3);
copy(h,g);
renumber(h,random(permutation(nv(h))));

[yn,p] = iso(g,h);

renumber(g,p)
if g==h
	disp('Graphs are the same')
else
	disp('Graphs are different')
end

%% Random tree
% We generate a random tree and a randomly renumbered
% copy, and then verify they are isomorphic.

random_tree(g,100);
copy(h,g);
renumber(h,random(permutation(nv(h))));

[yn,p] = iso(g,h);

renumber(g,p)
if g==h
	disp('Graphs are the same')
else
	disp('Graphs are different')
end

%% Two ways to make a grid
% We generate a grid graph in two different ways. The graphs are not equal,
% but they are isomorphic.

grid(g,6,3);
grid(h,3,6);
if g==h
	disp('Graphs are the same')
else
	disp('Graphs are different')
end

[yn,p] = iso(g,h)

renumber(g,p)
if g==h
	disp('Graphs are the same')
else
	disp('Graphs are different')
end

%% Mobius and non-Mobius ladders are not isomorphic
% The Mobius ladder on 20 vertices and the Cartesian product of a 10-cycle
% and an edge are both 3-regular graphs on 20 vertices. But they are not
% isomorphic.

mobius(g,20);
c10 = graph; cycle(c10,10);
k2 = graph; add(k2,1,2);
cartesian(h,c10,k2);

figure(1); clf; draw(g);
distxy(h);
figure(2); clf; draw(h);

[yn,p] = iso(g,h)

if yn
	disp('Graphs are isomorphic')
else
	disp('Graphs are not isomorphic')
end



%% Release storage

free(g); free(h); free(c10); free(k2);