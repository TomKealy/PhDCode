%% Creating a Cayley graph
%% Specifying the group by its generators
% We specify a finite group by giving a list of its generators. 
% The generators are a cell array of permutations of a common ground set.
% 
% In this example, we create a graph using two generators from the
% permutations of the set {1,2,3,4,5,6}.

p1 = permutation([1,3,6,4,5,2])
p2 = permutation([6,3,2,5,4,1])
perms = {p1,p2};

%% Create the graph
g = graph;
cayley(g,perms);

%% Draw the graph
distxy(g);
clf;
draw(g);

%% Another example
p1 = permutation([2 3 4 1])
p2 = swap(1,2)
cayley(g,{p1,p2},true) % verbose output
distxy(g);
clf;
draw(g);

%% Showing the labels on vertices
% The vertices of a Cayley graph are labeled by the group elements
% (permutations) they represent.
p1 = permutation([3 2 4 1]);
p2 = swap(1,2)*swap(3,4)
cayley(g,{p1,p2})
distxy(g);
clf;
ldraw(g)

%% The cube as a Cayley graph
p1 = permutation([2 1 3 4 5 6])
p2 = permutation([1 2 4 3 5 6])
p3 = permutation([1 2 3 4 6 5])
perms = {p1,p2,p3};
cayley(g,perms);
distxy(g);
clf;draw(g)

%% Free storage
free(g);
