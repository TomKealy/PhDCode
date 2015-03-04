%% Splitting the vertex set of a graph
% We demonstrate how to split the vertex set of a graph into two parts
% using the |split| method. The goal is to separate the vertices of the
% graph into two natural clusters. The splitting is based on the Fiedler
% vector of the Laplacian matrix. 
%% Create a random tree
g = graph;
random_tree(g,25);
%% Find an embedding for the tree
% Note: This requires the Optimization Toolbox
distxy(g);
%% Partition the vertices
p = split(g)
%% Draw the result
cdraw(g,p);
%% Partition a grid graph
grid(g,3,11);
p = split(g);
clf; % erase the previous drawing
cdraw(g,p)
%% Release storage
free(g);