%% Eulerian tours
% The |euler_trail| function finds Eulerian trails and tours in graphs.
%% Line graph of the Petersen graph
% The Petersen graph has no Eulerian trail or tour, but its line graph
% does. The Eulerian trail is reported as a list of edges.
g = graph;
petersen(g);
h = graph;
line_graph(h,g);  
distxy(h);
clf;ndraw(h);
euler_trail(h)
%% 4-cube minus an edge
% The four dimensional cube is a 4-regular graph. Deleting one edge from
% this graph leaves a connected graph with exactly two vertices of odd
% degree. We find an Euler trail in that graph.
k = graph;
cube(k,4);
distxy(k);
delete(k,1,2);
euler_trail(k)
clf;ndraw(k)
%% Release storage
free(g); free(h); free(k);