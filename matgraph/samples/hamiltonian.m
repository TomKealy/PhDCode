%% Finding and displaying Hamiltonian cycles in Matgraph
% We illustrate how to find a Hamiltonian cycle in a graph, convert that
% cycle into a subgraph, and then display the results.
%% Create a graph and find a Hamiltonian cycle
% For this illustration, we use the dodecahedron graph. 
%
% Note: The dodecahedron graph comes with a built-in embedding. We use this
% embedding later when we draw the graph and its Hamiltonian cycle.
g = graph;
dodecahedron(g); 
c = hamiltonian_cycle(g);
disp(c')
%% Convert the list of vertices into a cycle graph
n = nv(g);
h = graph(n); 
for k=1:n-1
	add(h,c(k),c(k+1))
end
add(h,c(1),c(n))
%% Copy the embedding of the original graph to the cycle
embed(h, getxy(g))
%% Create the drawing
% We first draw |g| with dotted lines and then overlay the result with |h|
% drawn with solid lines.
draw(g,':')
draw(h)
%% Release graphs 
free(g)
free(h)