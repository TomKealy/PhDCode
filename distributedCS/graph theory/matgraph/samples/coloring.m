%% Coloring graphs in Matgraph
% We illustrate how to find colorings of graphs and display the results.
%% Use the Buckyball graph
% We set |g| to be the Bucky graph.
g = graph;
bucky(g)
%% Greedy coloring
% First we find a greedy coloring of the graph and display the results.
% Notice that this coloring uses four colors.
c = color(g); % uses algorithm 'greedy' by default
disp(['This coloring uses ', int2str(np(c)),' colors']);
cdraw(g,c);
%% Optimal coloring
% An optimal coloring is slower (and may be intolerable on large graphs).
% However, for the Bucky graph, it uses fewer colors than the greedy
% algorithm.
tic; c = color(g,'optimal'); toc
disp(['This coloring uses ', int2str(np(c)),' colors']);
cdraw(g,c);
%% Release the graph
free(g);