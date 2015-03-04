%% Lower bound for the Ramsey number R(4)
% We show that the the Ramsey number R(4) is greater than 17 by presenting
% a graph with 17 vertices that has clique and independence numbers equal to
% three. The Matgraph |alpha| and |omega| functions require the
% Optimization Toolbox.
%% The Paley graph
% The graph we use is the Paley graph with n=17. The vertices of this graph
% correspond to the integers mod 17 with an edge between two vertices
% exactly when their difference is a quadratic residue mod 17. The |paley|
% function creates such a graph. 
g = graph;
paley(g,17);
%% Maximum clique
% We use |omega| to find the clique number and a maximum clique.
[w,S] = omega(g);
disp(['The clique number is ', int2str(w)])
notS = setdiff(1:17,S);
p = partition({S,notS});
clf;cdraw(g,p);
%% Maximum independent set
% We use |alpha| to find the independence number and a maximum independent
% set.
[a,S] = alpha(g);
disp(['The independence number is ', int2str(a)])
notS = setdiff(1:17,S);
p = partition({S,notS});
clf;cdraw(g,p);
%% Release storage
free(g)