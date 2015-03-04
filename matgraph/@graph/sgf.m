function A = sgf(g,M)
% sgf --- simple graph format: a 2-column matrix representation
% A = sgf(g) --- make the simple graph format matrix of g
% sgf(g,M)   --- overwrite g with the graph specified in M
%
% The simple graph format is a 2-column matrix representation of a graph.
% This format optionally includes embedding information.
% If A is the SGF of a graph, the rows of A are as follows:
%
% * First row is [n m] where n is the number of vertices and m is the number
%   of edges.
% * Rows 2 through m+1 specify the edges of the graph.
% * Optionally, rows (m+1)+1 through (m+1)+n give the x,y-coordinates of
%   the vertices.

if nargin==1
    A = sgf_make(g);
    return
end
    
if ~sgf_check(M)
    error('Matrix is not a properly formed sgf matrix');
end
 

n = M(1,1);
m = M(1,2);

resize(g,0);

if n>set_large
    sparse(g);
else
    full(g);
end

resize(g,n);

edges = M(2:m+1,:);
add(g,edges);

% see if there's an embedding
[r,c] = size(M);
if (r>m+1)
    xy = M(m+2:end,:);
    embed(g,xy);
end



function A = sgf_make(g)
% helper function to convert a graph to a sgf matrix
n = nv(g);
m = ne(g);
if hasxy(g)
    A = zeros(n+m+1,2);
else
    A = zeros(m+1,2);
end

A(1,:) = [n,m];
A(2:m+1,:) = edges(g);
if hasxy(g)
    A(m+2:end,:) = getxy(g);
end

function ok = sgf_check(M)
% check a matrix to be sure it's a valid sgf

ok = false;

[r,c] = size(M);
if (c ~= 2) || (r==0)
    return
end

row1 = M(1,:);
n = row1(1);
m = row1(2);

if (n < 0) || (m < 0)
    return
end

if ~((r==m+1) || (r==m+n+1))
    return
end

edges = M(2:m+1,:);
% make sure they are all integers
if any(any(edges ~= round(edges)))
    return
end
% make sure they are all in {1,2,...,n}
if any(any(edges < 1)) || any(any(edges>n))
    return
end
% make sure there are no loops
if any(edges(:,1)==edges(:,2))
    return
end

% all tests passed
    
ok = true;
