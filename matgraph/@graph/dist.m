function d = dist(g,v,w)
% dist(g,v,w) and dist(g,v) --- find distance(s) between vertices
% The form dist(g,v,w) finds the distance between vertices v and w.
% The form dist(g,v) returns a vector distance from v to all other vertices
% in the graph.
% The form dist(g) returns an n-by-n matrix whose ij entry is the distance
% between vertices i and j.
%
% The code for all pairs shortest distance was written by Michael Kleder
% and found on the Mathwork's MATLAB Central file exchange.


if nargin == 3
	dd = dist(g,v);
	d = dd(w);
    return
end

n = nv(g);


if nargin==1
    d = alltoall(g);
    return
end

d = onetoall(g,v);


% this subroutine was written by Michael Kleder and was found on The
% Mathwork's MATLAB Central file exchange. I think my new alltoall is
% faster, though.
% 
% function B = allspath(A)
% B=full(double(A));
% B(B==0)=Inf;
% C=ones(size(B));
% iter=0;
% while any(C(:))
%     C=B;
%     B=min(B,squeeze(min(repmat(B,[1 1 length(B)])+...
%         repmat(permute(B,[1 3 2]),[1 length(B) 1]),[],1)));
%     C=B-C;
% end
% B(logical(eye(length(B))))=0;
% return


% This is my replacement for allspath; I think this is faster for large
% sparse graphs.
function d = alltoall(g)
n = nv(g);
d = zeros(n,n);
for v=1:n
	d(:,v) = dist(g,v);
end
return


% This is a function that finds the shortest distance from v to all other
% vertices in a graph
function d = onetoall(g,v)
A = double(matrix(g));
n = nv(g);
vec = zeros(n,1);
vec(v) = 1;
d = inf*ones(1,n);
d(v) = 0;

count = 0;
while true
	count = count + 1;
	newv = (A*vec + vec)>0;
	morev = newv - vec;
	if (nnz(morev)==0)
		break
	end
	d(find(morev))=count;
	vec = newv;
end
return