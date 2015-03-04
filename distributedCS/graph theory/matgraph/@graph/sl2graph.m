function sl2graph(g,p)
% sl2graph(g,p) -- create an SL(2,p) graph
% g is the graph to be created
% p is a prime
% the vertices of g correspond to 2-by-2 matrices with determinant equal to
% 1 and entries modulo p
% Two vertices are adjacent if one can be obtained from the other by
% multiplication by either [1 1; 0 1] or [1 0; 1 1].
% In other words, g is a Cayley graph based on the group SL(2,p) with the
% two matrices as generators.

% Make sure p is prime

if ((p<1) || (length(factor(p)) > 1))
	error('p must be prime')
end

% Generate a table of all matrices [ a b ; c d ] in SL2(p). There are p^3
% of them that we store row-wise as [ a b c d ].

n = p^3-p;
vtable = zeros(n,4);
idx = 1;

% first the rows with a = 0
for c=1:p-1
	b = p-modinverse(c,p);
	for d=0:p-1
		row = [ 0, b, c, d ];
		vtable(idx,:) = row;
		idx = idx+1;
	end
end

% now the rows with a not equal to 0
for a=1:p-1
	for b=0:p-1
		for c=0:p-1
			aa = modinverse(a,p);
			d = aa*(b*c+1);
			d = mod(d,p);
			row = [ a, b, c, d ];
			vtable(idx,:) = row;
			idx = idx+1;
		end
	end
end


% DEBUG determinant check
for idx=1:size(vtable,1)
	row = vtable(idx,:);
	A = reshape(row,2,2);
	d = det(A);
	if (mod(d,p) ~= 1)
		disp([row,d])
	end
end

vtable = sortrows(vtable);

edgelist = [];
resize(g,0);
if (n>1000)
	sparse(g);
end
resize(g,n);

for v=1:n
	label(g,v,int2str(reshape(vtable(v,:),2,2)));
end

% Caley graph generators

S = [ 1 1 ; 0 1 ];
T = S';


for i = 1:n
	A = reshape(vtable(i,:),2,2);
	B = A*S;
	B = mod(B,p);
	j = findrow(B(:)',vtable);
	edgelist = [edgelist; [i,j]];
	
	B = A*T;
	B = mod(B,p);
	j = findrow(B(:)',vtable);
	edgelist = [edgelist; [i,j]];
end

add(g,edgelist);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = modinverse(a,p)
% computer inv(a) modulo p
[d,x,y] = gcd(a,p);
b = mod(x,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = findrow(row, A)
% find a row in a row sorted matrix, assume row is correct length
% return 0 if not found

n = size(A,1);
idx = 0;
if n==0
	return
end


mid = ceil(n/2);
midrow = A(mid,:);
if all(row == midrow)
	idx = mid;
	return
end

m = length(row); % presumably 4
for j=1:m
	if row(j) < midrow(j)
		idx = findrow(row,A(1:mid-1,:));
		return
	end
	if row(j) > midrow(j)
		idx = findrow(row,A(mid+1:end,:));
		if (idx>0)
			idx = idx+mid;
		end
		return
	end
end






