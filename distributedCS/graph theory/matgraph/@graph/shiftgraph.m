function shiftgraph(g,k,t)
% shiftgraph(g,k,t) -- create a shiftgraph g based on t-tuples of k symbols

resize(g,0);
n = k^t;
if (n>1000) 
	sparse(g)
end
resize(g,n);

edgelist = [];


for v=1:n
	vec = num2tuple(v-1,k,t);
	label(g,v,int2str(vec))
	for j=0:k-1
		wec = shift(vec,j);
		w = tuple2num(wec,k)+1;
		if (v~=w)
			edgelist = [edgelist; [v,w]];
		end
	end
end
add(g,edgelist);


function tup = num2tuple(n,k,t)
% convert a number into a base-k t-tuple
tup = zeros(1,t);
for i=1:t
	tup(i) = mod(n,k);
	n = n - tup(i);
	n = n/k;
end


function n = tuple2num(tup, k)
% convert a base-k tuple into a number
t = length(tup);
pows = k.^(0:t-1);
n = dot(pows,tup);


function tt = shift(t,s)
% shift type symbol s into the tuple t and drop the last element
n = length(t);
tt = [s, t(1:n-1)];
