function a = array(p)
% array(p) --- return a single array showing the part numbers of elements.
% If p is the partition { {1,2,4} {3,5,6} {7,8,9,10} } then array(p) gives 
% the array [1 1 2 1 2 2 3 3 3 3].

n = nv(p);
a = zeros(1,n);

c = parts(p);
for k=1:length(c)
    a(c{k}) = k;
end