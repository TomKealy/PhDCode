function yn = check(p)
% check(p) --- check that the datastructure holding p is a valid partition

s = full(sum(p.array,1));  % col sums should all be 1
t = full(sum(p.array,2));  % row sums should all be positive
yn = all(s == 1)  & all(t>0);
