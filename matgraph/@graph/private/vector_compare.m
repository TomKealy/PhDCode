function sg = vector_compare(v,w)
% vector_compare(v,w) -- lexicographically compare row vectors
% v, w are row vectors of the same length
% return -1, 0, 1 if v is lexicographically less, equal, or greater than w

if all(v==w)
	sg = 0;
	return
end
d=v-w;
a = find(d);
sg = sign(d(a(1)));
