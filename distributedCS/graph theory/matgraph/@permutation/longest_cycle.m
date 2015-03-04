function maxlen = longest_cycle(p)
% longest_cycle(p) -- return the length of the longest cycle in p
% where p is a permutation

clist = cycles(p);
maxlen = 0;
for k=1:length(clist)
	len = length(clist{k});
	if (len > maxlen)
		maxlen = len;
	end
end
