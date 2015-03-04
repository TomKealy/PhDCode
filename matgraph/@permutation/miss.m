function s = miss(p,dir)
% miss(p) -- maximum increasing subsequence
% Viewing the permutation p as a sequence (return from array(p)) find a
% longest increasing subsequence. May also be called with two arguments:
% miss(p,'up') -- find longest increasing subsequence
% miss(p,'down') -- find longest decreasing subsequence

if nargin < 2
    dir = 'up';
end

switch(lower(dir))
    case 'up'
        a = array(p);
        s = miss_array(a);
    case 'dn'
        a = array(p);
        a = fliplr(a);
        s = miss_array(a);
        s = fliplr(s);
    otherwise
        error('Second argument should be ''up'' or ''down''.');
end
end




function ss = miss_array(a)
n = length(a);
up = zeros(1,n);
up(n) = 1;

% create an array up whose i'th element gives the length of the longest
% increasing subsequence starting at position i.
for i=(n-1):-1:1
	up(i) = 1;
	for j=(i+1):n
		if (a(i) < a(j))
			if (up(i) <= up(j))
				up(i) = up(j)+1;
			end
		end
	end			
end

% now use the up array to find longest subseq
maxup = max(up);
si = find(up==maxup,1,'first');
for k=(maxup-1):-1:1
    tail = up(si(end)+1:end);
    i = find(tail==k,1,'first');
    si = [si, si(end)+i];
end
    
ss = a(si);
   
end
    
