function [elist, exists] = euler_trail(g)
% euler_trail(g) --- find an euler trail in g (if one exists)
% This can be called with two output arguments: 
%    [elist, exists] = euler_trail(g)
% the 2nd output argument is set to 0 if no trail exists.

n = nv(g);
d = deg(g);
d = mod(d,2);
if sum(d) > 2
    elist = []; % no trail possible, too many odd vertices
    exists = 0;
    return
end

% if there are no edges in the graph, return an empty list, too

if ne(g) == 0
    elist = [];
    exists = 1;
    return
end

% we have to make sure there is exactly one nontrivial component

c = components(g);
nc = np(c);
c = parts(c); % convert to cell array
counts = zeros(nc,1);
for k=1:nc
    counts(k) = length(c{k});
end

big = find(counts > 1);
if (length(big) > 1)
    elist = [];
    exists = 0;
    return
end


% if there are two odd vertices, we have to start on one of those.
% otherwise, we start at the first vertex in the big component.

if sum(d)==2
    odds = find(d);
    start = odds(1);
else
    start = c{big};
    start = start(1);
end


% we make an "erasable" copy of g 
h = graph;
copy(h,g);

% and now we start the tour

m = ne(g);
elist = zeros(m,2);
step = 0;
v = start;

while ne(h)
    step = step+1;
    next = neighbors(h,v);
    % if v has even degree (or only one neighbor) any neighbor will do
    if (mod(length(next),2) == 0) | (length(next) == 1)
        w = next(1);
    else
        for w=next
            % make sure vw is not a cut edge
            delete(h,v,w);
            if (has_path(h,v,w))
                break
            else
                add(h,v,w) % restore the edge
            end
        end
    end
    delete(h,v,w);
    elist(step,:) = [v,w];
    v = w;
end
free(h)