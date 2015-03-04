function p = color(g,algo)
% color(g,algo) --- color the graph g by a given algorithm
% At present, the only algorithm (and the default) is 'greedy'.


if nargin == 1
    algo = 'greedy';
end

if strcmp(algo,'greedy')
    p = greedy_color(g);
    return
end

error(['Algorithm "', algo, '" not implemented']);



function p = greedy_color(g)

n = nv(g);
d = deg(g);
[dd,vlist] = sort(-d);

colors = zeros(1,n);

% scan vertices by decreasing degree
for k=1:n
    v = vlist(k);
    vneig = neighbors(g,v);  % neighbors of v
    cneig = unique(colors(vneig));  % colors on neighbors
    if length(cneig) == 0
        mc = 0;
    else
        mc = cneig(length(cneig)); % max color seen on neighbors
    end
    avail = setdiff(1:mc,cneig); % colors available (if any)
    if isempty(avail)
        colors(v) = mc+1;
    else
        colors(v) = min(avail);
    end
end

mc = max(colors);
cpart = cell(1,mc);
for k=1:mc
    cpart{k} = find(colors==k);
end
p = partition(cpart);
    
    