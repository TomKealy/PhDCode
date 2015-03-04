function p = color(g,algo,max_time)
% color(g,algo) --- color the graph g by a given algorithm
% The algorithms are as follows:
%
% 'greedy':     the default, puts vertices in descending order of degree 
%               and runs sequential coloring
%
% 'rs':         random sequence, puts vertices in random order and 
%               runs sequential coloring
%
% 'repeat':     runs the 'rs' algorithm repeatedly until a fixed amount
%               of time passes (30 seconds); this can also be called as
%               color(g,'repeat',max_time)
%
% 'optimal':    always finds a minimum coloring (slow on large graphs)
%
% 'matlab':     use the COLOR method from the Optimization Toolbox
%               (this is experimental --- it appears to give the same
%               result as 'greedy')
%
% Algorithms rs, repeat, and optimal coded by Brian Towne.

if nargin <= 1
    algo = 'greedy';
end


switch lower(algo)
	case 'greedy'
		p = greedy_color(g);
	case 'rs'
		p = rs_color(g);
	case 'repeat'
		if nargin <=2
			max_time = 30;
		end
		p = repeat(g,max_time);
	case 'optimal'
		p = optimal_color(g);
	case 'matlab'
		p = matlab_color(g);
	otherwise
		error(['Algorithm "', algo, '" unimplemented']);
end
			





function p = greedy_color(g)

n = nv(g);
d = deg(g);
[dd,vlist] = sort(-d);
clear dd;
p = sequential_color(g,vlist,n);
return



function p = rs_color(g,max_colors)
% Random sequence

n = nv(g);
if nargin <= 1
    max_colors = n;
end
vlist = randperm(n);
p = sequential_color(g,vlist,max_colors);
return



function p = repeat(g,max_time)
% Starts w/ greedy coloring, then loops through random permutations 
% until it finds an optimal coloring.

tic;
% n = nv(g);
p = greedy_color(g);
max_colors = np(p)-1;
while (toc < max_time)
    temp = rs_color(g,max_colors);
    if ~ (np(temp)==0) % if rs_color returned w/o exhausting max_colors
        p = temp;
        max_colors = np(p)-1;
    end
end
return



function p = optimal_color(g)
% Finds a minimum coloring for g. Works by calling complete_extend.

n = nv(g);
p = greedy_color(g);
p = complete_extend(g,1,zeros(1,n),p);
return



function p = complete_extend(g,v,colors,best)
% Recursive function. Loops through every possible color for a given vertex.
% p is the current coloring. best is the best coloring found thus far.

n = nv(g);
vneig = neighbors(g,v);  % neighbors of v
cneig = unique(colors(vneig));  % colors on neighbors
avail = setdiff(1:min(v,np(best)-1),cneig); % colors available
for k=avail
    colors(v)=k;
    if v == n
        mc = max(colors);
        cpart = cell(1,mc);
        for k=1:mc
            cpart{k} = find(colors==k);
        end
        p = partition(cpart);
        if np(p) < np(best)
            best = p;
        end
    else
        best = complete_extend(g,v+1,colors,best);
    end
end
p = best;
return



function p = sequential_color(g,vlist,max_colors)
% Runs sequential coloring on vertices in order of vlist
% Exits if it uses more than max_colors. Setting max_colors = nv(g) 
% means this will never happen.

n = nv(g);
colors = zeros(1,n);

% scan vertices by order in vlist
for k=1:n
    v = vlist(k);
    vneig = neighbors(g,v);  % neighbors of v
    cneig = unique(colors(vneig));  % colors on neighbors
    if isempty(cneig) 
        mc = 0;
    else
        mc = cneig(length(cneig)); % max color seen on neighbors
    end
    avail = setdiff(1:mc,cneig); % colors available (if any)
    if isempty(avail)
        if mc+1 > max_colors % if max_colors is exhausted, exit
            p = partition(0);
            return
        else
            colors(v) = mc+1;
        end
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
return




function p = matlab_color(g)
% Uses the COLOR method from the optimization toolbox.
% That COLOR takes as input an interesection graph representation of g. The
% easiest way for us to do that is to use the incidence matrix, but perhaps
% if we used a smarter intersection graph representation, we could get
% better results from the Optimization Toolbox's COLOR.
M = incidence_matrix(g)';
i = color(M);
nc = max(i); 
parts = cell(nc,1);
for k=1:nc
	parts{k} = find(i==k);
end
p = partition(parts);
return


