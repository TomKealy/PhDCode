function [hlist, exists] = hamiltonian_cycle(g,h)
% hamiltonian_cycle(g) --- find a Hamiltonian cycle in g (if one exists)
% This can be called with two output arguments: 
%    [hlist, exists] = hamiltonian_cycle(g)
% the 2nd output argument is set to 0 if no such cycle exists.
% This can also be called hamiltonian_cycle(h,g). In this case, h is
% overwritten with a Hamiltonian cycle of g (if one exists).
%
% The important part of this function was written by Brian Towne.

if (nargin > 1) 
    [hlist, exists] = hamiltonian_cycle(h);
    if (exists)
        n = nv(h);
        copy(g,h);
        clear_edges(g);
        for k=1:n-1
            add(g,hlist(k),hlist(k+1));
        end
        add(g,hlist(1),hlist(n));
    end
    return;
end

n = nv(g);
used = false(n,1); % used(i)==true if i is in partially constructed cycle
path = zeros(n,1); % partially constructed cycle

% if there are no vertices in the graph, return an empty list, too
if nv(g) == 0
    hlist = [];
    exists = true;
    return
end


% We have to make sure there is exactly one component

c = components(g);
if np(c) > 1
    hlist = [];
    exists = false;
    return
end

% Since we are looking for a cycle, we can start w/ vertex 1 wolog:
used(1) = true;
path(1) = 1;

% Recursively try to extend path
[hlist,exists] = extend(g,1,n,used,path);
return



function [hlist, exists] = extend(g,i,n,used,path)

for j=neighbors(g,path(i))
    if used(j)==false
        path(i+1) = j;
        used(j) = true;
        if ((i+1==n) && (has(g,j,path(1))))
            hlist = path;
            exists = true;
            return
        else
            [hlist,exists] = extend(g,i+1,n,used,path);
            if exists
                return
            end
        end
        used(j) = false;
    end
end
hlist = [];
exists = false;
return

