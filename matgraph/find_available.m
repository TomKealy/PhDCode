function idx = find_available
% idx = find_available
% Find the first available slot in the GRAPH_MAGIC machinery or 0 if no
% slot is available. 
% NOTE: This is used by the @graph constructor. No reason the user should
% call this directly.



global GRAPH_MAGIC


% see if all the slots are taken
if (sum(GRAPH_MAGIC.in_use) == GRAPH_MAGIC.ngraphs)
    idx = 0;
    return
end

% If we made it this far, there's an available slot. Now we just find one. 
idx = find(GRAPH_MAGIC.in_use == 0);
idx = idx(1);