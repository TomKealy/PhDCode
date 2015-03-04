function x = q_pop_back
% q_pop_back --- pop the last element off the queue (like a stack)
% if the queue is empty, an empty array [] is returned

global GRAPH_MAGIC

if (GRAPH_MAGIC.Q.last < GRAPH_MAGIC.Q.first)
    x = [];   % empty -- cannot return an element
    return
end

x = GRAPH_MAGIC.Q.array(GRAPH_MAGIC.Q.last);
GRAPH_MAGIC.Q.last = GRAPH_MAGIC.Q.last - 1;
