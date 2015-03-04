function x = q_pop_front
% q_pop_front --- pop an element off the front of the queue (like a queue)
% if the queue is empty, an empty [] array is returned

global GRAPH_MAGIC


if (GRAPH_MAGIC.Q.last < GRAPH_MAGIC.Q.first)
    x = [];   % empty -- cannot return an element
    return
end

x = GRAPH_MAGIC.Q.array(GRAPH_MAGIC.Q.first);
GRAPH_MAGIC.Q.first = GRAPH_MAGIC.Q.first + 1;