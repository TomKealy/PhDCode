function q_init(n)
% q_init(n) --- initialize private queue with capacity n
% The queue supports the following modifying operations
% q_push: add an element to the end of the queue
% q_pop_front: pop an element from the front of the queue
% q_pop_back: pop an element from the end of the queue

global GRAPH_MAGIC


GRAPH_MAGIC.Q.array = zeros(n,1);
GRAPH_MAGIC.Q.first = 1;
GRAPH_MAGIC.Q.last = 0;