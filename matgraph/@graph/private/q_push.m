function q_push(x)
% q_push(x) --- add x to the back of the queue

global GRAPH_MAGIC

n = length(GRAPH_MAGIC.Q.array);

if GRAPH_MAGIC.Q.last == n
    error('GRAPH_MAGIC queue is full; cannot comply with q_push');
end

GRAPH_MAGIC.Q.last = GRAPH_MAGIC.Q.last + 1;
GRAPH_MAGIC.Q.array(GRAPH_MAGIC.Q.last) = x;