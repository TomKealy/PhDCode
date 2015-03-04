function qvals = q_get
% q_get --- get the values held in the queue

global GRAPH_MAGIC

if q_size == 0
    qvals = [];
    return
end

i = GRAPH_MAGIC.Q.first;
j = GRAPH_MAGIC.Q.last;

qvals = GRAPH_MAGIC.Q.array(i:j);