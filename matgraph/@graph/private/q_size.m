function n = q_size
% determine the number of elements held in the queue

global GRAPH_MAGIC

n = GRAPH_MAGIC.Q.last - GRAPH_MAGIC.Q.first + 1;