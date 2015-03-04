 function graph_init(capacity)
% graph_init(capacity) --- used to initialize the GRAPH_MAGIC data
% structure. 
% capacity gives the number of "slots" available for graphs.
% If capacity is omitted, a default capacity is used. See the
% code. 

DEFAULT_CAPACITY = 500;


if graph_system_exists
    disp('Graph system already initialized');
    disp('To start over, use graph_destroy');
    return;
end

global GRAPH_MAGIC;

% Fields in GRAPH_MAGIC:
%   ngraphs  number of graphs we can hold in this system
%   graphs   cell array to hold the graphs
%   in_use   flags to show which slots are available
%   large_size   threshold for sparse graphs
%   Q  a deque with subfields
%   Q.array   holds the queue values
%   Q.first   index of first element (or 0 if none)
%   Q.last    index of last element (or 0 if none)


if (nargin==0)
    capacity = DEFAULT_CAPACITY;
end


if (capacity < 1)
    capacity = DEFAULT_CAPACITY;
end

GRAPH_MAGIC.ngraphs = capacity;

GRAPH_MAGIC.graphs = cell(capacity,1);  % hold the graphs
GRAPH_MAGIC.in_use = zeros(capacity,1); % flag to show if slot is free

GRAPH_MAGIC.Q.array = [];
GRAPH_MAGIC.Q.first = 0;
GRAPH_MAGIC.Q.last = 0;


disp(['Graph system initialized. Number of slots = ', ...
    int2str(capacity),'.']);

set_large(1000);