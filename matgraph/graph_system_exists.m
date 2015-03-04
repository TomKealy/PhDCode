function tf = graph_system_exists
% graph_system_exists checks to see if the GRAPH_MAGIC global structure has
% been initialized.
tf = ~isempty(whos('global','GRAPH_MAGIC'));